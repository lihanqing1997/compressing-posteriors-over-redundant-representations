fail <- function(message) stop(message, call. = FALSE)

if (!requireNamespace("yaml", quietly = TRUE)) {
  fail("The release verifier requires the yaml package.")
}

expect_true <- function(value, message) {
  if (!isTRUE(value)) fail(message)
}

expect_near <- function(value, target, tolerance, message) {
  if (length(value) != 1L || !is.finite(value) ||
      abs(value - target) > tolerance) {
    fail(sprintf("%s: got %.16g, expected %.16g", message, value, target))
  }
}

read_required <- function(path) {
  if (!file.exists(path)) fail(paste("Missing required file:", path))
  utils::read.csv(path, check.names = FALSE)
}

expect_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    fail(paste(label, "is missing columns:", paste(missing, collapse = ", ")))
  }
}

parse_sha256 <- function(output, path) {
  compact <- gsub("[[:space:]]", "", output)
  selected <- compact[grepl("^[[:xdigit:]]{64}$", compact)]
  if (length(selected) != 1L) fail(paste("Could not hash", path))
  tolower(selected[[1L]])
}

sha256_file <- function(path) {
  normalized <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (.Platform$OS.type == "windows") {
    command <- Sys.which("certutil")
    if (!nzchar(command)) fail("certutil is required for manifest checking.")
    return(parse_sha256(
      system2(command, c("-hashfile", shQuote(normalized), "SHA256"),
              stdout = TRUE, stderr = TRUE),
      normalized
    ))
  }
  command <- Sys.which("sha256sum")
  arguments <- shQuote(normalized)
  if (!nzchar(command)) {
    command <- Sys.which("shasum")
    arguments <- c("-a", "256", shQuote(normalized))
  }
  if (!nzchar(command)) fail("sha256sum or shasum is required.")
  output <- system2(command, arguments, stdout = TRUE, stderr = TRUE)
  token <- strsplit(output[[1L]], "[[:space:]]+")[[1L]][[1L]]
  parse_sha256(token, normalized)
}

pdf_printable_text <- function(bytes) {
  values <- as.integer(bytes)
  values[values < 32L | values > 126L] <- 32L
  rawToChar(as.raw(values))
}

pdf_inflated_text <- function(bytes, printable) {
  stream_starts <- gregexpr("stream", printable, fixed = TRUE)[[1L]]
  stream_ends <- gregexpr("endstream", printable, fixed = TRUE)[[1L]]
  if (identical(stream_starts, -1L) || identical(stream_ends, -1L)) return("")
  output <- character()
  for (marker in stream_starts) {
    following <- stream_ends[stream_ends > marker]
    if (!length(following)) next
    first <- marker + nchar("stream", type = "bytes")
    last <- following[[1L]] - 1L
    while (first <= last && as.integer(bytes[[first]]) %in% c(10L, 13L)) {
      first <- first + 1L
    }
    while (last >= first && as.integer(bytes[[last]]) %in% c(10L, 13L)) {
      last <- last - 1L
    }
    if (first > last) next
    inflated <- suppressWarnings(try(
      memDecompress(bytes[first:last], type = "gzip"),
      silent = TRUE
    ))
    if (!inherits(inflated, "try-error")) {
      output <- c(output, pdf_printable_text(inflated))
    }
  }
  paste(output, collapse = "\n")
}

inspect_pdf <- function(path) {
  size <- file.info(path)$size
  expect_true(is.finite(size) && size >= 8, paste("PDF is empty:", path))
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  bytes <- readBin(connection, what = "raw", n = size)
  printable <- pdf_printable_text(bytes)
  metadata_text <- paste(
    printable,
    pdf_inflated_text(bytes, printable),
    sep = "\n"
  )
  header <- substr(printable, 1L, 16L)
  version_match <- regexec(
    "^%PDF-([0-9]+)\\.([0-9]+)", header, perl = TRUE
  )
  version_parts <- regmatches(header, version_match)[[1L]]
  expect_true(length(version_parts) == 3L, paste("Invalid PDF header:", path))
  major <- as.integer(version_parts[[2L]])
  minor <- as.integer(version_parts[[3L]])
  expect_true(
    major < 1L || (major == 1L && minor <= 5L),
    paste("PDF version exceeds 1.5:", path)
  )
  expect_true(
    !grepl("/CreationDate", metadata_text, fixed = TRUE),
    paste("PDF contains CreationDate metadata:", path)
  )
  expect_true(
    !grepl("/ModDate", metadata_text, fixed = TRUE),
    paste("PDF contains ModDate metadata:", path)
  )
  author_tokens <- regmatches(
    metadata_text,
    gregexpr(
      "/Author[[:space:]]*(\\([^)]*\\)|<[^>]*>|[^/<>[:cntrl:]]+)",
      metadata_text,
      perl = TRUE
    )
  )[[1L]]
  if (length(author_tokens) && !identical(author_tokens, "")) {
    author_values <- sub(
      "^/Author[[:space:]]*", "", author_tokens, perl = TRUE
    )
    empty_author <- grepl(
      "^(\\([[:space:]]*\\)|<[[:space:]]*>)$", author_values, perl = TRUE
    )
    expect_true(
      all(empty_author),
      paste("PDF contains nonempty Author metadata:", path)
    )
  }
  geographic_or_user_path <- paste(
    c(
      "[A-Za-z]+[[:space:]]+(Standard|Daylight)[[:space:]]+Time",
      "(America|Europe|Asia|Africa|Australia|Pacific|Atlantic)/",
      "[A-Za-z]:[/\\\\]Users[/\\\\]",
      "(^|[/\\\\])Users[/\\\\]"
    ),
    collapse = "|"
  )
  expect_true(
    !grepl(
      geographic_or_user_path, metadata_text, ignore.case = TRUE, perl = TRUE
    ),
    paste("PDF contains a geographic timezone or user-path token:", path)
  )

  pdfinfo <- Sys.which("pdfinfo")
  if (nzchar(pdfinfo)) {
    output <- system2(
      pdfinfo,
      shQuote(normalizePath(path, winslash = "/", mustWork = TRUE)),
      stdout = TRUE,
      stderr = TRUE
    )
    status <- attr(output, "status")
    expect_true(is.null(status) || status == 0L,
                paste("pdfinfo could not inspect:", path))
    author_line <- grep("^Author:", output, value = TRUE)
    if (length(author_line)) {
      expect_true(
        all(!nzchar(trimws(sub("^Author:", "", author_line)))),
        paste("PDF contains nonempty Author metadata:", path)
      )
    }
    expect_true(
      !any(grepl("^(CreationDate|ModDate):", output)),
      paste("PDF exposes creation/modification date metadata:", path)
    )
    expect_true(
      !any(grepl(
        geographic_or_user_path, output, ignore.case = TRUE, perl = TRUE
      )),
      paste("PDF metadata exposes a geographic timezone or user path:", path)
    )
  }
  invisible(TRUE)
}

manifest <- read_required("MANIFEST_SHA256.csv")
expect_columns(manifest, c("path", "bytes", "sha256"), "manifest")
expect_true(nrow(manifest) > 0L, "Manifest is empty.")
manifest_paths <- as.character(manifest$path)
expected_manifest_paths <- c(
  "REPRODUCTION.txt",
  "sim/config/constrained_exact_protocol.yml",
  "sim/src/region_restricted_smc.R",
  "sim/src/exact_gaussian_helpers.R",
  "sim/exact_bma.R",
  "sim/run_constrained_exact_study.R",
  "sim/merge_constrained_exact_chunks.R",
  "sim/audit_constrained_exact_production.R",
  "sim/summarize_study_results.R",
  "sim/make_study_figures.R",
  "sim/validate_study_protocol.R",
  "sim/verify_release.R",
  "sim/run_region_restricted_smc_unit_tests.R",
  "sim/run_constrained_region_search_unit_tests.R",
  "sim/run_factorized_region_mass_unit_tests.R",
  unlist(lapply(c(12L, 16L), function(p) {
    paste0(
      "sim/output/computational_study/constrained_exact_p", p, "_",
      c("aligned", "perturbed", "crossed"), ".csv"
    )
  }), use.names = FALSE),
  "sim/output/computational_study/constrained_exact_audit.csv",
  "sim/output/computational_study/constrained_exact_summary.csv"
)
expect_true(
  identical(manifest_paths, expected_manifest_paths),
  "Manifest does not exactly match the 23 released files."
)
expect_true(
  !anyNA(manifest_paths) && all(nzchar(manifest_paths)),
  "Manifest contains an empty path."
)
expect_true(
  identical(manifest_paths, gsub("\\\\", "/", manifest_paths)),
  "Manifest paths must use forward slashes."
)
expect_true(
  !any(grepl("^(?:[A-Za-z]:|/)", manifest_paths, perl = TRUE)),
  "Manifest contains an absolute path."
)
unsafe_component <- vapply(
  strsplit(manifest_paths, "/", fixed = TRUE),
  function(components) any(!nzchar(components) | components %in% c(".", "..")),
  logical(1)
)
expect_true(!any(unsafe_component), "Manifest contains an unsafe path component.")
expect_true(
  !any(manifest_paths == "MANIFEST_SHA256.csv"),
  "The manifest must not list itself."
)
duplicate_key <- if (.Platform$OS.type == "windows") {
  tolower(manifest_paths)
} else {
  manifest_paths
}
expect_true(!anyDuplicated(duplicate_key), "Manifest contains duplicate paths.")
expect_true(
  is.numeric(manifest$bytes) && all(is.finite(manifest$bytes)) &&
    all(manifest$bytes >= 0) && all(manifest$bytes == floor(manifest$bytes)),
  "Manifest byte counts must be nonnegative whole numbers."
)

release_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
comparison_root <- if (.Platform$OS.type == "windows") {
  tolower(release_root)
} else {
  release_root
}
for (index in seq_len(nrow(manifest))) {
  path <- manifest_paths[[index]]
  expect_true(file.exists(path), paste("Manifest path is missing:", path))
  expect_true(!isTRUE(file.info(path)$isdir), paste("Manifest path is not a file:", path))
  link_target <- Sys.readlink(path)
  expect_true(
    is.na(link_target) || !nzchar(link_target),
    paste("Manifest path is a symbolic link:", path)
  )
  resolved <- normalizePath(path, winslash = "/", mustWork = TRUE)
  comparison_resolved <- if (.Platform$OS.type == "windows") {
    tolower(resolved)
  } else {
    resolved
  }
  expect_true(
    startsWith(comparison_resolved, paste0(comparison_root, "/")),
    paste("Manifest path escapes the release root:", path)
  )
  expect_true(file.info(path)$size == manifest$bytes[[index]],
              paste("Byte count differs:", path))
  expect_true(sha256_file(path) == tolower(manifest$sha256[[index]]),
              paste("SHA-256 differs:", path))
}

filesystem_entries <- list.files(
  ".", recursive = TRUE, all.files = TRUE, full.names = TRUE,
  include.dirs = TRUE, no.. = TRUE
)
link_targets <- Sys.readlink(filesystem_entries)
expect_true(
  !any(!is.na(link_targets) & nzchar(link_targets)),
  "Release tree contains a symbolic link."
)
actual_paths <- list.files(
  ".", recursive = TRUE, all.files = TRUE, full.names = FALSE,
  include.dirs = FALSE, no.. = TRUE
)
actual_paths <- gsub("\\\\", "/", actual_paths)
actual_paths <- actual_paths[actual_paths != "MANIFEST_SHA256.csv"]
missing_from_tree <- setdiff(manifest_paths, actual_paths)
unmanifested <- setdiff(actual_paths, manifest_paths)
expect_true(
  !length(missing_from_tree),
  paste("Manifest entries are absent from the release tree:",
        paste(missing_from_tree, collapse = ", "))
)
expect_true(
  !length(unmanifested),
  paste("Release tree contains unmanifested files:",
        paste(unmanifested, collapse = ", "))
)

pdf_paths <- manifest_paths[grepl("\\.pdf$", manifest_paths, ignore.case = TRUE)]
for (path in pdf_paths) inspect_pdf(path)

text_paths <- list.files(
  ".", recursive = TRUE, full.names = TRUE,
  pattern = "\\.(R|csv|md|yml|yaml|txt|tex)$", ignore.case = TRUE
)
forbidden <- paste(
  c(
    paste0("[A-Za-z]:[/", "\\\\", "]Users[/", "\\\\", "]"),
    "/home/[[:alnum:]_.-]+/",
    "/Users/[[:alnum:]_.-]+/",
    paste0("github", "\\.com/"),
    paste0("[[:alnum:]._%+-]+", intToUtf8(64),
           "[[:alnum:].-]+\\.[A-Za-z]{2,}"),
    paste0("University", "[[:space:]]+of")
  ),
  collapse = "|"
)
for (path in text_paths) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (any(grepl(forbidden, lines, ignore.case = TRUE, perl = TRUE))) {
    fail(paste("Potential anonymity leak:", path))
  }
}

raw_directory <- file.path("sim", "output", "computational_study")
protocol <- yaml::read_yaml(file.path("sim", "config", "constrained_exact_protocol.yml"))

dimensions <- as.integer(unlist(
  protocol$exact_constrained_study$dimensions, use.names = FALSE
))
regimes <- names(protocol$exact_constrained_study$data_regimes)
epsilon_grid <- as.numeric(unlist(
  protocol$global_design$epsilon, use.names = FALSE
))
interval_modes <- as.character(unlist(
  protocol$exact_constrained_study$interval_modes, use.names = FALSE
))
replications <- seq_len(as.integer(
  protocol$global_design$production_replications_per_cell
))
protocol_id <- as.character(protocol$protocol_id)
score_budget <- as.integer(
  protocol$global_design$production_settings$score_budget
)
seed_base <- as.integer(protocol$seed_registry$production_base)
required_exact_columns <- c(
  "protocol_id", "run_mode", "registered_production", "p", "regime",
  "replication", "epsilon", "delta", "interval_mode", "score_oracle",
  "data_seed", "algorithm_seed", "score_budget", "status",
  "simultaneous_interval_coverage", "exact_optimal_cost",
  "certified_cost_lower", "certified_cost_upper", "certified_cost_gap",
  "selected_cost", "exact_cost_optimal", "true_omitted_mass",
  "omitted_mass_bound", "omitted_mass_covered", "epsilon_feasible",
  "cost_sandwich_covered", "budget_overshoot"
)

key_number <- function(value) format(value, digits = 15L, trim = TRUE)
all_exact <- list()
index <- 0L
for (p_index in seq_along(dimensions)) {
  p <- dimensions[[p_index]]
  for (regime_index in seq_along(regimes)) {
    regime <- regimes[[regime_index]]
    path <- file.path(
      raw_directory,
      paste0("constrained_exact_p", p, "_", regime, ".csv")
    )
    rows <- read_required(path)
    expect_columns(rows, required_exact_columns, basename(path))
    expect_true(nrow(rows) == 800L, paste(path, "does not contain 800 rows."))
    expect_true(
      all(rows$p == p) && all(rows$regime == regime),
      paste(path, "contains a mismatched dimension or regime.")
    )
    expect_true(
      all(rows$protocol_id == protocol_id) &&
        all(rows$run_mode == "production") &&
        all(rows$registered_production),
      paste(path, "contains an unregistered or mismatched protocol row.")
    )
    expected_grid <- expand.grid(
      replication = replications,
      epsilon = epsilon_grid,
      interval_mode = interval_modes,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    actual_keys <- with(
      rows,
      paste(replication, key_number(epsilon), interval_mode, sep = "\r")
    )
    expected_keys <- with(
      expected_grid,
      paste(replication, key_number(epsilon), interval_mode, sep = "\r")
    )
    expect_true(
      identical(sort(actual_keys), sort(expected_keys)),
      paste(path, "does not contain the complete registered cell grid.")
    )
    epsilon_index <- match(key_number(rows$epsilon), key_number(epsilon_grid))
    method_index <- match(rows$interval_mode, interval_modes)
    expected_data_seed <- seed_base + 100000L * p_index +
      10000L * regime_index + 100L * rows$replication
    expected_algorithm_seed <- expected_data_seed +
      10L * epsilon_index + method_index
    expect_true(
      !anyNA(epsilon_index) && !anyNA(method_index) &&
        all(rows$data_seed == expected_data_seed) &&
        all(rows$algorithm_seed == expected_algorithm_seed),
      paste(path, "does not follow the registered seed formulas.")
    )
    expect_true(
      all(rows$score_budget == score_budget) &&
        all(rows$delta == as.numeric(protocol$global_design$delta)),
      paste(path, "contains a changed score budget or failure allocation.")
    )
    index <- index + 1L
    all_exact[[index]] <- rows
  }
}

exact <- do.call(rbind, all_exact)
expect_true(nrow(exact) == 4800L, "Constrained exact-study row count is not 4,800.")
expect_true(
  all(exact$score_oracle == "precomputed_exact_table"),
  "An exact-study row does not use the frozen precomputed exact-table score oracle."
)
exact_key <- with(
  exact,
  paste(p, regime, replication, key_number(epsilon), interval_mode, sep = "\r")
)
expect_true(!anyDuplicated(exact_key), "Constrained exact-study keys are duplicated.")
expect_true(!anyDuplicated(exact$algorithm_seed), "Registered algorithm seeds collide.")
expect_true(
  all(exact$status %in% as.character(unlist(
    protocol$global_design$status_levels, use.names = FALSE
  ))),
  "An exact-study row has an unknown registered status."
)
expect_true(
  sum(exact$interval_mode == "exact") == 2400L &&
    sum(exact$interval_mode == "stochastic") == 2400L,
  "The interval-mode row counts are not balanced."
)
expect_true(
  all(exact$cost_sandwich_covered) && all(exact$omitted_mass_covered),
  "A recorded cost sandwich or omitted-mass bound fails its exact check."
)
expect_true(
  !any(exact$status == "implementation-error") &&
    !any(exact$budget_overshoot > 0),
  "An implementation error or score-budget overshoot remains."
)

exact_mode <- exact$interval_mode == "exact"
stochastic_mode <- exact$interval_mode == "stochastic"
expect_true(
  all(exact$status[exact_mode] == "bounds-met") &&
    all(exact$exact_cost_optimal[exact_mode]) &&
    all(exact$certified_cost_gap[exact_mode] == 0),
  "An exact-mode row does not recover and certify the exact optimum."
)
expect_true(
  sum(exact$simultaneous_interval_coverage[stochastic_mode]) == 2394L,
  "The recorded stochastic-mode simultaneous-coverage count is not 2,394."
)
aligned_stochastic <- stochastic_mode & exact$regime == "aligned"
perturbed_stochastic <- stochastic_mode & exact$regime == "perturbed"
crossed_stochastic <- stochastic_mode & exact$regime == "crossed"
expect_true(
  sum(exact$status[aligned_stochastic] == "bounds-met") == 800L,
  "The aligned stochastic-mode cost-gap count is not 800."
)
expect_true(
  sum(is.finite(exact$certified_cost_upper[perturbed_stochastic])) == 798L &&
    sum(exact$status[perturbed_stochastic] == "bounds-met") == 1L &&
    sum(exact$status[perturbed_stochastic] == "budget-exhausted") == 795L,
  "The perturbed stochastic-mode status counts changed."
)
expect_true(
  sum(is.finite(exact$certified_cost_upper[crossed_stochastic])) == 622L &&
    sum(exact$status[crossed_stochastic] == "bounds-met") == 0L,
  "The crossed stochastic-mode status counts changed."
)

recompute_derived <- function(script) {
  output <- tempfile(fileext = ".csv")
  on.exit(unlink(output), add = TRUE)
  executable <- file.path(
    R.home("bin"),
    if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  )
  command_output <- system2(
    executable,
    c(shQuote(script), shQuote(output)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(command_output, "status")
  expect_true(
    is.null(status) || status == 0L,
    paste("Could not recompute derived artifact with", script)
  )
  read_required(output)
}

compare_derived <- function(saved, recomputed, keys, label) {
  expect_true(
    identical(names(saved), names(recomputed)),
    paste(label, "columns differ from their recomputation.")
  )
  order_frame <- function(data) {
    ordering <- do.call(order, data[keys])
    rownames(data) <- NULL
    data[ordering, , drop = FALSE]
  }
  saved <- order_frame(saved)
  recomputed <- order_frame(recomputed)
  expect_true(
    isTRUE(all.equal(
      saved, recomputed, tolerance = 1e-12, check.attributes = FALSE
    )),
    paste(label, "does not match recomputation from the production rows.")
  )
}

summary_path <- file.path(
  raw_directory, "constrained_exact_summary.csv"
)
summary <- read_required(summary_path)
expect_true(
  nrow(summary) == 48L && all(summary$replications == 100L),
  "The exact-study summary does not have 48 complete cells."
)
summary_recomputed <- recompute_derived(
  file.path("sim", "summarize_study_results.R")
)
compare_derived(
  summary, summary_recomputed,
  c("p", "regime", "epsilon", "interval_mode"),
  "Exact-study summary"
)

audit_path <- file.path(raw_directory, "constrained_exact_audit.csv")
audit <- read_required(audit_path)
expected_audit_checks <- c(
  "expected_rows",
  "deterministic_exactness_and_sandwich",
  "stochastic_simultaneous_coverage_rate",
  "aligned_finite_incumbent_rate_epsilon_0.10",
  "perturbed_finite_incumbent_rate_epsilon_0.10",
  "budget_overshoots",
  "implementation_errors"
)
expect_true(
  identical(as.character(audit$check), expected_audit_checks) &&
    nrow(audit) == 7L && all(audit$passed),
  "The seven-check protocol audit is absent, reordered, or did not pass."
)
audit_recomputed <- recompute_derived(
  file.path("sim", "audit_constrained_exact_production.R")
)
compare_derived(audit, audit_recomputed, "check", "Protocol audit")

message(sprintf(
  paste(
    "Release verification passed: %d hashed files;",
    "%d registered result rows and all derived fields checked."
  ),
  nrow(manifest),
  nrow(exact)
))
