parse_arguments <- function() {
  raw <- commandArgs(trailingOnly = TRUE)
  out <- list()
  for (argument in raw) {
    if (!startsWith(argument, "--") || !grepl("=", argument, fixed = TRUE)) {
      stop("Arguments must have the form --name=value.", call. = FALSE)
    }
    pieces <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1L]]
    out[[pieces[1L]]] <- paste(pieces[-1L], collapse = "=")
  }
  out
}

arguments <- parse_arguments()
p <- as.integer(arguments$p)
regime <- arguments$regime
sources <- strsplit(arguments$sources, ",", fixed = TRUE)[[1L]]
output <- arguments$output
if (!p %in% c(12L, 16L) ||
    !regime %in% c("aligned", "perturbed", "crossed") ||
    length(sources) < 1L || any(!file.exists(sources)) ||
    is.null(output) || !nzchar(output)) {
  stop("Invalid chunk-merge arguments or missing source files.", call. = FALSE)
}

rows <- lapply(sources, read.csv, stringsAsFactors = FALSE)
all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
rows <- lapply(rows, function(row) {
  missing <- setdiff(all_names, names(row))
  for (name in missing) row[[name]] <- NA
  row[all_names]
})
rows <- do.call(rbind, rows)
rows <- rows[rows$p == p & rows$regime == regime, , drop = FALSE]

key <- with(rows, paste(replication, epsilon, interval_mode, sep = "/"))
duplicate_key <- unique(key[duplicated(key) | duplicated(key, fromLast = TRUE)])
if (length(duplicate_key) > 0L) {
  reproducible_fields <- setdiff(names(rows), "elapsed_seconds")
  for (duplicate in duplicate_key) {
    duplicated_rows <- rows[key == duplicate, reproducible_fields, drop = FALSE]
    signatures <- apply(duplicated_rows, 1L, paste, collapse = "|")
    if (length(unique(signatures)) != 1L) {
      stop(
        paste("Conflicting duplicate registered key:", duplicate),
        call. = FALSE
      )
    }
  }
  keep <- !duplicated(key, fromLast = TRUE)
  rows <- rows[keep, , drop = FALSE]
  key <- key[keep]
}
expected <- expand.grid(
  replication = seq_len(100L),
  epsilon = c(0.01, 0.05, 0.10, 0.20),
  interval_mode = c("exact", "stochastic"),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
expected_key <- with(expected, paste(replication, epsilon, interval_mode, sep = "/"))
missing_key <- setdiff(expected_key, key)
extra_key <- setdiff(key, expected_key)
if (length(missing_key) > 0L || length(extra_key) > 0L || nrow(rows) != 800L) {
  stop(
    sprintf(
      "Chunk merge is incomplete: rows=%d, missing=%d, extra=%d.",
      nrow(rows), length(missing_key), length(extra_key)
    ),
    call. = FALSE
  )
}
if (any(rows$status == "implementation-error")) {
  stop("Chunk merge contains an implementation error.", call. = FALSE)
}

epsilon_order <- match(rows$epsilon, c(0.01, 0.05, 0.10, 0.20))
interval_order <- match(rows$interval_mode, c("exact", "stochastic"))
rows <- rows[order(rows$replication, epsilon_order, interval_order), , drop = FALSE]
write.csv(rows, output, row.names = FALSE, na = "")
cat(sprintf(
  "Merged %d p=%d %s registered rows after validating %d duplicate keys.\n",
  nrow(rows), p, regime, length(duplicate_key)
))
