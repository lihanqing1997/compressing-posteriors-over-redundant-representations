arguments <- commandArgs(trailingOnly = TRUE)
figure_directory <- if (length(arguments)) arguments[[1L]] else file.path("..", "figures")
output_directory <- file.path("sim", "output", "computational_study")
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

read_required <- function(name, columns) {
  path <- file.path(output_directory, name)
  if (!file.exists(path)) stop("Missing canonical result: ", path, call. = FALSE)
  rows <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  missing <- setdiff(columns, names(rows))
  if (length(missing)) {
    stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "), call. = FALSE)
  }
  rows
}

open_pdf <- function(name, width, height) {
  grDevices::pdf(
    file.path(figure_directory, name), width = width, height = height,
    family = "Helvetica", useDingbats = FALSE, onefile = TRUE
  )
}

regime_label <- c(
  aligned = "Aligned", perturbed = "Perturbed", crossed = "Crossed"
)

# Exact optimum and stochastic-mode certified bracket from the saved summary.
exact <- read_required(
  "tmlr_constrained_exact_summary.csv",
  c(
    "p", "regime", "epsilon", "interval_mode", "replications",
    "finite_incumbent_rate", "mean_exact_optimal_cost",
    "mean_certified_cost_lower", "mean_certified_cost_upper"
  )
)
exact <- exact[exact$interval_mode == "stochastic", , drop = FALSE]
open_pdf("certified_cost_frontier.pdf", 7.2, 4.4)
old <- par(
  mfrow = c(2, 3), mar = c(3.1, 3.4, 2.8, 0.8),
  oma = c(2.2, 2.3, 0.5, 0.3),
  mgp = c(2.0, 0.65, 0), tcl = -0.25, cex = 0.78
)
for (p_value in c(12L, 16L)) {
  for (regime in c("aligned", "perturbed", "crossed")) {
    cell <- exact[exact$p == p_value & exact$regime == regime, , drop = FALSE]
    cell <- cell[order(cell$epsilon), , drop = FALSE]
    if (nrow(cell) != 4L) stop("Incomplete exact-summary plotting cell.", call. = FALSE)
    values <- c(
      cell$mean_exact_optimal_cost, cell$mean_certified_cost_lower,
      cell$mean_certified_cost_upper
    )
    values <- values[is.finite(values) & values > 0]
    if (!length(values)) stop("No finite positive cost value in plotting cell.", call. = FALSE)
    plot(
      cell$epsilon, cell$mean_exact_optimal_cost, type = "o", pch = 16,
      col = "#202020", lwd = 1.5, log = "y",
      ylim = c(min(values), max(values) * 1.28),
      xaxt = "n", xlab = "", ylab = "",
      main = paste0("p = ", p_value, ", ", regime_label[[regime]]),
      cex.main = 1.02
    )
    axis(1, at = cell$epsilon, labels = format(cell$epsilon, trim = TRUE))
    lines(cell$epsilon, cell$mean_certified_cost_lower, type = "o", pch = 1,
          col = "#2B6CB0", lwd = 1.4)
    finite_upper <- is.finite(cell$mean_certified_cost_upper) &
      cell$mean_certified_cost_upper > 0
    lines(cell$epsilon[finite_upper], cell$mean_certified_cost_upper[finite_upper],
          type = "o", pch = 2, col = "#C53030", lwd = 1.4)
    label_x <- cell$epsilon
    label_x[[1L]] <- label_x[[1L]] + 0.008
    label_x[[length(label_x)]] <- label_x[[length(label_x)]] - 0.008
    text(
      label_x, rep(max(values) * 1.13, nrow(cell)),
      labels = paste0(round(100 * cell$finite_incumbent_rate), "%"),
      cex = 0.65, col = "#555555"
    )
    if (p_value == 12L && regime == "aligned") {
      legend(
        "bottomleft", bty = "n", cex = 0.66,
        legend = c("Exact optimum", "Certified lower", "Certified upper"),
        col = c("#202020", "#2B6CB0", "#C53030"),
        pch = c(16, 1, 2), lwd = 1.3
      )
    }
  }
}
mtext("Omitted-mass tolerance", side = 1, outer = TRUE, line = 0.7)
mtext("Mean declared cost (log scale)", side = 2, outer = TRUE, line = 0.7)
par(old)
dev.off()

cat("Wrote the constrained exact-study figure to", figure_directory, "\n")
