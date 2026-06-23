#!/usr/bin/env Rscript
# check_audit_dependencies.R
# ---------------------------------------------------------------------------
# Report which R packages required / suggested by the statistical-audit suite
# are installed. Prints ONE valid install.packages(...) command for whatever
# is missing. Returns a nonzero exit status if any REQUIRED package is absent.
# It never installs anything silently.
# ---------------------------------------------------------------------------

REQUIRED_PKGS <- c(
  "jsonlite", "gmp", "dplyr", "tidyr", "purrr", "readr", "stringr",
  "tibble", "ggplot2", "scales", "viridisLite"
)

OPTIONAL_PKGS <- c(
  # reporting / figures
  "patchwork", "ggridges", "ggrepel", "viridis", "gt", "ragg", "svglite",
  "quarto",
  # statistical tooling
  "boot", "broom", "rsample", "future", "future.apply", "progressr"
)

check_pkgs <- function(pkgs) {
  vapply(pkgs, function(p) requireNamespace(p, quietly = TRUE), logical(1))
}

main <- function() {
  req <- check_pkgs(REQUIRED_PKGS)
  opt <- check_pkgs(OPTIONAL_PKGS)

  cat("NIST statistical-audit dependency check\n")
  cat("=======================================\n\n")
  cat("Required packages:\n")
  for (p in REQUIRED_PKGS) cat(sprintf("  [%s] %s\n", if (req[p]) "ok " else "MISS", p))
  cat("\nOptional packages (suite degrades gracefully if missing):\n")
  for (p in OPTIONAL_PKGS) cat(sprintf("  [%s] %s\n", if (opt[p]) "ok " else " -- ", p))

  missing_req <- REQUIRED_PKGS[!req]
  missing_opt <- OPTIONAL_PKGS[!opt]
  # 'quarto' is an R package; the quarto CLI is checked separately below.
  to_install <- unique(c(missing_req, setdiff(missing_opt, "quarto")))

  cat("\n")
  if (length(to_install) > 0L) {
    cat("To install missing packages, run:\n\n")
    cat(sprintf('  install.packages(c(%s))\n\n',
                paste(sprintf('"%s"', to_install), collapse = ", ")))
  } else {
    cat("All required and optional R packages are installed.\n\n")
  }

  # quarto CLI (for report rendering) is independent of the R package
  quarto_cli <- nzchar(Sys.which("quarto"))
  cat(sprintf("Quarto CLI on PATH: %s\n", if (quarto_cli) "yes" else "no (HTML render will be skipped)"))

  if (length(missing_req) > 0L) {
    cat(sprintf("\nFAIL: %d required package(s) missing.\n", length(missing_req)))
    quit(status = 1L)
  }
  cat("\nOK: all required packages present.\n")
  quit(status = 0L)
}

.invoked_file <- sub("^--file=", "", commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))])
if (length(.invoked_file) > 0L && grepl("check_audit_dependencies\\.R$", .invoked_file)) {
  main()
}
