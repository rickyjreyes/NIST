# helper-audit.R
# Locate and source the statistical-audit utilities for testing. nist_scan_lib.R
# is already sourced by helper-lib.R, so audit_utils.R will not re-source it.
.find_R <- function(fname) {
  candidates <- c(
    file.path("..", "..", "R", fname),  # wd = tests/testthat
    file.path("R", fname),               # wd = repo root
    file.path("..", "R", fname)
  )
  for (p in candidates) if (file.exists(p)) return(normalizePath(p))
  stop(sprintf("could not locate R/%s", fname))
}
# Ensure the canonical library is loaded first (helper load order is not
# guaranteed), so audit_utils.R will not attempt to re-source it.
if (!exists("scan_k", mode = "function")) source(.find_R("nist_scan_lib.R"))
AUDIT_PATH <- .find_R("audit_utils.R")
source(AUDIT_PATH)
DATA_DIR_T <- normalizePath(file.path(dirname(AUDIT_PATH), "..", "data"))
