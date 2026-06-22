# helper-lib.R
# Source the R scanner library from whichever working directory testthat uses.
.find_lib <- function() {
  candidates <- c(
    file.path("..", "..", "R", "nist_scan_lib.R"),  # wd = tests/testthat
    file.path("R", "nist_scan_lib.R"),               # wd = repo root
    file.path("..", "R", "nist_scan_lib.R")
  )
  for (p in candidates) {
    if (file.exists(p)) return(normalizePath(p))
  }
  stop("could not locate R/nist_scan_lib.R")
}

LIB_PATH <- .find_lib()
REPO_ROOT_T <- normalizePath(file.path(dirname(LIB_PATH), ".."))
source(LIB_PATH)
