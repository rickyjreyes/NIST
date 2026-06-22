# Standard testthat runner.
#   Rscript tests/testthat.R
# or
#   Rscript -e 'testthat::test_dir("tests/testthat")'
library(testthat)
test_dir("tests/testthat", reporter = "check")
