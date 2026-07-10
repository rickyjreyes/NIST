FROM rocker/r-ver:4.4.2

ENV DEBIAN_FRONTEND=noninteractive

# Rocker versioned images pin both R and the CRAN snapshot. Install only system
# libraries with apt; install R packages into the Rocker R library with
# install2.r so packages cannot be pulled from Ubuntu's separate R runtime.
RUN apt-get update \
    && apt-get install -y --no-install-recommends libgmp-dev \
    && rm -rf /var/lib/apt/lists/*

RUN install2.r --error --skipinstalled --ncpus -1 \
      jsonlite \
      gmp \
      testthat

RUN Rscript -e "required <- c('jsonlite', 'gmp', 'testthat'); missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing)) stop('missing R packages: ', paste(missing, collapse=', '))"

WORKDIR /app
COPY . .

CMD ["Rscript", "tests/testthat.R"]
