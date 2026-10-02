FROM rocker/r-ver:4.4.2

ENV DEBIAN_FRONTEND=noninteractive

# Rocker versioned images pin both R and the CRAN snapshot. Install only system
# libraries with apt; install R packages into the Rocker R library with
# install2.r so packages cannot be pulled from Ubuntu's separate R runtime.
RUN apt-get update \
    && apt-get install -y --no-install-recommends libgmp-dev \
    && rm -rf /var/lib/apt/lists/*

# jsonlite/gmp/testthat cover the core computation and tests. ggplot2/viridis
# are required by the audit figures; future/future.apply enable the declared
# parallel path. If future is absent outside this image, the audit intentionally
# falls back to the statistically identical sequential implementation.
RUN install2.r --error --skipinstalled --ncpus -1 \
      jsonlite \
      gmp \
      testthat \
      ggplot2 \
      viridis \
      future \
      future.apply

RUN Rscript -e "required <- c('jsonlite','gmp','testthat','ggplot2','viridis','future','future.apply'); missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing)) stop('missing R packages: ', paste(missing, collapse=', '))"

WORKDIR /app
COPY . .

CMD ["Rscript", "tests/testthat.R"]
