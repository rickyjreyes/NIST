FROM rocker/r-ver:4.4.2

ENV DEBIAN_FRONTEND=noninteractive

# Use distribution-built R packages instead of compiling the complete testthat
# dependency tree from CRAN during every image build. This removes transient
# source-build failures while retaining the repository-pinned R interpreter.
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libgmp-dev \
        r-cran-gmp \
        r-cran-jsonlite \
        r-cran-testthat \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "required <- c('jsonlite', 'gmp', 'testthat'); missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]; if (length(missing)) stop('missing R packages: ', paste(missing, collapse=', '))"

WORKDIR /app
COPY . .

CMD ["Rscript", "tests/testthat.R"]
