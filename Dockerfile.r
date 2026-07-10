FROM rocker/r-ver:4.4.2

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       build-essential \
       libcurl4-openssl-dev \
       libgmp-dev \
       libssl-dev \
       libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages(c('jsonlite', 'gmp', 'testthat'), repos='https://cloud.r-project.org', Ncpus=2); stopifnot(requireNamespace('jsonlite', quietly=TRUE), requireNamespace('gmp', quietly=TRUE), requireNamespace('testthat', quietly=TRUE))"

WORKDIR /app
COPY . .

CMD ["Rscript", "tests/testthat.R"]
