FROM rocker/r-ver:4.4.2

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends libgmp-dev \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages(c('jsonlite', 'gmp', 'testthat'), repos='https://cloud.r-project.org', Ncpus=2)"

WORKDIR /app
COPY . .

CMD ["Rscript", "tests/testthat.R"]
