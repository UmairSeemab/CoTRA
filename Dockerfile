# syntax=docker/dockerfile:1

ARG R_VERSION=4.5.2
FROM rocker/r-ver:${R_VERSION}

LABEL org.opencontainers.image.title="CoTRA" \
      org.opencontainers.image.description="Comprehensive Toolbox for RNA Sequencing Data Analysis" \
      org.opencontainers.image.source="https://github.com/UmairSeemab/CoTRA" \
      org.opencontainers.image.licenses="GPL-3.0"

ENV DEBIAN_FRONTEND=noninteractive \
    R_LIBS_USER=/usr/local/lib/R/site-library \
    COTRA_DATA_DIR=/data \
    COTRA_RESULTS_DIR=/results \
    COTRA_HOST=0.0.0.0 \
    COTRA_PORT=3838

# System libraries required by CoTRA's CRAN/Bioconductor dependency stack.
RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    curl \
    ca-certificates \
    build-essential \
    cmake \
    pkg-config \
    pandoc \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libxslt1-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff-dev \
    libcairo2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libxt-dev \
    libglpk-dev \
    libgmp3-dev \
    libmpfr-dev \
    libgsl-dev \
    libhdf5-dev \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libmagick++-dev \
    libpoppler-cpp-dev \
    libv8-dev \
    libicu-dev \
    libbz2-dev \
    liblzma-dev \
    libpcre2-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/CoTRA
COPY . /opt/CoTRA

# Install the dependency stack first. The script exits non-zero if a required
# dependency is still missing, so dependency problems occur at image build time.
RUN R -q -e 'options(timeout=2000, repos=c(CRAN=Sys.getenv("CRAN", unset="https://cloud.r-project.org"))); install.packages(c("remotes", "BiocManager", "hdf5r"))' \
    && Rscript /opt/CoTRA/docker/install-dependencies.R \
    && R -q -e 'remotes::install_local("/opt/CoTRA", dependencies=FALSE, upgrade="never")' \
    && R -q -e 'd <- CoTRA::check_cotra_dependencies(quiet=TRUE); if (!isTRUE(d$ok)) stop(paste("Missing CoTRA dependencies:", paste(d$missing, collapse=", "))); if (!requireNamespace("hdf5r", quietly=TRUE)) stop("hdf5r is required for 10x HDF5 input")'

RUN mkdir -p /data /results /opt/cotra-runtime \
    && chmod -R 0777 /data /results /opt/cotra-runtime

COPY docker/start-cotra.sh /usr/local/bin/start-cotra
COPY docker/healthcheck.R /opt/cotra-runtime/healthcheck.R
RUN chmod +x /usr/local/bin/start-cotra

EXPOSE 3838
VOLUME ["/data", "/results"]

HEALTHCHECK --interval=30s --timeout=10s --start-period=120s --retries=5 \
  CMD Rscript /opt/cotra-runtime/healthcheck.R || exit 1

ENTRYPOINT ["/usr/local/bin/start-cotra"]
