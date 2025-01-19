FROM rocker/rstudio:4.3.3

# Install necessary system libraries
RUN apt-get update \
&& apt-get install -y --no-install-recommends \
    clang \
    libgit2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libxml2-dev \
    libcairo2-dev \
    libsqlite-dev \
    libmariadbd-dev \
    libpq-dev \
    libcurl4-openssl-dev \
    libv8-dev \
    libssl-dev \
    cmake && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Set CRAN repository
RUN echo 'options(repos = c(CRAN = "https://packagemanager.rstudio.com/cran/2023-04-15"), download.file.method = "libcurl")' >> $(Rscript -e 'cat(R.home())')/etc/Rprofile.site

# Install R packages directly in Dockerfile
RUN Rscript -e "install.packages(c('testit', 'MASS', 'gsubfn', 'TOSTER', 'lavaan', 'tidyverse', 'meta', 'knitr', 'metafor', 'rJava', 'openxlsx', 'metaviz', 'readxl', 'gplots', 'pheatmap', 'RColorBrewer', 'ggplot2', 'dplyr', 'svglite', 'GGally', 'here', 'reshape2', 'ggpubr', 'psych', 'MBESS', 'xtable', 'kableExtra', 'gt', 'pagedown'))"

# Set working directory
WORKDIR /home/rstudio/analysis

# Copy scripts and data into the container
COPY . /home/rstudio/analysis
