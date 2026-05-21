FROM rocker/r-ver:4.5.1

ENV DEBIAN_FRONTEND=noninteractive
ENV RENV_PATHS_CACHE=/opt/renv/cache

WORKDIR /home

# 1. Install system dependencies FIRST
RUN apt-get update && apt-get install -y \
    cmake \
    default-jdk \
    git \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libglpk-dev \
    libgsl-dev \
    libicu-dev \
    libpng-dev \
    libssl-dev \
    libx11-dev \
    libxml2-dev \
    pandoc \
 && rm -rf /var/lib/apt/lists/*
 
# 2. Copy lockfile only
COPY renv.lock .

# 3. Install renv
RUN R -e "install.packages('renv', repos='https://cloud.r-project.org')"

# 4. Restore exact package versions
RUN R -e "renv::activate(); renv::restore(prompt = FALSE)"

# Optional: sanity check
RUN R -e "renv::status()"

CMD ["R"]