FROM rocker/r-ver:4.4.1

LABEL maintainer="jon.brate@fhi.no"

# Install system dependencies for building R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# Install only the required R packages from CRAN
RUN Rscript -e "install.packages(c('dplyr','readr','stringr'), repos='https://cloud.r-project.org/')"

# Set working directory
WORKDIR /app

# Copy the R script to the container
COPY create_samplesheet.R /app/create_samplesheet.R

# Make the script executable
RUN chmod +x /app/create_samplesheet.R

# Default command to run the script
CMD ["Rscript", "/app/create_samplesheet.R"]
