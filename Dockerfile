FROM ubuntu:18.04
FROM rocker/tidyverse:4.2.1

MAINTAINER jon.brate@fhi.no

# Create script dir
RUN mkdir /home/scripts
WORKDIR /home/scripts

# Copy R script to the container
COPY create_samplesheet.R /home/scripts/create_samplesheet.R

# Run the script
CMD R -e "source('/home/scripts/create_samplesheet.R')"
