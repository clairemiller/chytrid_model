# Use the r-base image with R version 4.3.1
FROM r-base:4.3.1
  
# Install the remotes package (to install specific versions)
RUN R -e "install.packages('remotes')"

# Install specific versions of the required packages
RUN R -e "remotes::install_version('ggplot2', version = '3.5.1')"
RUN R -e "remotes::install_version('tidyr', version = '1.3.0')"
RUN R -e "remotes::install_version('dplyr', version = '1.1.1')"
RUN R -e "remotes::install_version('scales', version = '1.3.0')"
RUN R -e "remotes::install_version('GillespieSSA', version = '0.6.2')"
RUN R -e "remotes::install_version('assertthat', version = '0.2.1')"
RUN R -e "remotes::install_version('doParallel', version = '1.0.17')"
RUN R -e "remotes::install_version('tictoc', version = '1.2.1')"
RUN R -e "remotes::install_version('knitr', version = '1.42')"

# Set the working directory
WORKDIR /project

# Set the default command to run bash
CMD ["bash"]