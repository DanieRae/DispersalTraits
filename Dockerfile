# get the base image, the rocker/verse has R, RStudio and pandoc
FROM rocker/geospatial:4.3.2

# required
MAINTAINER Your Name <valentin.lucet@gmail.com>

COPY . /DispersalTraits

# go into the repo directory
RUN . /etc/environment \
  # Install linux depedendencies here
  # e.g. need this for ggforce::geom_sina
  && sudo apt-get update \
  && sudo apt-get install libudunits2-dev -y \
  # build this compendium package
  && R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))" \
  && R -e "remotes::install_github('rstudio/renv')" \
  # install pkgs we need
  && cd /DispersalTraits \
  && R -e "renv::restore()" \
  && R -e "remotes::install_github('quarto-dev/quarto-r')" \
  # render the manuscript into a docx, you'll need to edit this if you've
  # customised the location and name of your main qmd file
  && R -e "quarto::quarto_render('/DispersalTraits/analysis/paper/paper.qmd')"
