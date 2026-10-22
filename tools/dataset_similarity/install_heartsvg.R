repos <- c(CRAN = "https://cloud.r-project.org")
install.packages(c(
  "spatstat",
  "Seurat",
  "plotly",
  "plot3D",
  "TTR",
  "stringr",
  "fpp2",
  "seasonal",
  "cluster",
  "poolr",
  "reshape2",
  "data.table",
  "plyr",
  "dplyr",
  "ggplot2",
  "MASS",
  "EnvStats",
  "BreakPoints",
  "gprofiler2",
  "devtools"
), repos = repos)
devtools::install_github("cz0316/HEARTSVG", force = TRUE)
