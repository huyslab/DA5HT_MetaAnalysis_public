# Following instructions from 
# https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Linux

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
	platform = "linux/arm64/v8"
}else{
	platform = args[1]
}

print("Running install_r.R")


# Install packages 
install.packages(c(
  "tidyverse", 
  "meta",
  "knitr",
  "metafor",
  "rJava",
  "openxlsx",
  "metaviz",
  "readxl",
  "gplots",
  "pheatmap",
  "RColorBrewer",
  "ggplot2",
  "dplyr",
  "svglite",
  "testit",
  "GGally",
  "here",
  "reshape2",
  "ggpubr",
  "psych",
  "MBESS",
  "xtable",
  "kableExtra",
  "gt",
  "pagedown"
))
