## Load packages

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(rstanarm))
suppressPackageStartupMessages(library(loo))
suppressPackageStartupMessages(library(brms))
suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(bayesdfa))
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(mgcv))
suppressPackageStartupMessages(library(raster))
suppressPackageStartupMessages(library(colorspace))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(rjags))
suppressPackageStartupMessages(library(R2jags))
suppressPackageStartupMessages(library(codatools))
suppressPackageStartupMessages(library(bayesplot))
library(readxl)
library(ggsimple) ## see: https://github.com/michaelmalick/ggsimple
library(MARSS)
library(rnaturalearth)
library(rnaturalearthhires)
library(grid)
library(stringr)

dir.create("./figures/", showWarnings = FALSE)
dir.create("./outputs/", showWarnings = FALSE)

source("./functions.R")
source("./stan_utils.R")

M1 <- c("#00A0ED", "#DF1E7F", "#EAA301", "#6A4EE2", "#FF6000", "#000000")
M1d <- colorspace::darken(M1, amount = 0.2)

## Load output data
for(i in list.files(path = "./outputs/", pattern = "*.RData$")) {
    load(paste("./outputs/", i, sep = ""))
}

proj_wgs <- 4326  ## WGS84
