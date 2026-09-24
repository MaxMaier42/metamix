# Pre-computes the vignettes.
#
# The vignette sources with live code are the *.Rmd.orig files in this
# directory. Fitting the models takes several minutes, which is too long for
# CRAN and for every R CMD check, so this script knits each source once into a
# static *.Rmd with the output and figures baked in. CRAN, R CMD check and
# pkgdown then render the static files without running Stan or RoBMA.
#
# Run from the package root after installing the current version of metamixr
# (devtools::install()); RoBMA (>= 4.0.0) and JAGS are needed for
# compareRoBMA. Takes about 10-15 minutes:
#   source("vignettes/precompute.R")

old_wd <- setwd("vignettes")

for (src in list.files(pattern = "[.]Rmd[.]orig$")) {
  out <- sub("[.]orig$", "", src)
  message("Knitting ", src, " -> ", out)
  knitr::knit(src, output = out, envir = new.env(), quiet = TRUE)
}

setwd(old_wd)
