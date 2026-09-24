## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a new submission.
* "unable to verify current time" is a local network/clock quirk of the
  machine the check was run on.

## Notes for the reviewers

* The package contains Stan models compiled at install time (rstantools
  package structure); installation takes several minutes.
* The vignettes are pre-computed (static R Markdown with the output included),
  because fitting the models takes minutes to hours. Their executable sources
  (`vignettes/*.Rmd.orig`) are excluded from the build. Consequently the
  RoBMA package and JAGS, which one vignette compares against, are not
  dependencies.
* Tests that run MCMC for more than a few seconds are skipped on CRAN
  (`skip_on_cran()`); the tests that run on CRAN take under a minute.
