# Prepares the `mertens_nudge` data set shipped in data/mertens_nudge.rda.
#
# Source: supplementary data of Mertens, Herberz, Hahnel & Brosch (2022),
# "The effectiveness of nudging: A meta-analysis of choice architecture
# interventions across behavioral domains", PNAS 119(1), e2107346118,
# doi:10.1073/pnas.2107346118, as used in the re-analysis by Maier et al.
# (2022), PNAS 119(31), e2200300119, doi:10.1073/pnas.2200300119.
#
# The object was originally read with readr::read_csv(), which attaches a
# column specification and a `problems` external pointer to the tibble. This
# script converts it to a plain data frame without those attributes and
# re-saves it. Run from the package root.

load("data/mertens_nudge.rda")

mertens_nudge <- as.data.frame(mertens_nudge)
attr(mertens_nudge, "spec") <- NULL
attr(mertens_nudge, "problems") <- NULL

stopifnot(identical(class(mertens_nudge), "data.frame"),
          nrow(mertens_nudge) == 447L, ncol(mertens_nudge) == 24L)

save(mertens_nudge, file = "data/mertens_nudge.rda", compress = "bzip2", version = 2)
