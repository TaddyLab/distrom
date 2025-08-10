

# Current version: 0.1.0

# Update release process:
# Alter code, documentation, testthat, vignette, and NEWS appropriately.
# Pass expectation tests, local R-CMD-Check, and remote R-CMD-Checks.
# Update cran-comments.md, release to CRAN with devtools::release()
# After code has been accepted to CRAN, push new code to GitHub, update GitHub release page.

library(devtools)
library(available)
library(rhub)
library(testthat)

setwd("/Users/nelsonrayl/Desktop/init/lmForc/lmForc")

#===============================================================================
# Compile Updated Package
#===============================================================================

# Load package.
devtools::load_all("/Users/nelsonrayl/Desktop/init/lmForc/lmForc")

# Compile package documentation.
devtools::document()

# Build vignette.
# Note* to edit the vignette, edit the file /vignettes/lmForc.Rmd
devtools::build_vignettes()
vignette("lmForc", package = "lmForc")

# Re-build package.
devtools::build("/Users/nelsonrayl/Desktop/init/lmForc/lmForc")

#===============================================================================
# Test Package
#===============================================================================

# Run expectation tests.
devtools::test()

# Run R-CMD-Check locally.
devtools::check()

# Run R-CMD-Check on external platforms.
platforms <- c(
  "fedora-clang-devel", 
  "ubuntu-gcc-release", 
  "windows-x86_64-devel",
  "macos-m1-bigsur-release"
)

rhub::check(
  path = "/Users/nelsonrayl/Desktop/init/lmForc/lmForc",
  platform = platforms,
  email = "nelsonrayl14@gmail.com"
)

# Check R-CMD-Check on win_devel.
devtools::check_win_devel()

# Test reverse dependency for downstream packages.
#?devtools::revdep

#===============================================================================
# Release Package
#===============================================================================

# Final spell check.
devtools::spell_check()

# Release package to CRAN.
devtools::release()

# Re-submit package to CRAN without going through all release() questions.
#devtools::submit_cran()

#===============================================================================
# GitHub
#===============================================================================

# Add R-CMD-Check badge to GitHub.
# If the R-CMD-Check process updates, run this command again to update the
# file /.github/workflows/R-CMD-check.yaml
usethis::use_github_action("check-standard")

#===============================================================================
# Monitor CRAN Downloads
#===============================================================================

library(cranlogs)

downloads <- cran_downloads("lmForc", from = "2021-10-11", to = "2022-12-20")
sum(downloads$count)






install.packages("badger")
library(badger)

badge_cran_download("badger", "grand-total", "blue")


