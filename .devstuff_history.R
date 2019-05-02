# the goal of this file is to have a trace of all devtools/usethis
# call you make for yout project

usethis::use_build_ignore(".devstuff_history.R")
usethis::use_build_ignore("inst/Covar_Selection/Covar_Selection.html")
usethis::use_build_ignore("inst/SDM_Selection/SDM_Selection.html")
usethis::use_git_ignore("inst/Covar_Selection/Covar_Selection.html")
usethis::use_git_ignore("inst/SDM_Selection/SDM_Selection.html")

# Depends: GlobalOptions, R (>= 3.3.0)

# Build vignettes
# 1/ Run vignettes that are in 'inst'
# 2/ Compress png images
system("mogrify -strip inst/SDM_Selection/*.png")
system("mogrify -strip inst/Covar_Selection/*.png")
# 3/ Build vignettes that are in vignettes (not run show images only for CRAN)
devtools::build_vignettes()

# Get dependencies
attachment::att_to_description(extra.suggests = c("pkgdown", "maptools"))
attachment::create_dependencies_file()

# Test site build
pkgdown::build_site()
