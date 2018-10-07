# SDMSelect 0.1.4.9000

# SDMSelect 0.1.4

* Added pkgdown page on github
* Allow for factors in presence-absence models (transformed as 0/1)
* Correct options with version 0.1.0 of {GlobalOptions}
* fixed `NaN` could be returned when comparing distribution with itself in `best_distri`.
* Reduce number of outputs saved at each iteration to only the best ones. With a lot of covariates tested, saving all model tested was too huge. Too keep more models, `modelselect_opt$lim_pvalue` can be lowered.
* Save outputs of model selection as `rds`, not `RData`.
* Added `gplot_data`, a modified version of `rasterVis::gplot` to retrieve data of raster for visualisation instead of directly plotting it in ggplot.