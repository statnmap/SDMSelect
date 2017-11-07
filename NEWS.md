# SDMSelect 0.1.3.9000

- fixed `NaN` could be returned when comparing distribution with itself in `best_distri`.
- Reduce number of outputs saved at each iteration to only the best ones. With a lot of covariates tested, saving all model tested was too huge. Too keep more models, `modelselect_opt$lim_pvalue` can be lowered.
- Save outputs of model selection as `rds`, not `RData`.