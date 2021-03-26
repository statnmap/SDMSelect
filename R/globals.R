# To get rid of "no visible binding for global variable" when CHECK
# ModelResults
utils::globalVariables(c(".", "pred", "na.row", "pred.absmean", "cov.val",
                         "fit", "se.fit", "dataY"))
# Param_corr
utils::globalVariables(c("term", "Corr", "Var1_num", "Var2_num"))
# Fig_split
utils::globalVariables(c("ntot", "diff.x", "diff.sum"))

# All variables that are called using load.
# Will be removed when all variables are called from rds.
utils::globalVariables(c("Allinfo_all", "lcc_proj", "saveAlea", "resAIC_save",
                         "resParam_save", "test_models", "Output_models",
                         "BestTHD_crossV", "DiffSelSpeTHD_crossV", "Model",
                         "saveAlea.nb", "ToutAUC_crossVMean", "fix.lambda",
                         "BestTHD_crossV_unlist",
                         "MinUn_crossV", "MaxZero_crossV"))