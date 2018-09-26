## ------------------------------------------------------------------------
library(SDMSelect)
set.seed(50)

# Temp directory for saving all outputs
tmpdir <- paste0(tempdir(), "/out_CovarSelection")
dir.create(tmpdir)

## ------------------------------------------------------------------------
data <- dplyr::mutate_at(mtcars, 8:11, as.character)

## ------------------------------------------------------------------------
data.new <- Prepare_dataset(
  x = data, var = 1, cov = 2:ncol(data),
  datatype = "Cont", na.rm = TRUE
)

## ------------------------------------------------------------------------
corSpearman <- Param_corr(
  x = data.new, rm = 1, thd = 0.7, visual = FALSE,
  plot = TRUE, saveWD = tmpdir, img.size = 5)


## ------------------------------------------------------------------------
modelselect_opt(RESET = TRUE)
modelselect_opt$Max_nb_Var <- 5
modelselect_opt$datatype <- "ContPosNull"
modelselect_opt$modeltypes <- modelselect_opt$modeltypes[c(1, 5, 11)]

## ---- eval=FALSE---------------------------------------------------------
#  res.file <- findBestModel(
#    x = data.new, datatype = "ContPosNull",
#    corSpearman = corSpearman, saveWD = tmpdir,
#    verbose = 1)

## ---- eval=FALSE---------------------------------------------------------
#  BestModels <- ModelOrder(saveWD = res.file, plot = TRUE)

## ---- echo=FALSE---------------------------------------------------------
BestModelsFile <- system.file("Covar_Selection", "BestModels.rds", package = "SDMSelect")
BestModels <- readr::read_rds(BestModelsFile)

## ---- message=FALSE, results='hide', eval=FALSE--------------------------
#  Num.Best <- BestModels$VeryBestModels_crossV$Num[1]
#  res.file <- ModelResults(saveWD = tmpdir, plot = TRUE,
#                   Num = Num.Best, Marginals = TRUE)

