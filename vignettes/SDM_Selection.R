## ------------------------------------------------------------------------
library(SDMSelect)
library(ggplot2)
require(dismo)

set.seed(50)
# Temp directory for saving all outputs
tmpdir <- paste0(tempdir(), "/out_SDM")
dir.create(tmpdir)

## ------------------------------------------------------------------------
# World map for plots
data(wrld_simpl, package = "maptools")

# Vectors of paths to rasters of covariates
covariates.paths <- 
  list.files(file.path(system.file(package = "dismo"), 'ex'),
             pattern = 'grd$',
             full.names = TRUE)

# Here covariates <- stack(covariates.paths) 
# would have been enough for this specific case.
is.factor <- grep("biome.grd", covariates.paths)

covariates <- Prepare_covarStack(cov.paths = covariates.paths, 
                              is.factor = is.factor)


## ------------------------------------------------------------------------
print(covariates)

## ------------------------------------------------------------------------
# Presence data
file <- file.path(system.file(package = "dismo"), "ex/bradypus.csv")
bradypus <- read.table(file,  header = TRUE,  sep = ',')[,-1]
# Random-absence data
randAbs <- dismo::randomPoints(covariates, 400)
colnames(randAbs) <- c("lon", "lat")
# Combine presence and absence
data.sp <- rbind(data.frame(randAbs, Obs = 0),
                 data.frame(bradypus, Obs = 1)
              )
# Transform data as SpatialPointsDataFrame
sp::coordinates(data.sp) <- ~lon+lat
sp::proj4string(data.sp) <- "+proj=longlat +datum=WGS84"

# Extract covariates, combine with dataset and set factor covariate
data <- CovarExtract(x = data.sp, cov.paths = covariates.paths,
                     is.factor = is.factor)
  

## ---- out.width='50%'----------------------------------------------------
# Show observations positions
par(mar = c(1,1,1,1))
raster::plot(covariates, "bio1")
sp::plot(wrld_simpl, add = TRUE)
sp::plot(data, 
     col = c("red", "blue")[data@data$Obs + 1],
     pch = 20, cex = c(0.5, 1)[data@data$Obs + 1],
     add = TRUE)

## ------------------------------------------------------------------------
data.prepared <- Prepare_dataset(
  x = data, var = 1, cov = 2:ncol(data),
  datatype = "PA", na.rm = TRUE
)

## ------------------------------------------------------------------------
thd <- spatialcor_dist(
  x = data.prepared, longlat = !is.projected(data),
  binomial = TRUE, saveWD = tmpdir,
  plot = TRUE,
  figname = "Correlogram"
)

## ---- out.width='50%'----------------------------------------------------
# Create a regular grid from dataset
RefRaster <- RefRasterize(x = data, res = round(thd[2]))
# Use rectangular grid to resample dataset
data.new <- Prepare_dataset(
  x = data.prepared, var = 1, cov = 2:ncol(data),
  RefRaster = RefRaster, datatype = "PA", na.rm = TRUE
)

# Plot data.new
par(mar = c(1,1,1,1))
raster::plot(covariates, "bio1")
sp::plot(wrld_simpl, add = TRUE)
sp::plot(data.new, 
     col = c("red", "blue")[data.new@data$dataY + 1],
     pch = 20, cex = c(0.5, 1)[data@data$Obs + 1],
     add = TRUE)

## ------------------------------------------------------------------------
corSpearman <- Param_corr(
  x = data.new, rm = 1, thd = 0.7, visual = FALSE,
  plot = TRUE, img.size = 8, saveWD = tmpdir)


## ------------------------------------------------------------------------
modelselect_opt(RESET = TRUE)
modelselect_opt$Max_nb_Var <- 2
modelselect_opt$datatype <- "PA"

## ------------------------------------------------------------------------
res.file <- findBestModel(x = data.new, datatype = "PA", 
                          corSpearman = corSpearman, 
                          saveWD = tmpdir, 
                          verbose = 1)

## ------------------------------------------------------------------------
# Less discriminant test to select the best models
modelselect_opt$lim_pvalue_final <- 1e-3

# Order models and find the bests
BestModels <- ModelOrder(saveWD = tmpdir, plot = TRUE)

## ---- message=FALSE, results='hide'--------------------------------------
Num.Best <- BestModels$VeryBestModels_crossV$Num[1]
res.file <- ModelResults(saveWD = tmpdir, plot = TRUE, 
                 Num = Num.Best)

## ------------------------------------------------------------------------
pred.r <- Map_predict(object = covariates, saveWD = tmpdir, Num = Num.Best)

## ---- out.width='60%', fig.height=6--------------------------------------
rasterVis::gplot(raster::raster(pred.r, "resp.fit")) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient("Probability", low = 'yellow', high = 'blue') +
  coord_equal()

## ---- out.width='90%', fig.width=9, fig.height=4.5-----------------------
rasterVis::gplot(raster::dropLayer(pred.r, which(!names(pred.r) %in% c("Q5", "Q95")))) +
  geom_tile(aes(fill = value)) +
  facet_grid(~variable) +
  scale_fill_gradient("Probability", low = 'yellow', high = 'blue') +
  coord_equal()

## ---- out.width='60%', fig.height=6--------------------------------------
rasterVis::gplot(raster::dropLayer(pred.r, which(!names(pred.r) %in% c("IQR")))) +
  geom_tile(aes(fill = value)) +
  facet_grid(~variable) +
  scale_fill_gradient("Absolute\nDispersion", low = 'white', high = 'red') +
  coord_equal()

## ---- out.width='60%', fig.height=6--------------------------------------
rasterVis::gplot(raster::dropLayer(pred.r, which(!names(pred.r) %in% c("IQR.M")))) +
  geom_tile(aes(fill = value)) +
  facet_grid(~variable) +
  scale_fill_gradient("Relative\nDispersion\nto median", low = 'white', high = 'red') +
  coord_equal()

## ------------------------------------------------------------------------
model_selected <- model_select(
  saveWD = tmpdir,
  new.data = data.new,
  Num = Num.Best)
BestThd <- model_selected$Seuil

## ---- out.width='60%', fig.height=6--------------------------------------
rasterVis::gplot(raster::raster(pred.r, "resp.fit")) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2("Probability\nof\nPresence",
                       low = 'white', mid = 'yellow',  high = 'blue',
                       midpoint = BestThd) +
  coord_equal()

## ---- out.width='60%', fig.height=6--------------------------------------
rasterVis::gplot(raster::raster(pred.r, "ProbaSup")) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient2("Probability\nto be over\nThreshold", 
                      low = 'white', mid = 'yellow', high = 'forestgreen',
                      midpoint = 50) +
  coord_equal()

## ---- out.width='90%', fig.width=9, fig.height=4.5-----------------------
rasterVis::gplot(raster::dropLayer(pred.r, which(!grepl("mask", names(pred.r))))) +
  geom_tile(aes(fill = factor(value))) +
  facet_grid(~variable) +
  scale_fill_manual("mask", values = c("0" = "red", "1" = "forestgreen")) +
  coord_equal()

