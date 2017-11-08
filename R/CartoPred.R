#' @title Predictions in response and link scales
#' @description Internal function to be coupled with \code{\link[raster]{predict}}
#' and \code{\link[raster]{clusterR}}
#'
#' @param object fitted model as output of \code{\link{glm}} or \code{\link{gam}}
#' @param newdata a data frame in which to look for variables with which to predict

predR <- function(object, newdata)
{
  v <- v2 <- list()
  v$fit <- v$se.fit <- v2$fit <- v2$se.fit <- NA
  if (grepl("glm.fit", object$method)) {
    try(v <- raster::predict(object, newdata, se.fit=TRUE, type="response"))
    try(v2 <- raster::predict(object, newdata, se.fit=TRUE, type="link"))
  } else {
    try(v <- mgcv::predict.gam(object, newdata, se.fit=TRUE, type="response"))
    try(v2 <- mgcv::predict.gam(object, newdata, se.fit=TRUE, type="link"))
  }
  cbind(p = as.vector(v$fit), se = as.vector(v$se.fit),
        p_l = as.vector(v2$fit), se_l = as.vector(v2$se.fit))
}


#' Create Multilayer raster for predictions and uncertainty analysis
#'
#' @param object Multi-layer Raster with at least covariates of the model selected
#' @param mask Optional Raster object to be used as mask for output raster: NA values
#' of the mask raster induce NA values in the predictions.
#' @inheritParams ModelResults
#'
#' @export

Map_predict <- function(object, saveWD, mask = NULL, Num = NULL,
                        model = NULL, modeltype = NULL, powerXI = NULL,
                        zip.file = TRUE, cl = NULL)
{

  if (utils::file_test("-f", saveWD) & file.exists(saveWD) & grepl("zip", saveWD)) {
    if (zip.file == TRUE) {zip.file <- saveWD}
    utils::unzip(saveWD, exdir = gsub(".zip", "", saveWD))
    saveWD <- gsub(".zip", "", saveWD)
  } else if (utils::file_test("-d", saveWD)) {
    if (zip.file == TRUE) {zip.file <- paste0(saveWD, ".zip")}
  } else {
    stop("saveWD is neither an existing directory nor a zip file")
  }

  if (is.null(cl)) {
    cl_inside <- TRUE
  } else {
    cl_inside <- FALSE
    options(rasterClusterObject = cl)
    options(rasterClusterCores = length(cl))
    options(rasterCluster = TRUE)
    options(rasterClusterExclude = NULL)
  }

  if (!is.null(Num) & !is.numeric(Num)) {
    Num <- as.numeric(as.character(Num))
  }

  saveWD <- normalizePath(saveWD)

  # Load information
  # load(paste0(saveWD, "/Allinfo_all.RData"))
  Allinfo_all <- readr::read_rds(paste0(saveWD, "/Allinfo_all.rds"))
  datatype <- Allinfo_all$datatype

  # Load variables of model configuration
  modelselect_opt(readr::read_rds(paste0(saveWD, "/modelselect_opt.save.rds")))

  # Options ----
  Y.max <- modelselect_opt$Y.max
  nbclust <- modelselect_opt$nbclust

  if (grepl("KrigeGLM", datatype)) {
    Y_data_sample_lcc <- spTransform(Allinfo_all$data, CRS(lcc_proj))
  } else {
    Y_data_sample_lcc <- NA
  }
  if (is(Allinfo_all$data)[1] == "SpatialPointsDataFrame") {
    Y_data_sample <- dplyr::as.tbl(Allinfo_all$data@data)
  } else {
    Y_data_sample <- dplyr::as.tbl(Allinfo_all$data)
  }

  # Define the model on which to produce outputs
  model_selected <- model_select(saveWD = saveWD,
                                 new.data = Y_data_sample,
                                 Num = Num, model = model,
                                 modeltype = modeltype, powerXI = powerXI)

  Num <- model_selected$Num
  modeltype <- model_selected$modeltype
  modelX <- model_selected$modelX
  Param <- model_selected$Param
  if (grepl("PA|Tweed", model_selected$modeltype)) {
    Seuil <- model_selected$Seuil
    # SeuilMinUn <- model_selected$SeuilMinUn
    # SeuilMaxZero <- model_selected$SeuilMaxZero
  }

  maxY <- max(Y_data_sample$dataY) * Y.max

  # Search for factors in modelX to modify names of object
  w.not <- which(!all.vars(modelX$formula) %in% names(object) & all.vars(modelX$formula) != "dataY")
  if (length(w.not) != 0) {
    for (w in w.not) {
      # Test for Name included in like factor_*
      w.r <- which(unlist(lapply(names(object), function(x) grepl(x, all.vars(modelX$formula)[w]))))
      if (length(w.r) == 1) {
        names(object)[w.r] <- all.vars(modelX$formula)[w]
      }
      if (length(w.r) > 1) {
        stop("Multiple similar names in data do not allow for object rename")
      }
      if (length(w.r) == 0) {
        stop(paste(all.vars(modelX$formula)[w], "or similar do not exist in object"))
      }
    }
  }

  # Remove unnecessary layers in object
  object <- raster::dropLayer(object,
                                   i = which(!names(object) %in%
                                               all.vars(modelX$formula)[-1]))
  ## Create Raster of masks for prediction from Y_data_sample


  ## Select thd with validation dataset

  # Predictions in response and link scales ----
  if (cl_inside) {
    raster::beginCluster(nbclust)
  }
  pred_map_x1 <- raster::clusterR(object, raster::predict,
                          args = list(model = modelX,
                                      fun = predR, index = 1:4), #ModelSelect:::predR
                          m = 3)
  if (cl_inside) {
    raster::endCluster()
  }
  gc()
  names(pred_map_x1) <- c("resp.fit", "resp.se.fit", "link.fit", "link.se.fit")

  # Back transform log data by "exp" with Laurent transformation ------------
  # if(Family %in% c("LogCont","LogContGLM")){
  if (grepl("Log", model_selected$modeltype)) {
    raster::values(pred_map_x1)[,1] <- exp(
      raster::values(pred_map_x1)[,1] + 0.5 * var(residuals(modelX)))
  }

  # Probability of presence as estimated by tweed model ---------------------
  if (grepl("Tweed", model_selected$modeltype)) {
    Tweed.proba.of.pres <- raster::raster(pred_map_x1, 1)
    names(Tweed.proba.of.pres) <- "Tweed.proba.of.pres"
    if (cl_inside) {
      cl <- parallel::makePSOCKcluster(nbclust)
    }
    w.noNA <- which(!is.na(raster::values(pred_map_x1)[,1]))
    predNoNA <- t(raster::values(pred_map_x1)[w.noNA, 1])
    # maxY <- max(Y_data_sample$dataY) * 2

    # Calculate probability of presence (y=0 in Tweedie is probability of absence)
    proba <- 1-(snow::parApply(cl, predNoNA, 2, function(x, Param, modelX, maxY) {
      if (x > maxY) {
        # res <- 1 # Too long calculation for too high values
        res <- tweedie::dtweedie(y = 0, xi = Param, mu = maxY, phi = summary(modelX)$dispersion)
      } else {
        #library(statmod)
        #library(tweedie)
        # x <- v$fit
        #if(is.na(x) == FALSE){
        ## v$fit <- object$family$linkinv(v2$fit)
        # res <- 1-ptweedie(q=0, xi=Param, mu=x, phi=summary(modelX)$dispersion, power=NULL)
        res <- tweedie::dtweedie(y = 0, xi = Param, mu = x, phi = summary(modelX)$dispersion)
        # tweedie::dtweedie(y = 0, xi = Param, mu = 0.0000001, phi = summary(modelX)$dispersion)
      }
      res
    }, Param = Param, modelX = modelX, maxY = maxY))
    if (cl_inside) {
      parallel::stopCluster(cl)
    }
    gc()
    raster::values(Tweed.proba.of.pres)[w.noNA] <- proba
    pred_map_x1 <- raster::addLayer(pred_map_x1, Tweed.proba.of.pres)
    rm(Tweed.proba.of.pres); gc()
  }

  # CV calculation -----------------------------------------------------------
  if (grepl("PA", model_selected$modeltype)) {
    pred_map_CV <- abs(raster::raster(pred_map_x1, layer = 4) /
                         raster::raster(pred_map_x1, layer = 3))
    names(pred_map_CV) <- "CV.link.scale"
  } else {
    pred_map_CV <- abs(raster::raster(pred_map_x1, layer = 2) /
                         raster::raster(pred_map_x1, layer = 1))
    names(pred_map_CV) <- "CV"
  }
  pred_map_x1 <- raster::addLayer(pred_map_x1, pred_map_CV)
 rm(pred_map_CV); gc()

  # Quantiles thus uncertainty of predictions for each raster cell ----
  # To be seen later...

  # title Calculate quantiles of output predictions to evaluate uncertainty
  # description Predict function for a glm or gam outputs mean and sd. To get an
  # idea of the uncertainty of the assessment, this function calculates the
  # min (5\%) and max (95\%) quantiles of the estimation of a prediction from
  # estimates in the scale of the link function.
  #
  # param x dataframe of two columns of mean and sd estimates in the scale of
  # the link function
  # param args list of additional arguments like the model used, threshold for
  # presence-absence data. This parameter is supplied as a list to make the function
  # compatible with \code{\link[raster]{calc}} and \code{\link[raster]{clusterR}}
  # export

  # predinc <- function(x, args) {
  #   if (!missing(args)) {
  #     try(Seuil <- args$Seuil)
  #     try(modelX <- args$modelX)
  #     try(Family <- args$Family)
  #
  #   }
  #   if (missing(Family)) {
  #     if (!missing(modelX)) {
  #       Family <- "gaussian"
  #       if (sum(grepl("binomial", modelX$family)) != 0) {
  #         Family <- "binomial"
  #       }
  #       if (sum(grepl("Tweed", modelX$family)) != 0) {
  #         Family <- "binomial"
  #       }
  #     }
  #   }
  #   b <- rep(NA, 3)
  #   if (Family == "binomial") {
  #     try(b <- 1/(1 + exp(-qnorm(c(0.05, 0.5, 0.95), mean = x[1], sd = x[2]))))
  #   }
  #
  # #    if (exists(Seuil) ){# & Family %in% c("binomial", "Tweed")) {
  # #   if (Family %in% c("binomial", "Tweed")) {
  #     if (!missing(Family)) {
  #       ProbaSup <- NA
  #       try(ProbaSup <- 100*(1-pnorm(log(Seuil/(1-Seuil)), mean = x[1], sd = x[2])))
  #     }
  # #   }
  #
  #   # pnorm(log(Seuil/(1-Seuil)), mean = -1.4, sd = 4)
  # #   res <- cbind(b[1], b[2], b[3])
  # #   if (exists(ProbaSup)) {
  #     res <- cbind(b[1], b[2], b[3], ProbaSup)
  # #   }
  #   return(res)
  # }
  #
  #       link.rast.tmp <- stack(raster(pred_map_x1, 3), raster(pred_map_x1, 4))
  #
  #       if (cl_inside) {beginCluster(nbclust)}
  #         z1 <- clusterR(link.rast.tmp, calc, args=list(fun=predinc), export=c('Seuil', 'modelX'))
  #       if (cl_inside) {endCluster()}
  #

  if (grepl("PA", model_selected$modeltype)) {

    predincPA <- function(y, S) raster::calc(y, function(x) {
      b <- rep(NA, 5)
      ProbaSup <- NA
      try(b <- 1/(1 + exp(-qnorm(c(0.05, 0.25, 0.5, 0.75, 0.95), mean = x[3], sd = x[4]))))
      try(ProbaSup <- 100*(1-pnorm(log(S/(1-S)), mean = x[3], sd = x[4])))
      cbind(b[1], b[2], b[3], b[4], b[5], ProbaSup)
    })

    if (cl_inside) {raster::beginCluster(nbclust)}
    pred_map_inc <- raster::clusterR(pred_map_x1, predincPA, m = 3,
                             args = list(S = Seuil))
    if (cl_inside) {raster::endCluster()}

    # Probleme of sd calculation appearing with GAM
    if (!grepl("GLM", model_selected$modeltype)) {
      for (i in 1:raster::nlayers(pred_map_inc)) {
        raster::values(pred_map_inc[[i]])[
          which(raster::values(pred_map_x1[[1]]) < 1E-10 |
                  raster::values(pred_map_x1[[4]]) >= 5E4)] <- 0
      }
    }
  }

  if (!grepl("PA|Tweed", model_selected$modeltype)) {
    # predictions in the scale of the link function as uncertainty is gaussian
    predinc <- function(y) raster::calc(y, function(x) {
      b <- rep(NA,5)
      try(b <- qnorm(c(0.05, 0.25, 0.5, 0.75, 0.95), mean = x[3], sd = x[4]))
      cbind(b[1], b[2], b[3], b[4], b[5])
    })

    if (cl_inside) {raster::beginCluster(nbclust)}
    system.time(pred_map_inc <- raster::clusterR(pred_map_x1, predinc, m = 3))
    if (cl_inside) {raster::endCluster()}

    # Back transform by "exp"
    if (grepl("Log", model_selected$modeltype)) {
      # With Laurent Correction as data are log(dataY)
      for (i in 1:raster::nlayers(pred_map_inc)) {
        raster::values(pred_map_inc)[,i] <- exp(
          raster::values(pred_map_inc)[,i] + 0.5*var(residuals(modelX)))
      }
    }
    if (grepl("Gamma", model_selected$modeltype)) {
      # values(pred_map) <- (exp(values(pred_map))-1)
      # No Need of Laurent Correction as data are not log(dataY)
      # but only link function is log and values are quantiles
      for (i in 1:raster::nlayers(pred_map_inc)) {
        raster::values(pred_map_inc)[,i] <- exp(raster::values(pred_map_inc)[,i])
      }
    }
  }

  if (grepl("Tweed", model_selected$modeltype)) {

    predincTweed <- function(x, Param, modelX, maxY) {
      if(x > maxY){
        x <- maxY # Too long calculation for too high values
      }
      # library(statmod);library(tweedie)
      tweedie::qtweedie(p = c(0.05, 0.25, 0.5, 0.75, 0.95), xi = Param, mu = x,
                        phi = summary(modelX)$dispersion)
      # tweedie::qtweedie(p=c(0.05,0.5,0.95),xi=Param,mu=1,phi=summary(modelX)$dispersion)
    }

    pred_map_inc <- raster::raster(pred_map_x1,1)
    for (i in 2:5) {
      pred_map_inc <- raster::stack(pred_map_inc, raster::raster(pred_map_x1,1))
    }
    if (cl_inside) {
      cl <- parallel::makePSOCKcluster(nbclust)
    }
    predNoNA <- t(raster::values(pred_map_x1)[which(
      !is.na(raster::values(pred_map_x1)[,1])),1])
    pred_map_inc_tmp <- snow::parApply(
      cl, predNoNA, 2,
      function(x, Param, modelX, maxY) {
        if (x > maxY) {
          x <- maxY
          # Too long calculation for too high values
        }
        res <- tweedie::qtweedie(p = c(0.05, 0.25, 0.5, 0.75, 0.95), xi = Param, mu = x,
                                 phi = summary(modelX)$dispersion)
        # tweedie::qtweedie(p=c(0.05,0.5,0.95),xi=Param,mu=0.057,phi=summary(modelX)$dispersion)
        return(res)
      }, Param = Param, modelX = modelX, maxY = maxY)
    if (cl_inside) {
      parallel::stopCluster(cl)
    }
    raster::values(pred_map_inc)[
      which(!is.na(raster::values(pred_map_x1)[,1])),1:5] <-
      t(pred_map_inc_tmp)[,1:5]

    # Not ready
    # 	      predinc <- function(y,G,P){
    # 		library(statmod)
    # 		library(tweedie)
    # 		raster::calc(y, function(x){
    # 		  b <- rep(NA,3)
    # 		  # ProbaSup <- NA
    # 		  # try(b <- 1/(1 + exp(-qnorm(c(0.05,0.5,0.95),mean=x[3], sd=x[4]))))
    # 		  try(b <- qtweedie(p=c(0.05,0.5,0.95),xi=P,mu=x[1],phi=G))
    # 		  # try(ProbaSup <- 100*(1-pnorm(log(S/(1-S)),mean=x[3], sd=x[4])))
    # 		  cbind(b[1],b[2],b[3])#,ProbaSup)
    # 		})}
    #     if (cl_inside) {
    # 	    raster::beginCluster(nbclust)
    #     }
    # 		system.time(
    # 		  pred_map_inc <- clusterR(pred_map,predinc,m=3,args=list(G=summary(modelX)$dispersion,P=Param)) # S=Seuil,
    # 		)
    #      if (cl_inside) {
    # 	      raster::endCluster()
    #       }
  }
  # IQR and IQR/median

  pred_map_inc <- raster::addLayer(pred_map_inc,
                                   raster::raster(pred_map_inc, 4) -
                                     raster::raster(pred_map_inc, 2),
                           (raster::raster(pred_map_inc, 4) -
                              raster::raster(pred_map_inc, 2)) /
                             raster::raster(pred_map_inc, 3))
  if (!grepl("PA", model_selected$modeltype)) {
    names(pred_map_inc) <- c("Q5", "Q25", "Q50", "Q75", "Q95", "IQR", "IQR.M")
  } else {
    names(pred_map_inc) <- c("Q5", "Q25", "Q50", "Q75", "Q95", "ProbaSup", "IQR", "IQR.M")
  }
  pred_map_x1 <- raster::addLayer(pred_map_x1, pred_map_inc)
  rm(pred_map_inc); gc()


  # Calculate sure Presence and absence according to SeuilMinUn and SeuilMaxZero
  # if (sum(grepl("binomial", modelX$family)) != 0 |
  # sum(grepl("Tweed", modelX$family)) != 0) {
  #       raster::values(pred_map5a)[which(values(pred_map_inc)[,3] <= SeuilMinUn)] <- -1
  #       raster::values(pred_map5a)[which(values(pred_map_inc)[,1] >= SeuilMaxZero)] <- 1
  # }

  # Calculate mask for parameters in the scope of the data ---------------------
  # as predictions should not be calculated out of the scope of the data
  for (mask.i in all.vars(as.formula(model_selected$formulaX))[-1]) {

    mask.tmp <- raster::raster(pred_map_x1, 1) * NA
    #     mask.tmp <- raster(object, 1) * NA

    if (is.numeric(data.frame(model_selected$new.data)[,mask.i])) {
      raster::values(mask.tmp)[which(
        raster::values(object)[,which(names(object) == mask.i)] >=
          min(data.frame(model_selected$new.data)[,mask.i]) &
          raster::values(object)[,which(names(object) == mask.i)] <=
          max(data.frame(model_selected$new.data)[,mask.i]))] <- 1
    } else {
      raster::values(mask.tmp)[which(
        raster::values(object)[,
                            which(names(object) == mask.i)] %in%
          unique(data.frame(model_selected$new.data)[,mask.i]))] <- 1
      #           &
      #           !is.na(values(object)[,
      #             which(paste0("factor_", names(object)) == mask.i)]))] <- 1
    }
    names(mask.tmp) <- paste0("mask_", mask.i)
    pred_map_x1 <- raster::addLayer(pred_map_x1, mask.tmp)
    # pred_map_x1 <- dropLayer(pred_map_x1, 14:27)
    rm(mask.tmp)
    gc()
  }

  # Apply mask on prediction maps
  # --------------------------
  if (!is.null(mask)) {
    pred_map <- raster::mask(pred_map_x1, mask)
  } else {
    pred_map <- pred_map_x1
  }

  raster::writeRaster(pred_map, filename = paste0(saveWD, "/prediction_raster.grd"),
              overwrite = TRUE)
  pred_map <- raster::stack(paste0(saveWD, "/prediction_raster.grd"))

  # Zip all outputs in a unique file.
  if (zip.file != FALSE) {
    # Options j allows to save directly files without path
    utils::zip(zip.file, files = saveWD, flags = "-urj9X", extras = "",
        zip = Sys.getenv("R_ZIPCMD", "zip"))
    # Return the file path of the zip file
    message(paste("all results saved in", zip.file))
  }
  return(pred_map)
} # end of function

#' Transform raster as data.frame to be later used with ggplot
#' Modified from rasterVis::gplot
#'
#' @param x A Raster* object
#' @param maxpixels Maximum number of pixels to use
#'
#' @details rasterVis::gplot is nice to plot a raster in a ggplot but
#' if you want to plot different rasters on the same plot, you are stuck.
#' If you want to add other information or transform your raster as a
#' category raster, you can not do it. With `SDMSelect::gplot_data`, you retrieve your
#' raster as a data.frame that can be modified as wanted using `dplyr` and
#' then plot in `ggplot` using `geom_tile`.
#'
#' @export

gplot_data <- function(x, maxpixels = 50000)  {
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  coords <- raster::xyFromCell(x, seq_len(raster::ncell(x)))
  ## Extract values
  dat <- raster::stack(as.data.frame(raster::getValues(x)))
  names(dat) <- c('value', 'variable')

  dat <- dplyr::as.tbl(cbind(coords, dat))
  dat
}