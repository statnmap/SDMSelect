#' Prepare dataset to be included in the model selection procedure
#'
#' This specifies data, covariates and allows resampling dataset in a regular grid for
#' autocorrelation purposes
#'
#' @param x data.frame or SpatialPointsDataFrame of observations with covariates
#' @param coords name or number of column with coordinates x,y.
#' Or a 2-column matrix of x-y coordinates. If missing, output object will be
#' a simple data.frame except if x is SpatialPointsDataFrame.
#' @param proj4string projection of the dataset. Default to x projection if exist or
#' '+proj=longlat +datum=WGS84'.
#' @param var column name or number of variable to be predicted
#' @param cov column name or numbers to be used as covariates
#' @param RefRaster raster as in library(raster). Raster used as smallest grid for data
#' gridding procedure. Resampling dataset into regular grid to decrease
#' spatial-autocorrelation. See details.
#' @param datatype string. Choose among "Cont" (continuous data, even if only positive),
#' "PA" (presence-absence data) or "Count" (count data). Required when RefRaster is
#' not NULL. dataY average for regular grid resampling are rounded for count models.
#' @param na.rm logical. If TRUE, removes complete dataset rows with NA values.
#' For further analysis, in particular with cross-validation procedure,
#' NA values in the dataset is a big problem. This will bias cross-validation
#' indices as they won't be calculated on the same amount of data
#' depending on covariates with NA in the model tested. Thus, although this is
#' drastic, rows of data with NA values are removed from the dataset with a warning.
#'
#' @details RefRaster recommendation : If there is spatial autocorrelation in the data,
#' use the higher resolution covariate raster as a reference or use \code{\link{spatialcor_dist}} to
#' determine smallest resolution to choose. If sampling plan is supposed not unbalanced
#' regarding covariates distribution, set RefRaster to NULL. Set to NULL if distribution
#' tested in the models will be KrigeGLM or KrigeGLM.dist (see \code{\link{AIC_indices}} or
#' \code{\link{crossV_indices}}).
#'
#'
#' @return return a dataset where variable of interest is called 'dataY'
#' (for compatibility with other functions. To be changed later to allow user defined name)
#' and where only covariates of interest are kept. 'factor_' is added to covariates
#' names that are factors (for compatibility with other functions).
#' @export

Prepare_dataset <- function(x, var = 1, cov = 2:ncol(x), coords,
                            proj4string = NULL, RefRaster = NULL,
                            datatype = "Cont", na.rm = TRUE
) {

  # Remove space in names
  names(x) <- make.names(names(x))
  if (length(grep("_P[[:digit:]]*$", names(x))) > 0) {
    names(x)[grep("_P[[:digit:]]*$", names(x))] <- paste0(names(x), "x")
    message(
      c("For compatibility reasons, names of covariates cannot finish by ",
        "'_P[:digit:]*'. Thus 'x' was added at the end of the problematic names"))
  }

  # Get information from a SpatialPointsDataFrame
  if (is(x)[1] == "SpatialPointsDataFrame") {
    proj4string <- sp::proj4string(x)
    coords <- x@coords #coordinates(x)
    test.na.x <- x@data
  } else if (grepl("MULTI", is(x)[1])) {
    stop("Spatial Object of library sf are not yet implemented,
         please provide a SpatialPointsDataFrame")
  } else {
    test.na.x <- x
    if (!missing(coords)) {
      if (is.vector(coords) & length(coords) == 2) {
        coords <- data.frame(x)[,coords]
      } else {
        coords <- data.frame(coords)
      }
    }
  }

  # Test for NA values in covariates and remove rows
  if (sum(is.na(test.na.x[,cov])) != 0) {
    w.na <- unique(row(test.na.x[,cov])[which(is.na(test.na.x[,cov]))])
    w.na <- w.na[order(w.na)]
    if (na.rm) {
      message(c("Your dataset contained covariates with NA values. ",
                length(w.na), " row(s) deleted for future cross-validation modelling"))
      if (!missing(coords)) {
        coords <- coords[-w.na,]
      }
      x <- x[-w.na,]
    } else {
      warning(c("Your dataset contains covariates with NA values : ", length(w.na), " row(s). ",
                "If you continue the analysis with the cross-validation procedure,",
                "you should set na.rm = TRUE in this function. See help."))
    }
  }


  #  Verify coordinates system when spatial data or coords supplied
  if (is.null(proj4string) & !missing(coords)) {
    warning(c("proj4string is not specified. Non projected system is supposed:",
              " +proj=longlat +datum=WGS84"))
    proj4string <- "+proj=longlat +datum=WGS84"
  }
  # If output is not spatial data
  if (missing(coords)) {
    message("coords empty, output dataset will not be a SpatialPointsDataFrame")
    Y_data_sp <- as.data.frame(cbind(x[,var], x[,cov]))
  } else {
    # If output is spatial data
    Y_data_sp <- sp::SpatialPointsDataFrame(coords = coords,
                                            data = as.data.frame(x)[,c(var, cov)],
                                            proj4string = sp::CRS(proj4string))
  }
  # First column is the y variable to be fitted
  names(Y_data_sp)[1] <- "dataY"

  if (grepl("PA", datatype) & any(!Y_data_sp$dataY %in% c(0,1))) {
    # Transform data as 0/1
    if (is.character(Y_data_sp$dataY) | is.factor(Y_data_sp$dataY)) {
      if (length(unique(Y_data_sp$dataY)) != 2) {
       stop("Only two levels are possible for presence-absence models")
      }
      levels <- levels(as.factor(Y_data_sp$dataY))
      message(levels[1], " was changed as 0, and ",
              levels[2], " was changed as 1")
      Y_data_sp$dataY <- as.numeric(as.factor(Y_data_sp$dataY)) - 1
    } else if (is.numeric(Y_data_sp$dataY)) {
    message("dataY was not 0/1, all numeric values above 0 were transformed as 1 for Presence-absence models")
    Y_data_sp$dataY <- 1 * (Y_data_sp$dataY != 0)
    } else {
      stop("Values in dataY were not numeric, character or factor with two levels.")
    }
  } else if (any(!is.numeric(Y_data_sp$dataY))) {
   stop("Values in dataY should be numeric for non PA models")
  }

  # List names of covariates
  namesCovariates <- names(Y_data_sp)[-which(names(Y_data_sp) == "dataY")]

  # All covariates that are factor are renamed 'factor' for further script
  # simplifications
  # IsData <- apply(t(namesCovariates), 2, function(x) is(Y_data_sp@data[, x])[1])
  if (missing(coords)) {
    IsNum <- unlist(lapply(Y_data_sp[,namesCovariates], is.numeric))
  } else {
    IsNum <- unlist(lapply(Y_data_sp@data[,namesCovariates], is.numeric))
  }
  if (length(which(!IsNum &
                   !grepl("factor_", namesCovariates))) > 0) {
    namesCovariates[which(!IsNum &
                            !grepl("factor_", namesCovariates))] <-
      paste0("factor_", namesCovariates[which(
        !IsNum & !grepl("factor_", namesCovariates))])
  }
  names(Y_data_sp)[-which(names(Y_data_sp) == "dataY")] <- namesCovariates

  # Resample the dataset using a raster grid
  if (is.null(RefRaster)) {
    Y_data_sample <- Y_data_sp
  } else {
    if (missing(coords)) {
      stop(c("Cannot resample the dataset into grid if coords is not specified",
             " or x is not SpatialPointsDataFrame"))
    }
    if (proj4string != sp::proj4string(RefRaster)) {
      warning(c("Dataset projection may be different than RefRaster one. ",
                "Resampling in the grid may not work"))
    }

    raster::values(RefRaster) <- 1:raster::ncell(RefRaster)
    # In which cells are samples
    Y_data_Grid <- raster::extract(RefRaster, Y_data_sp)

    # Mean position
    coordX_tmp <- tapply(sp::coordinates(Y_data_sp)[, 1], Y_data_Grid, mean)
    coordY_tmp <- tapply(sp::coordinates(Y_data_sp)[, 2], Y_data_Grid, mean)

    if (grep("PA", datatype)) {
      # Presence if at least one presence in cell
      dataY <- tapply(Y_data_sp$dataY, Y_data_Grid, max)
      modelselect_opt$datatype <- "PA"
    } else {
      # Mean of y variable in each cell
      dataY <- tapply(Y_data_sp$dataY, Y_data_Grid, mean)
      if (grep("Count", datatype)) {
        Y_data_sample$dataY <- round(Y_data_sample$dataY)
        message(c("Variable has been rounded to be used as count data ",
                  "with the regular grid resampling procedure"))
        if (!grepl("Count", modelselect_opt$datatype)) {
          warning(c("You should also choose one of the Count datatype for ",
                    "modelselect_opt$datatype"))
        }
      }
    }

    # Mean of covariates in each cell
    phys_tmp <- lapply(t(namesCovariates), function(j) {
      if (!grepl("factor_", j)) {
        phys_tmp_a <- tapply(
          as.numeric(as.character(Y_data_sp@data[,which(names(Y_data_sp) == j)])),
          Y_data_Grid, function(x) mean(x, na.rm = TRUE))
      } else {
        # Class covariates - The first of the most frequent class in the grid cell is
        # chosen
        phys_tmp_a <- tapply(
          as.character(Y_data_sp@data[,which(names(Y_data_sp) == j)]),
          Y_data_Grid,
          function(x) {
            if (length(x) != 0) {
              # factors are transformed as numeric factors for model selection
              # Most frequent factor
              res <- utils::tail(names(sort(table(x))), 1)
            } else {
              res <- NA
            }
            res
          }) %>% as.character() %>% as.factor()
      }
      phys_tmp_a
    })
    phys_tmp <- as.data.frame(phys_tmp)
    names(phys_tmp) <- namesCovariates

    for (i in grep("factor_", namesCovariates)) {
      phys_tmp[, i] <- as.factor(phys_tmp[, i])
    }

    # Save new dataset for the grid size chosen
    Y_data_sample <- sp::SpatialPointsDataFrame(
      coords = as.data.frame(cbind(coordX_tmp, coordY_tmp)),
      data = as.data.frame(cbind(dataY, phys_tmp)),
      proj4string = CRS(sp::proj4string(Y_data_sp)))
  }

  if (!missing(coords)) {
    if (!is.null(geoR::dup.coords(sp::coordinates(Y_data_sample)))) {
      # Message to me:
      # First thought that using as.data.frame on dataset led to rounded coordinates
      # This is false. Original data is kept safe.
      # However, running dup.coords on a dataset with format data.tbl does not work.
      # Thus, non need to round coordinates to number of digits = 15 (maximum number of digits
      # without computing problems : see print.default)
      warning(c("There are duplicated coordinates, ",
                "this will prevents the use of Kriging if one of your plans"))
    }
  }
  attr(Y_data_sample, "obs.col") <- 1
  return(Y_data_sample)
}


#' Create a RefRaster for regular grid dataset resampling
#'
#' @param x SpatialPointsDataFrame
#' @param res raster resolution targeted. Vector of one or two values for x and y resolutions.
#' @param degree Unit of the resolution. TRUE for degrees, FALSE for meters. Default to FALSE
#' Size one or two, if different resolution for x and y.
#'
#'
#' @importFrom raster projection xmin xmax ymin ymax
#'
#' @export

RefRasterize <- function(x, res, degree = FALSE) {

  # Two values resolution
  res2 <- (length(res) == 1) * c(res, res) + (length(res) != 1) * res

  if (grepl("Spatial", is(x)[1])) {
    ext <- raster::extent(x)
    crs <- raster::projection(x)
    # Calculate resolution according to projection
    if (is.na(raster::projection(x))) {
      warning("x has no projection, res is supposed to be in the same unit than x")
    } else {
      # If Same unit, no change res is res
      # Different unit
      if (!sp::is.projected(x) & !degree) {
        # Calculate average distance
        pt <- c(xmin(x), (ymin(x) + ymax(x)) / 2)
        pts <- cbind(xmax(x), (ymin(x) + ymax(x)) / 2)
        dist.max.x <- sp::spDistsN1(pts, pt, longlat = TRUE) * 1000
        if (res2[1] > dist.max.x) {
          stop(c("x resolution is bigger than extent in x-coordinates. ",
                 "Extent in x in meters is ", dist.max.x))
        }
        pt <- c((xmin(x) + xmax(x)) / 2, ymin(x))
        pts <- cbind((xmin(x) + xmax(x)) / 2, ymax(x))
        dist.max.y <- sp::spDistsN1(pts, pt, longlat = TRUE) * 1000
        if (res2[2] > dist.max.y) {
          stop(c("y resolution is bigger than extent in y-coordinates. ",
                 "Extent in y in meters is ", dist.max.y))
        }

        # nb of columns
        n.col <- floor(dist.max.x / res2[1])
        res2[1] <- (xmax(x) - xmin(x)) / n.col
        # nb of rows
        n.row <- floor(dist.max.y / res2[2])
        res2[2] <- (ymax(x) - ymin(x)) / n.row

        warning(c("Your dataset is in geographic coordinates and res is in meter, ",
                  "resolution has been recalculated in degrees. res = ",
                  paste(round(res2, digits = 7), collapse = ", ")))
      }
      if (sp::is.projected(x) & degree) {
        x.m <- data.frame(c(xmin(x),xmax(x)), c(ymin(x),ymax(x)))
        names(x.m) <- c("x", "y")
        sp::coordinates(x.m) <- ~x+y
        sp::proj4string(x.m) <- sp::proj4string(x)
        x.deg <- sp::spTransform(x.m, CRSobj = sp::CRS("+proj=longlat"))

        dist.max.x <- xmax(x.deg) - xmin(x.deg)
        dist.max.y <- ymax(x.deg) - ymin(x.deg)
        if (res2[1] > dist.max.x) {
          stop(c("x resolution is bigger than extent in x-coordinates. ",
                 "Extent in x in meters is ", dist.max.x))}
        if (res2[2] > dist.max.y) {
          stop(c("y resolution is bigger than extent in y-coordinates. ",
                 "Extent in y in meters is ", dist.max.y))}

        # nb of columns
        n.col <- round(dist.max.x / res2[1])
        res2[1] <- (xmax(x) - xmin(x)) / n.col
        # nb of rows
        n.row <- round(dist.max.y / res2[2])
        res2[2] <- (ymax(x) - ymin(x)) / n.row

        warning(c("Your dataset is projected with coordinates in meters,",
                  "resolution has been recalculated in meters. res = ",
                  paste(res2, collapse = ", ")))
      }
    }
  } else {
    names(x) <- c("x", "y")
    sp::coordinates(x) <- ~x+y
    crs <- NA
    warning(c("x is a data.frame, no projection will be assigned. ",
              "res is supposed to be in the same unit than x"))
  }

  ext <- raster::extent(x)
  nCol <- trunc((ext[2] - ext[1])/res2[1]) + 2
  nRow <- trunc((ext[4] - ext[3])/res2[2]) + 2
  r <- raster::raster(xmn = ext[1] - res2[1],
                      xmx = ext[1] - res2[1] + nCol * res2[1],
                      ymn = ext[3] - res2[2],
                      ymx = ext[3] - res2[2] + nRow * res2[2], ncol = nCol,
                      nrow = nRow, crs = sp::proj4string(x), vals = 1:(nRow * nCol))
  r
}


#' @title Figure to see spatial autocorrelation that may allow to define the grid
#' size for gridded procedure
#'
#' @description Figure to see spatial autocorrelation that may allow to define the grid
#' size for gridded procedure. Two figures are produced with different distance step.
#'
#' @param x A SpatialPointsDataFrame
#' @param y The column number on which to calculate correlation.
#' Otherwise a column should be names "dataY"
#' @param longlat logical. longlat data or not. Default to FALSE
#' @param max1 Numeric. maximum distance of panel 1 in meters.
#' Default to max distance divided by 3.
#' @param lag1 Numeric. step distance for calculating autocorrelation in meters.
#' Default to max1/100.
#' @param max2 Numeric. maximum distance of panel 1 in meters. Default to max1/10.
#' @param lag2 Numeric. step distance for calculating autocorrelation in meters.
#' Default to max2/100.
#' @param binomial Logical. Presence-absence data (TRUE) or continuous data (FALSE).
#' @param thd logical. If TRUE, the function suggest a distance threshold according to
#' breakpoint in correlation. This is only a suggestion. nls function is used.
#' @param plot Logical.
#' @param saveWD directory where to save the output figure. If null, figure appears on screen.
#' @param figname character. If saveWD is not empty, you can specify the name of
#' the output figure (without extension). Default to "Correlogram".
#' @param simplify.grid logical. If dataset is too big, it will be difficult to calculate
#' all distances between all points. A simplification is to divide the area into a grid
#' and calculate distances and correlation in each cell of the grid separately.
#' Results are then merged. Grid used depends on the largest of max1 and max2 values.
#' If less than 10 values, there are randomly merged with another cell
#'
#' @importFrom grDevices png jpeg
#'
#' @return A figure of correlation is showed or saved (in saveWD).
#' A threshold is suggested if thd=TRUE
#' @export

spatialcor_dist <- function(x, y, longlat = FALSE, max1, lag1, max2,
                            lag2, binomial = TRUE, thd = TRUE,
                            plot = TRUE, saveWD, figname,
                            simplify.grid = FALSE) {

  if (is(x)[1] != "SpatialPointsDataFrame") {
    stop("x must be a SpatialPointsDataFrame")
  }
  if (sum(grepl("dataY", names(x))) == 0 & missing(y)) {
    stop("One column should be named 'dataY' or y should be specified,\nplease use 'Prepare.dataset' function if needed")
  }
  if (!missing(y)) {names(x)[y] <- "dataY"}

  if (simplify.grid) {
    if (missing(max1) & missing(max2)) {
      stop("In the simplify.grid procedure, max1 or max2 must be specified !")
    } else {
      if (missing(max1)) {res <- max2}
      if (missing(max2)) {res <- max1}
      if (!missing(max1) & !missing(max2)) {res <- max(max1, max2)}
    }
    x.r <- RefRasterize(x = x, res = res)
    x.cells <- raster::extract(x.r, x)
    # Find cells with less than 10 values
    w.less <- which(x.cells %in% names(table(x.cells))[which(table(x.cells) <= 10)])
    x.cells[w.less] <- sample(names(table(x.cells))[which(table(x.cells) > 10)],
                              size = length(w.less), replace = TRUE)
  } else {
    x.cells <- 1
  }

  weight <- numeric(0)
  meanbyDist_all <- numeric(0)
  meanbyDistB_all <- numeric(0)

  for (grid in unique(x.cells)) { # grid <- unique(x.cells)[10]

    if (simplify.grid) {
      x.grid <- x[which(x.cells == grid),]
    } else {
      x.grid <- x
    }

      if (!sp::is.projected(x.grid)) {
        longlat <- TRUE
      }
      xyDist <- sp::spDists(x.grid, longlat = longlat)
      if (longlat) {xyDist <- xyDist * 1000} # in meters

    diag(xyDist) <- NA
    xyDist2 <- as.data.frame.table(xyDist)
    xyDist2[, 1] <- as.numeric(as.factor(xyDist2[, 1]))
    xyDist2[, 2] <- as.numeric(as.factor(xyDist2[, 2]))

    if (missing(max1)) {
      max1 <- max(xyDist, na.rm = TRUE)/3
    }
    if (missing(lag1)) {
      lag1 <- max1/100
    }
    if (missing(max2)) {
      max2 <- max(xyDist, na.rm = TRUE)/10
    }
    if (missing(lag2)) {
      lag2 <- max2/100
    }

    # Divide distances into classes
    seqdist <- seq(0, max1, lag1)
    cutxyDist <- cut(xyDist2[, 3], breaks = seqdist, labels = seqdist[-length(seqdist)])
    seqdistB <- seq(0, max2, lag2)
    cutxyDistB <- cut(xyDist2[, 3], breaks = seqdistB, labels = seqdistB[-length(seqdistB)])

    # test if pres is the same between stations %in% c('PA','PASeuil','PAGLM')
    if (binomial) {
      testequal <- (x.grid$dataY[xyDist2[, 1]] == x.grid$dataY[xyDist2[, 2]])
    } else {
      testequal <- ((x.grid$dataY[xyDist2[, 1]] - x.grid$dataY[xyDist2[, 2]]) ^ 2) /
        length(x.grid$dataY)
    }
    # Calculate the similarity into each class
    meanbyDist <- tapply(testequal, cutxyDist, function(x) mean(x, na.rm = TRUE))
    meanbyDistB <- tapply(testequal, cutxyDistB, function(x) mean(x, na.rm = TRUE))
    # sdbyDistB <- tapply(testequal, cutxyDistB, function(x) sd(x, na.rm = TRUE))

    if (simplify.grid) {
      weight <- c(weight, length(which(x.cells == grid)))

      meanbyDist_all <- cbind(meanbyDist_all, meanbyDist)
      meanbyDistB_all <- cbind(meanbyDistB_all, meanbyDistB)

    }
  }
  rm(xyDist, xyDist2); gc()
  # Regroup outputs of loop
  # Average according to number of points
  if (simplify.grid) {
    meanbyDist <- apply(meanbyDist_all, 1, function(x) {
      if (sum(is.na(x)) != length(x)) {
        sum(x * weight, na.rm = TRUE) / sum(weight[!is.na(x)])
      } else {
        NA
      }
    })
    meanbyDistB <- apply(meanbyDistB_all, 1, function(x) {
      if (sum(is.na(x)) != length(x)) {
        sum(x * weight, na.rm = TRUE) / sum(weight[!is.na(x)])
      } else {
        NA
      }
    })
  }

  if (binomial) {
    ymin <- -0.01
    ymax <- 1.05
    laby <- "Correlation"
  } else {
    ymin <- min(c(meanbyDist, meanbyDistB), na.rm = TRUE)
    ymax <- max(c(meanbyDist, meanbyDistB), na.rm = TRUE)
    laby <- "MSE"  # Mean squared error
  }

  # Find the recommended distance threshold according to break in correlation
  if (thd) {
    # thd 1
    dist <- as.numeric(as.character(names(meanbyDist))) + 0.5 * lag1
    lm.init <- lm(meanbyDist ~ dist)
    res <- apply(t(dist), 2, function(z) {
      aic <- NA
      try(aic <- AIC(nls(meanbyDist ~
                           (dist <= z) * (A * dist + B) +
                           (dist > z) * (D * (dist - z) + A * z + B),
                         start = list(A = coef(lm.init)[2],
                                      B = coef(lm.init)[1],
                                      D = 0))), silent = TRUE)
      aic
    })
    dist.min <- dist[which(res == min(res, na.rm = TRUE))]
    lm1 <- nls(meanbyDist ~ (dist <= dist.min) * (A * dist + B) + (dist > dist.min) *
                 (D * (dist - dist.min) + A * dist.min + B),
               start = list(A = coef(lm.init)[2],
                            B = coef(lm.init)[1], D = 0))
    # thd 2
    dist2 <- as.numeric(as.character(names(meanbyDistB))) + 0.5 * lag2
    lm.init <- lm(meanbyDistB ~ dist2)
    res <- apply(t(dist2), 2, function(z) {
      aic <- NA
      try(aic <- AIC(nls(meanbyDistB ~
                           (dist2 <= z) * (A * dist2 + B) +
                           (dist2 > z) * (D * (dist2 - z) + A * z + B),
                         start = list(A = coef(lm.init)[2],
                                      B = coef(lm.init)[1], D = 0))),
          silent = TRUE)
      aic
    })
    dist2.min <- dist2[which(res == min(res, na.rm = TRUE))]
    lm2 <- nls(meanbyDistB ~ (dist2 <= dist2.min) * (A * dist2 + B) +
                 (dist2 > dist2.min) * (D * (dist2 - dist2.min) +
                                          A * dist2.min + B),
               start = list(A = coef(lm.init)[2],
                            B = coef(lm.init)[1], D = 0))
  }


  # Plot the correlogram
  if (plot) {
    if (!missing(saveWD)) {
      if (!dir.exists(saveWD)) {dir.create(saveWD)}
      if (missing(figname)) {figname <- "Correlogram"}
      jpeg(filename = paste0(saveWD, "/", figname, ".jpg"), width = 16, height = 8, units = "cm",
           pointsize = 5, quality = 100, bg = "white", res = 300)
    } else {
      dev.new()
    }
    par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

    plot(as.numeric(as.character(names(meanbyDist))) + 0.5 * lag1, meanbyDist,
         pch = 20, cex = 1, ylim = c(ymin, ymax), xlim = c(0, max1), xaxs = "i",
         xlab = "distance (m)", ylab = laby,
         main = paste("Correlogram : lag=", round(lag1, digits = 2), "m", sep = ""))
    abline(v = c(lag1, 2 * lag1, 5 * lag1, 10 * lag1),
           col = "grey", lty = "dashed")
    axis(1, at = c(lag1, 2 * lag1, 5 * lag1, 10 * lag1),
         labels = c("lag1", "2x", "5x", "10x"), line = 1.5,
         col.axis = "grey40", col = "grey40")
    if (thd) {
      # thd1
      lines(dist, (dist <= dist.min) * (coef(lm1)[1] * dist + coef(lm1)[2]) +
              (dist > dist.min) * (coef(lm1)[3] * (dist - dist.min) + coef(lm1)[1] *
                                     dist.min + coef(lm1)[2]),
            col = "orange", lwd = 1.5)
      abline(v = dist.min)
      axis(1, at = dist.min, labels = paste("thd1 =", round(dist.min, digits = 1)),
           line = -0.8,
           col.axis = "grey40",
           col = "grey40")
    }

    plot(as.numeric(as.character(names(meanbyDistB))) + 0.5 * lag2, meanbyDistB,
         pch = 20, cex = 1, ylim = c(ymin, ymax), xlim = c(0, max2), xaxs = "i",
         xlab = "distance (m)", ylab = laby,
         main = paste("Correlogram : lag=", round(lag2, digits = 2), "m", sep = ""))
    abline(v = c(lag2, 2 * lag2, 5 * lag2, 10 * lag2), col = "grey", lty = "dashed")
    axis(1, at = c(lag2, 2 * lag2, 5 * lag2, 10 * lag2),
         labels = c("lag2", "2x", "5x", "10x"), line = 1.5,
         col.axis = "grey", col = "grey")
    if (thd) {
      # thd2
      lines(dist2, (dist2 <= dist2.min) *
              (coef(lm2)[1] * dist2 + coef(lm2)[2]) +
              (dist2 > dist2.min) *
              (coef(lm2)[3] * (dist2 - dist2.min) +
                 coef(lm2)[1] * dist2.min + coef(lm2)[2]),
            col = "orange", lwd = 1.5)
      abline(v = dist2.min)

      axis(1, at = dist2.min,
           labels = paste("thd2 =", round(dist2.min, digits = 2)),
           line = -0.8, col.axis = "grey40",
           col = "grey40")
    }
    if (!missing(saveWD)) {
      dev.off()
    }
  }
  if (thd) {
    thd.out <- c(dist.min, dist2.min)
    names(thd.out) <- paste0("thd",1:2)
    return(thd.out)
  }
}

#' Merge covariate rasters from different extents and projections into a multi-layer Raster
#'
#' @param cov.paths vector of paths of all raster files to be added in the multi-layer raster
#' @param r.ref Raster used as a reference for extent, resolution, crs.
#' Default is the first one of the list.
#' @param ext Optional object of class extent to crop r.ref before stack.
#' @param outWD Directory where to save new reprojected, cropped rasters.
#' Default to "tempdir()/extent".
#' @param names.cov Optional. vector of names of multi-layer raster. Same length as x.
#' @param is.factor position of covariates that are factors or categories.
#' Resampling (if required) is bilinear by default but "ngb" (nearest neighbour) for factors.
#' @param verbose level of messages shown
#' @param cl a cluster as made with \code{\link[parallel]{makeCluster}}.
#' If cl is empty, nbclust in \code{\link{modelselect_opt}} will be used
#' to create a cluster.
#'
#' @export
#'
#' @details Rasters are merged in a unique RasterStack. If the projection of
#' a raster is different than the reference (r.ref), it is reprojected. If the
#' extent is different, it is resampled. If the resolution is smaller than
#' 0.9*res, it is first aggregated, then resampled.

Prepare_covarStack <- function(cov.paths, r.ref, ext, outWD,
                                names.cov, is.factor, verbose = 1,
                                cl = NULL) {
  if (missing(outWD)) {
   outWD <- paste0(tempdir(), "/extent")
  }
  if (!dir.exists(outWD)) {
    dir.create(outWD)
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
  nbclust <- modelselect_opt$nbclust
  if (missing(r.ref)) {
    r.ref <- raster::raster(cov.paths[1])
  }
  if (!missing(ext)) {
    r.ref <- raster::crop(r.ref, ext)
  }
  method <- rep("bilinear", length(cov.paths))
  if (!missing(is.factor)) {
    method[is.factor] <- "ngb"
  }

  for (i in 1:length(cov.paths)) {
    r.x <- raster::raster(cov.paths[i])
    name.file <- strsplit(basename(cov.paths[i]),
                          raster::extension(cov.paths[i]), fixed = TRUE)[[1]][1]
    if (!missing(names.cov)) {
      name.x <- names.cov[i]
    } else {
      name.x <- name.file
    }
    file.ext <- paste0(outWD, "/", gsub("_extent", "", name.file), "_extent.tif")
    # Re-write projection if same proj written differently
    if (sp::proj4string(r.x) != sp::proj4string(r.ref)) {
      #' Method depends if projection is really different,
      #' otherwise this is only a change in name of proj
      if (is_ProjEqual(x = sp::proj4string(r.x), y = sp::proj4string(r.ref))) {
        sp::proj4string(r.x) <- sp::proj4string(r.ref)
        if (verbose >= 1) {
          message(paste("Equal proj modified", i, ":", name.x))
        }
      } else {
        # Use resample as shown below
        # compareRaster will return FALSE
        if (verbose >= 1) {
          message(paste("projectRaster", i, "=", name.x))
        }
        if (cl_inside) {
          raster::beginCluster(nbclust)
        }
        r.ext <- raster::projectRaster(
          r.x,
          crs = sp::proj4string(r.ref),
          method = method[i]
        )
        r.x <- r.ext
        rm(r.ext); gc()
        if (cl_inside) {
          raster::endCluster()
        }
      }
    }

    # Resample rasters if needed
    if (!raster::compareRaster(r.x, r.ref, stopiffalse = FALSE)) {
      if (file.exists(file.ext)) {
        print(paste("load resampled", i, "=", name.x))
        r.ext <- raster::raster(file.ext)
      } else {
        if (cl_inside) {
          raster::beginCluster(nbclust)
        }
        if (raster::xres(r.x) < 0.9 * raster::xres(r.ref) |
            raster::yres(r.x) < 0.9 * raster::yres(r.ref)) {
          # Aggregate before resampling if resolution is far from the target
          if (verbose >= 1) {
            message(paste("aggregate", i, "=", name.x,
                          "with", method[i]))
          }
          r.x <- raster::aggregate(
            r.x,
            fact = trunc(c(raster::xres(r.ref)/raster::xres(r.x),
                           raster::yres(r.ref)/raster::yres(r.x))),
            fun = function(x, na.rm = na.rm) {
              if (method[i] == "bilinear") {
                res <- mean(x, na.rm = na.rm)
              } else {
                # Most frequent factor
                res <- as.numeric(as.character(
                  utils::tail(names(sort(table(x))), 1)))
                if (length(res) == 0) {
                  res <- NA
                }
              }
              res
            },
            na.rm = TRUE)
        }
        if (verbose >= 1) {
          message(paste("resample", i, "=", name.x, "with", method[i]))
        }
        r.ext <- raster::resample(r.x, r.ref, method = method[i],
                                  filename = file.ext,
                                  overwrite = TRUE)
        if (cl_inside) {
          raster::endCluster()
        }
      }
    } else {
      r.ext <- r.x
    }
    if (i == 1) {
      r.multi <- r.ext
    } else {
      r.multi <- raster::stack(r.multi, r.ext)
    }
  }
  if (!missing(names.cov)) {
    names(r.multi) <- names.cov
  }
  return(r.multi)
}


#' @title Function to extract covariates from different rasters
#'
#' @description Original covariates raster may have different projections.
#' CovarExtract allows to extract in these different projections.
#'
#' @param x SpatialPointsDataFrame
#' @inheritParams Prepare_covarStack
#'
#' @export

CovarExtract <- function(x, cov.paths, names.cov, is.factor) {
  if (is(x)[1] != "SpatialPointsDataFrame") {
    stop("x should be a SpatialPointsDataFrame")
  }
  res.list <- lapply(1:length(cov.paths), function(i) {
    r.tmp <- raster::raster(cov.paths[i])
    if (!is_ProjEqual(sp::proj4string(r.tmp), sp::proj4string(x))) {
      x.proj <- sp::spTransform(x, raster::crs(r.tmp))
    } else {x.proj <- x}
    res <- data.frame(raster::extract(r.tmp, x.proj))
    names(res) <- names(r.tmp)
    res
  })
  res.tab <- do.call("cbind", res.list)
  if (!missing(names.cov)) {
    names(res.tab) <- names.cov
  }
  if (!missing(is.factor)) {
    for (i in is.factor) {
      res.tab[,i] <- as.factor(as.character(res.tab[,i]))
    }
  }
  sp::cbind.Spatial(x, res.tab)

}


#' @title Function to identify covariates that are correlated in terms of Spearman's rank
#'
#' @description Correlation engender problems of identifiability.
#' Correlated parameters in the dataset will be separated in tested models.
#' Test for the Spearman factor for non linear correlation between covariates
#' of all the stations.
#' Complete the test with a visual test if needed.
#'
#'
#' @param x dataframe of covariates only, on which to test correlation of covariates.
#' This can also be dataset as issued from \code{\link{Prepare_dataset}}.
#' @param rm vector of column numbers to be removed from the analysis. Default to NULL.
#' If you specify your complete dataset as \code{x}, you may define \code{rm=1} to
#' remove your observations from the columns.
#' @param visual logical. Whether to define visually if data are considered
#' correlated or not. See details. Be careful, all previous figures are
#' closed with graphics.off() before running visual analysis.
#' @param thd numeric. Correlation (absolute) value above which to consider that covariates
#' are correlated and should not remain in the same model. This value is necessary
#' when visual = FALSE. See details.
#' @param plot logical. Whether to plot the figure of correlation values between
#' covariates. Set to FALSE if using a MPI cluster.
#' @param saveWD path to directory where to save a jpg figure file. If NULL and
#' plot = TRUE, then figure is shown on screen.
#' @param figname character. The name (w/o extension) of the figure to be saved in saveWD.
#' @param img.size size in cm of the output (square) image. May increase labels size as well.
#'
#' @importFrom magrittr %>%
#' @importFrom methods is
#' @importFrom grDevices graphics.off dev.new dev.set dev.prev dev.off png
#' @import graphics
#'
#' @details
#' \itemize{
#' \item Correlation of covariates is to be tested in the dataset itself. It is not
#' important if covariates are correlated in real life. What affect model fits
#' is the correlation within the dataset. Therefore is this function...
#' \item Spearman's rank correlation coefficient has been chosen as it allows to
#' test for correlation between categorical and continuous variables. This may
#' look as a non-sense to test for correlation with categorical covariates.
#' In some cases, the order of classes inside a category may have a sense, e.g.
#' continuous variable that has been categorised for any reason. In this function
#' classes are temporary turned into numbers to calculate correlation.
#' This implies that the alphanumeric order of classes has a sense.
#' The Spearman test is also less sensitive to non-gaussian distributions of data
#' \item In some cases, correlation value is high because of one extreme rank value.
#' Visual verification allows to define for each couple of covariates if the user
#' may consider correlation or not. Figures are grouped by correlation values
#' which accelerates the visual verification.
#' \item It is suggested to do a visual verification the first time data are
#' processed or to do a complete real exploration of your dataset before
#' running any model. This exploration may suggest to remove or combine
#' covariates before starting the modelling procedure. After that, you may be
#' able to define a thd value for covariate correlation limit.
#' \item Defining a high thd value may let some correlated covariates to appear in
#' the same model in the following of the procedure. Nevertheless, because the
#' procedure of this package is based on cross-validation, if two variables are
#' correlated they may likely have a low score as they will be as efficient as
#' covariate alone. Thus, if you hesitate between two thd values, choose the
#' higher one, as it will allow more models to pass through the selection procedure.
#' }
#'
#' @export

Param_corr <- function(x, rm = NULL, visual = FALSE, thd = 0.7, plot = TRUE, saveWD = NULL,
                       figname = "Covariate_correlation", img.size = 12)
{

  if (!is.null(attr(x, "obs.col"))) {
    rm <- unique(c(rm, attr(x, "obs.col")))
  }

  # if (visual)  {require(tcltk)}

  if (is(x)[1] == "SpatialPointsDataFrame") {
    if (!is.null(rm)) {
      physParam_NewSp <- x@data[,-rm]
    } else {
      physParam_NewSp <- x@data
    }
  } else {
    if (!is.null(rm)) {
      physParam_NewSp <- x[,-rm]
    } else {
      physParam_NewSp <- x
    }
  }

  # All combinations of 2 covariates
  # Change names so that names order reflect column order
  names.save <- names(physParam_NewSp)
  names(physParam_NewSp) <-
    paste0("p", formatC(1:ncol(physParam_NewSp),
                        width = nchar(ncol(physParam_NewSp)),
                        flag = "0"), names(physParam_NewSp))
  # Change factor as numeric for correlation only
  if (any(sapply(physParam_NewSp, function(x) !is.numeric(x)))) {
    message("Variables are transformed as numeric for correlation tests only. ",
            "This may not be optimal, but this will not alter the rest of the analysis")
  }

  fig.tmp <- physParam_NewSp %>%
    dplyr::mutate_at(.vars = grep("^factor_", names(physParam_NewSp)),
                     .funs = function(x) as.numeric(as.factor(x))) %>%
    dplyr::mutate_if(sapply(., function(x) !is.numeric(x)), as.numeric)# %>%
    # dplyr::mutate_if(sapply(., is.factor), as.numeric)

  # Calculate correlation with corrr
  corSpearman.tmp <- fig.tmp %>%
    corrr::correlate(method = "spearman") %>%
    dplyr::mutate_all(function(x) ifelse(is.na(x), 1, x))
  class(corSpearman.tmp) <- c("cor_df", class(corSpearman.tmp))

  # dplyr 0.8.0
  corSpearman <- corSpearman.tmp %>%
    corrr::shave() %>%
    tidyr::gather(key = "Var2", value = "Corr", -term) %>%
    dplyr::rename(Var1 = term) %>%
    #fashion(leading_zeros = TRUE) %>%
    # as.data.frame.table() %>%
    # filter(Freq != "") %>%
    dplyr::mutate_at(
      .vars = 1:2,
      # .funs = dplyr::funs("num" = as.numeric(as.factor(.)))) %>%
      .funs = list("num" = ~as.numeric(as.factor(.x)))) %>%
    dplyr::filter(!is.na(Corr)) %>%
    as.data.frame()

  # Test visually each couple of covariates using tcltk
  if (visual)
  {
    warning("Visual analysis uses library(tcltk) which may not work
            properly using Rstudio.")
    graphics.off()
    cutSpearman <- cut(abs(corSpearman$Corr), breaks = seq(0,1,0.1))
    keep <- logical(0)
    placeX <- logical(0)

    # Test correlations in decreasing correlation
    for (val in length(levels(cutSpearman)):4) {
      # 4 is for no test under correlation < 0.3 which is too low for any corr
      val1 <- levels(cutSpearman)[val]
      X <- which(cutSpearman == val1)
      placeX <- c(placeX, X)
      if (length(X) != 0) {
        seq.fig <- c(seq(1, length(X), 25), length(X))
        for (n.fig in 1:(length(seq.fig) - 1)) {
          dev.new()
          # Draw bowplot of couples of covariates
          if (length(X) <= 6) {par(mfrow = c(2,3), ask = TRUE)}
          if (length(X) > 6) {par(mfrow = c(3,4), ask = TRUE)}
          if (length(X) > 12) {par(mfrow = c(5,5), ask = TRUE)}
          for (ii in seq.fig[n.fig]:seq.fig[n.fig + 1]) {
            if (length(X) > 25 & ii == seq.fig[n.fig] + 2) {
              mtext(text = paste("Figure", n.fig), side = 3,
                    line = 1, outer = TRUE)
            }
            boxplot(rank(physParam_NewSp[,unlist(corSpearman[ii,4])],
                         na.last = "keep") ~
                      round(rank(physParam_NewSp[,unlist(corSpearman[ii,5])],
                                 na.last = "keep")/20),
                    main = val1,
                    xlab = names(physParam_NewSp)[unlist(corSpearman[ii,4])],
                    ylab = names(physParam_NewSp)[unlist(corSpearman[ii,5])])
          }
          if (length(X) > 12) {
            dev.set(dev.prev())
          }
        } # n.fig
        # Fonction to ask if corr is retained or not
        tt <- tcltk::tktoplevel()
        tl <- tcltk::tklistbox(tt, height = 4, selectmode = "single",
                               background = "white")
        tcltk::tkgrid(tcltk::tklabel(tt, text = "Is it overall visually correlated ?"))
        tcltk::tkgrid(tl)
        fruits <- c("yes_all", "no_none", "yes_somes", "")
        for (i in (1:4))
        {
          tcltk::tkinsert(tl, "end", fruits[i])
        }
        tcltk::tkselection.set(tl,3)

        OnOK <- function()
        {
          rbVal <- fruits[as.numeric(tcltk::tkcurselection(tl)) + 1]
          tcltk::tkdestroy(tt)
          if (rbVal == "yes_all") {tcltk::tkmessageBox(
            message = "All is considered correlated")
            assign("truc", 1)}
          if (rbVal == "no_none") {tcltk::tkmessageBox(
            message = "Nothing is considered correlated")
            assign("truc", 0)}
          if (rbVal == "yes_somes") {tcltk::tkmessageBox(
            message = "You need to choose which are correlated,
            put 1 in the table")
            assign("truc", -1)}
        }
        truc <- NA
        OK.but <- tcltk::tkbutton(tt,text = "   OK   ", command = OnOK)
        tcltk::tkgrid(OK.but)
        tcltk::tkfocus(tt)

        repeat {if (is.na(truc) == FALSE) {break}}

        # Assign data chosen to the list of data
        if (truc == -1) {
          ind <- cbind(names(physParam_NewSp)[unlist(corSpearman[X,4])],
                       names(physParam_NewSp)[unlist(corSpearman[X,5])], 1)
          utils::fix(ind)
          keep <- c(keep, ind[,3])
        } else {
          keep <- c(keep, rep(truc, length(X)))
        }
        } # if X != 0
      }

    par(ask = FALSE)
  } # end of visual == TRUE


  # Store whether the parameters are correlated or not
  corSpearman$valid <- 0
  if (visual) {
    corSpearman$valid[placeX] <- keep
  } else {
    # Parameters with spearman rank correlation higher than 0.8 wont be in the same model
    message(c("thd value is set to "), thd)
    corSpearman$valid[which(abs(corSpearman$Corr) >= thd)] <- 1
  }

  if (plot) {
    # Figure of correlations coefficients
    fig.cor <- fig.tmp %>% cor(method = "spearman")
    fig.cor[which(is.na(fig.cor))] <- 1

    colnames(fig.cor) <- rownames(fig.cor) <- names.save

    # Test with ggplot
    # library(ggplot2)
    # ggplot(data = corSpearman, aes(Var2, Var1, fill = Corr)) +
    #   geom_tile(color = "white") +
    #   scale_fill_gradient2(low = "blue", high = "red", mid = "white",
    #                        midpoint = 0, limit = c(-1,1), space = "Lab",
    #                        name = "Pearson\nCorrelation") +
    #   theme_minimal()+
    #   theme(axis.text.x = element_text(angle = 45, vjust = 1,
    #                                    size = 12, hjust = 1)) +
    #   coord_fixed()

    # Add shade lines for couples removed from models
    fig.cor.0 <- fig.cor * 0
    for (x in 1:nrow(corSpearman)) {
      fig.cor.0[corSpearman[x, "Var1_num"], corSpearman[x,"Var2_num"]] <-
        fig.cor.0[corSpearman[x,"Var2_num"], corSpearman[x,"Var1_num"]] <-
        corSpearman$valid[x]
    }
    # Do a loop to use suggested dev size of corrplot (less white margins)
    ratio <- 1
    # for (loop.fig in 1:2) {
    # To reset par()
    par(mar = c(0,0,0,0)); plot.new(); dev.off()
    if (!is.null(saveWD)) {
      png(filename = paste0(saveWD, "/", figname, ".png"),
          width = img.size, height = img.size * ratio, units = "cm", pointsize = 5,
          #quality = 100,
          res = 500
      )
    } else {
      dev.new()
    }
    # Correlation plot
    corrplot::corrplot(
      fig.cor, #mar = c(0, 0.5, 0.5, 0.5),
      tl.cex = 2.5 * min(img.size, img.size * ratio) / ncol(physParam_NewSp), method = "square",
      diag = FALSE,
      tl.srt = 45, # cex.main = 2.5,
      title = "",
      type = "lower",
      # tl.pos = "n",
      tl.col = "black",
      cl.cex = 1.5 * min(img.size, img.size * ratio) / ncol(physParam_NewSp) + ncol(physParam_NewSp) / 150,
      cl.ratio = 0.1,
      # h.ratio = TRUE,
      win.asp = ratio)
    # Keep or not keep couples
    corrplot::corrplot(
      fig.cor.0,
      tl.cex = img.size / ncol(physParam_NewSp),
      method = "shade",
      col = c("transparent"), addshade = "positive",
      shade.col = "grey",
      shade.lwd = 2 * min(img.size, img.size * ratio) / ncol(physParam_NewSp),
      diag = FALSE, add = TRUE, bg = "transparent",
      tl.pos = "n", cl.pos = "n",
      type = "lower", tl.col = "black")
    if (!is.null(saveWD)) {
      dev.off()

      # Remove white margins
      # if (loop.fig == 2) {
      Rm_WhiteMargins(paste0(saveWD, "/", figname, ".png"), pixel = 5)
      # }
    }
    #} # end of loop.fig
  } # end of plot
  corSpearman <- corSpearman %>%
    dplyr::mutate(Var1 = names.save[Var1_num]) %>%
    dplyr::mutate(Var2 = names.save[Var2_num])

  return(corSpearman)
} # End of ParamCorr function


