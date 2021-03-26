#' Function to fit a model among all possible
#'
#' @param x data.frame or SpatialPointsDataFrame of observations with covariates
#' @param x_lcc Projected version of x for kriging
#' @param model model formula written as character
#' @param modeltype sub-model type
#' @param MaxDist Maximum distance for variogram ("Krige*" modeltype only)
#' @param Phi Range for Phi fitting ("Krige*" modeltype only)
#' @param fixXI Power of the Tweedie model (\code{\link[statmod]{tweedie}})
#' @param Model Model type for variogram ("Krige*" modeltype only)
#' @param fix.lambda logical, indicating whether the Box-Cox transformation parameter
#' lambda should be regarded as fixed (fix.lambda = TRUE) or should be be estimated
#' (fix.lambda = FALSE). Defaults to TRUE.
#' @param lambda value of the Box-Cox transformation parameter lambda.
#' Regarded as a fixed value if fix.lambda = TRUE otherwise as the initial value
#' for the minimisation algorithm. Defaults to 1. Two particular cases are
#' lambda = 1 indicating no transformation and lambda = 0 indicating log-transformation.
#'
#' @importFrom splines ns
#'
#' @export

fit_model <- function(x, model, fixXI,
                  modeltype,
                  x_lcc,
                  MaxDist, Phi, Model,
                  fix.lambda, lambda) {
  modelX <- NA
  res1 <- NA
  res2 <- NA
  if (!grepl("GLM", modeltype)) {
    if (grepl("PA", modeltype)) {
      try(modelX <- mgcv::gam(as.formula(model), method = "REML", data = x,
                      family = "binomial", control = list(maxit = 200)), silent = TRUE)
    } else if (grepl("Cont", modeltype)) {
      try(modelX <- mgcv::gam(as.formula(model), method = "REML", data = x,
                      control = list(maxit = 500)), silent = TRUE)
    } else if (grepl("Gamma", modeltype)) {
      try(modelX <- mgcv::gam(as.formula(model), method = "REML", data = x,
                      family = Gamma(link = "log"),
                      control = list(maxit = 500)), silent = TRUE)
    } else if (grepl("Count", modeltype)) {
      # try(modelX <- gam(as.formula(model), method = "REML", data = x,
      #                 family = "poisson",
      #                 control = list(maxit = 500)),silent = TRUE)
      # try(modelX <- gam(as.formula(model),method='REML',
      # data=Y_data_adjust,family='poisson',control=list(maxit=500)),silent=TRUE)
      try(modelX <- mgcv::gam(as.formula(model), method = "REML", data = x,
                      family = "poisson"), silent = TRUE)
      Maxit <- 100
      Epsilon <- 1e-07
      cnt <- 0
      if (is.na(modelX) == TRUE) {
        modelX$converged <- FALSE
      }
      while (modelX$converged == FALSE & cnt <= 5) {
        cnt <- cnt + 1
        # print(cnt) print(Maxit) print(Epsilon)
        Maxit <- Maxit + 500
        Epsilon <- Epsilon * 5
        try(modelX <- mgcv::gam(as.formula(model), method = "REML", data = x,
                        family = "poisson", control = list(epsilon = Epsilon, maxit = Maxit)),
            silent = TRUE)
        # print(modelX$converged)
      }
      if (cnt > 5) {
        modelX <- NA
      }
    }
  } else {
    if (grepl("PA", modeltype)) {
      try(modelX <- glm(as.formula(model), data = x, family = "binomial"),
          silent = TRUE)
    } else if (grepl("Cont", modeltype)) {
      try(modelX <- glm(as.formula(model), data = x, family = "gaussian"),
          silent = TRUE)
    } else if (grepl("Gamma", modeltype)) {
      try(modelX <- glm(as.formula(model), data = x, family = Gamma(link = "log")),
          silent = TRUE)
    } else if (grepl("Count", modeltype)) {
      # try(modelX <- glm(as.formula(model), data = x, family = "poisson"),
      #     silent = TRUE)
      try(modelX <- glm(as.formula(model), data = x, family = "poisson"),
          silent = TRUE)
      Maxit <- 25
      Epsilon <- 1e-08
      cnt <- 0
      if (is.na(modelX) == TRUE) {
        modelX$converged <- FALSE
      }
      while (modelX$converged == FALSE & cnt <= 5) {
        cnt <- cnt + 1
        # print(cnt) print(Maxit) print(Epsilon)
        Maxit <- Maxit + 500
        Epsilon <- Epsilon * 5
        try(modelX <- glm(as.formula(model), data = x, family = "poisson",
                        control = list(epsilon = Epsilon, maxit = Maxit)), silent = TRUE)
        # print(modelX$converged)
      }
      if (cnt > 5) {
        modelX <- NA
      }
    } else if (grepl("TweedGLM", modeltype)) {
      if (fixXI == 0) {
        # Test if convergence is possible on simple case before tweedie.profile
        a <- try(modelX <- glm(as.formula(model), data = x,
                             family = statmod::tweedie(var.power = 1.5, link.power = 0),
                             control = list(maxit = 100)),
                 silent = TRUE)
        if (is.na(modelX)[1]) {
          res1 <- NA
        } else {
          if (!modelX$converged) {
            res1 <- NaN
            res2 <- fixXI
          } else {
            try(modelX1 <- tweedie::tweedie.profile(as.formula(model), data = x,
                                         xi.vec = seq(1.1, 1.9, 0.1), do.plot = FALSE, fit.glm = FALSE,
                                         do.ci = FALSE, method = "series"), silent = TRUE)
            a <- try(modelX <- glm(
              as.formula(model), data = x,
              family = statmod::tweedie(var.power = modelX1$xi.max, link.power = 0),
              control = list(maxit = 100)),
              silent = TRUE)
            try(res1 <- tweedie::AICtweedie(modelX), silent = TRUE)
            try(res2 <- modelX1$xi.max[1])
          }
        }
      } else {
        # a <- try(modelX <- glm(as.formula(model), data = x,
        #                      family = tweedie(var.power = fixXI, link.power = 0)), silent = TRUE)
        options(warn = -1)
        try(modelX <- glm(as.formula(model), data = x,
                        family = statmod::tweedie(var.power = fixXI, link.power = 0)),
            silent = TRUE)
        Maxit <- 25
        Epsilon <- 1e-08
        cnt <- 0
        if (is.na(modelX)) {
          modelX$converged <- FALSE
        }
        while (!modelX$converged & cnt <= 5) {
          cnt <- cnt + 1
          Maxit <- Maxit + 500
          Epsilon <- Epsilon * 5
          try(modelX <- glm(as.formula(model), data = x,
                          family = statmod::tweedie(var.power = fixXI, link.power = 0),
                          control = list(epsilon = Epsilon, maxit = Maxit)),
              silent = TRUE)
        }
        if (cnt > 5) {
          modelX <- NA
        }
        options(warn = 0)

        if (is.na(modelX)[1]) {
          res1 <- NA
        } else {
          if (!modelX$converged) {
            res1 <- NA
          } else {
            try(res1 <- tweedie::AICtweedie(modelX), silent = TRUE)
          }
        }
        res2 <- fixXI
      }

    } else if (grepl("KrigeGLM", modeltype) &
               !grepl("KrigeGLM.dist", modeltype)) {
      datageo <- geoR::as.geodata(
        x_lcc,
        data.col = which(names(x_lcc) == "dataY"),
        covar.col = c(1:(ncol(x_lcc)))[-which(names(x_lcc) == "dataY")]
      )
      # Modification to make it work with poly()
      datageo$covariate <- tibble::as_tibble(datageo$covariate)
      try(data.v <- geoR::variog(datageo, max.dist = MaxDist, trend = as.formula(model),
                           breaks = seq(0, MaxDist, length = 20), lambda = lambda, messages = FALSE),
          silent = TRUE)
      try(data.mlIni <- geoR::variofit(data.v, cov.model = Model, limits = geoR::pars.limits(phi = Phi),
                                 messages = FALSE), silent = TRUE)
      try(Ini <- summary(data.mlIni)$estimated.pars[-1], silent = TRUE)
      try(Nug <- data.mlIni$nugget, silent = TRUE)
      try(modelX <- geoR::likfit(datageo, trend = as.formula(model), cov.model = Model,
                         ini = Ini, fix.nugget = FALSE, nugget = Nug, limits = geoR::pars.limits(phi = Phi),
                         lik.met = "REML", lambda = lambda, fix.lambda = fix.lambda, messages = FALSE),
          silent = TRUE)
    }
  } # end of GLM
  list(modelX = modelX, res1 = res1, res2 = res2)
}

#' Function that returns the AIC of a list of models
#' In the case of a tweedie model (TweedGLM), it also returns the XI value to be used in cross-validation if not fixed previously
#'
#' @param x the model number to be fitted. From 1 to length(Models_tmp_nb).
#' @param Y_data_sample data.frame or SpatialPointsDataFrame of observations with covariates
#' @param Models_tmp_nb matrix with column of model formulas as character
#' @param Y_data_sample_lcc dataset to be fit on, with projected CRS
#' ("Krige*" modeltype only)
#' @inheritParams fit_model
#'
#' @importFrom stats AIC
#'
#' @export

AIC_indices <- function(x, Y_data_sample, Models_tmp_nb, modeltype, # libfile = .libPaths(),
                        fixXI, Y_data_sample_lcc = NA, MaxDist = NA, Phi = NA,
                        Model, fix.lambda, lambda) {
  options(error = expression(NULL))

  # For box-cox or Log transformation, data should be strictly positive Add +1 to
  # data if needed
  Sup <- 0
  if ((!is.na(lambda) & !identical(lambda, 1)) |
      grepl("Log|Gamma", modeltype)) {
    if (min(Y_data_sample$dataY) < 0) {
      stop("For Box-Cox or Log transformation or Gamma model,\n data should be strictly positive")
    }
    if (min(Y_data_sample$dataY) == 0) {
      Sup <- 1
      warning("For box-cox or Log transformation or Gamma model, +1 was added to data")
    }
  }
  Y_data_sample$dataY <- Y_data_sample$dataY + Sup
  if (!is.na(Y_data_sample_lcc)) {
    Y_data_sample_lcc$dataY <- Y_data_sample_lcc$dataY + Sup
  }

  # data.frame has been modified so that poly is not anymore usable ------------
  Y_data_sample <- tibble::as_tibble(as.data.frame(Y_data_sample))

  formula_tested <- as.character(Models_tmp_nb[x, 1])
  res <- NA
  res1 <- NA
  res2 <- NA

  # Fit model ----
  model.out <- fit_model(
    x = Y_data_sample,
    model = formula_tested, fixXI = fixXI,
    modeltype = modeltype,
    x_lcc = Y_data_sample_lcc,
    MaxDist = MaxDist, Phi = Phi, Model = Model,
    fix.lambda = fix.lambda, lambda = lambda)

  modelX <- model.out$modelX
  res1 <- model.out$res1
  res2 <- model.out$res2

  try(res <- AIC(modelX), silent = TRUE)
  if (grepl("Log", modeltype)) {
    try(res <- AIC(modelX) + 2 * sum(log(Y_data_sample$dataY)), silent = TRUE)
  }
  if (!is.na(lambda) & !is.na(modelX[1])) {
    if (fix.lambda) {
      try(res1 <- AIC(modelX) - 2 * sum(log((Y_data_sample$dataY)^(lambda - 1))),
          silent = TRUE)
      res2 <- lambda
    } else {
      try(res1 <- AIC(modelX) - 2 * sum(log((Y_data_sample$dataY)^(modelX$lambda -
                                                                   1))), silent = TRUE)
      try(res2 <- modelX$lambda)
    }
  }

  options(error = NULL)
  if (grepl("TweedGLM", modeltype) | !is.na(lambda)) {
    res <- c(res1, res2)
  }
  res
}


#' Function that return indices of quality of fit for a list of models with a
#' sub-sample of data used as validation sample in a cross-validation procedure
#'
#' @param MC index of the cross-validation subset to achieve
#' @param formulas vector of model formulas to be tested
#' @param saveAlea list of indices of data observations used for validation
#' @param Y_data_sample dataset on which to run cross-validation
#' @param resParam_save Vector of length of formulas with special parameter for
#' Tweedie or Krige* models as calculated in \code{\link{AIC_indices}}
#' @param Y_data_sample_lcc Projected dataset in meters for Krige* models
#' @param seqthd Sequence of thresholds tested to cut between 0 and 1 for PA data.
#' @inheritParams fit_model
#'
#'
#' @return
#' if: modeltype in 'PA','PAGLM','PASeuil','TweedGLM'
#'  then: resAIC,resUBRE,resDev,ResDev_crossV,minUn,maxZero, MeanTHD,DiffSelSpe,ROC_crossV,MSE_crossV,Logl)
#' if: in Cont, Density, KrigeGLM
#'  then: resAIC,resUBRE,resDev,ResDev_crossV,MeanPercentError,MSE_crossV,CorPearson,Logl

#' @export

crossV_indices <- function(MC, formulas, modeltype, saveAlea, Y_data_sample,
                           seqthd, resParam_save, Y_data_sample_lcc, #Species, , libfile = .libPaths()
                           MaxDist, Phi, model, lambda) { #nb,

  # For box-cox or Log transformation, data should be strictly positive Add +1 to
  # data if needed
  Sup <- 0
  if ((!is.na(lambda) & !identical(lambda, 1)) |
      grepl("Log|Gamma", modeltype)) {
    if (min(Y_data_sample$dataY) < 0) {
      stop("For Box-Cox or Log transformation or Gamma model, data should be strictly positive")
    }
    if (min(Y_data_sample$dataY) == 0) {
      Sup <- 1
      warning("For Box-Cox or Log transformation or Gamma model, +1 was added to data")
    }
  }
  Y_data_sample$dataY <- Y_data_sample$dataY + Sup
  if (!is.na(Y_data_sample_lcc)) {
    Y_data_sample_lcc$dataY <- Y_data_sample_lcc$dataY + Sup
  }
  # Cut data into 2 random datasets 75-25% (or 90-10%) for cross-validation
  crossV <- na.omit(c(saveAlea[, MC]))  # remove NA when 10-fold crossV is not a multiple of data length
  Y_data_crossV <- Y_data_sample[crossV, ]
  Y_data_adjust <- Y_data_sample[-crossV, ]

  # Null Deviance calculation Deviance is calculated on original data to allow
  # comparison with models without transformation
  Y_data_crossV_dataY <- Y_data_crossV$dataY - Sup

  # data.frame has been modified so that poly is not anymore usable ------------
  Y_data_crossV <- tibble::as_tibble(as.data.frame(Y_data_crossV))
  Y_data_adjust <- tibble::as_tibble(as.data.frame(Y_data_adjust))

  ## Not sure if this is the right calculation but at least no errors with zeros
  # BTW, model selection is based on RMSE
  Zero <- which(Y_data_crossV_dataY == 0)
  Un <- which(Y_data_crossV_dataY != 0)
  DevNull <- rep(NA, length(Y_data_crossV_dataY))

  try(DevNull[Zero] <- -(Y_data_crossV_dataY[Zero] - mean(Y_data_crossV_dataY,
                                                          na.rm = TRUE)))
  try(DevNull[Un] <- Y_data_crossV_dataY[Un] *
        log(Y_data_crossV_dataY[Un]/mean(Y_data_crossV_dataY, na.rm = TRUE)) -
        (Y_data_crossV_dataY[Un] - mean(Y_data_crossV_dataY, na.rm = TRUE)))
  SumDevNull <- 2 * sum(DevNull)
  #   }

  # SSQNull <- (sum((Y_data_crossV_dataY -
  # mean(Y_data_crossV_dataY,na.rm=TRUE))^2))/length(Y_data_crossV_dataY) #
  # deviance(lm(Y_data_crossV_dataY ~ 1))

  y <- 1:nrow(formulas)

  # Apply is for fitting each model one by one ----
  resTOT <- apply(t(y), 2, function(x) {
    formula_tested <- as.character(formulas[x, 1])

    resAIC <- NA
    resUBRE <- NA
    resDev <- NA
    ResDev_crossV <- NA
    predcrossV <- rep(NA, nrow(Y_data_crossV))
    proba <- rep(NA, nrow(Y_data_crossV))
    MeanPercentError <- NA
    # modelX <- NA
    minUn <- NA
    maxZero <- NA
    MeanTHD <- NA
    DiffSelSpe <- rep(NA, length(seqthd))
    ROC_crossV <- NA
    RMSE_crossV <- NA
    CorPearson <- NA
    Logl <- NA

    if (exists("resParam_save")) {
      if (!is.null(resParam_save)) {
        fixXI <- resParam_save[x] #[[nb]][x]
        lambda <- resParam_save[x] #[[nb]][x]
      }
    }

    # Fit model ----
    model_output <- fit_model(x = Y_data_adjust, model = formula_tested,
                              fixXI = fixXI,
                              modeltype = modeltype,
                              x_lcc = Y_data_sample_lcc[-crossV, ],
                              MaxDist = MaxDist, Phi = Phi, Model = model,
                              fix.lambda = fix.lambda,
                              lambda = lambda)
    modelX <- model_output$modelX

    if (!is.na(modelX)[1]) {
      try(resAIC <- AIC(modelX), silent = TRUE)
      # In case of LogCont, AIC is increased so that it is comparable to the AIC of
      # Normal distributions
      if (grepl("TweedGLM", modeltype)) {
        try(resAIC <- tweedie::AICtweedie(modelX))
      }

      if (grepl("Log", modeltype)) {
        try(resAIC <- AIC(modelX) + 2 * sum(log(Y_data_adjust$dataY)),
            silent = TRUE)
      }
      if (!is.na(lambda)) {
        try(resAIC <- AIC(modelX) + 2 *
              sum(log((Y_data_adjust$dataY)^(resParam_save[x] - 1))), ##[[nb]][x]
            silent = TRUE)
      }

      # UBRE is only calculated for GAM models
      if (!grepl("GLM", modeltype)) {
        try(resUBRE <- modelX$gcv.ubre, silent = TRUE)
      }
      if (!grepl("KrigeGLM", modeltype)) {
        try(resDev <- (modelX$null.deviance - modelX$deviance)/modelX$null.deviance,
            silent = TRUE)
      }

      # Prediction for the crossV dataset
      if (grepl("GLM", modeltype)) {
        if (!grepl("KrigeGLM", modeltype)) {
          try(predcrossV <- predict.glm(modelX, newdata = Y_data_crossV,
                                        type = "response"), silent = TRUE)
        } else {
          trend.sim <- trend.pred <- NA
          datageo_fit <- geoR::as.geodata(
            Y_data_sample_lcc[-crossV, ],
            data.col = which(names(Y_data_sample_lcc) == "dataY"),
            covar.col = 4:(ncol(Y_data_sample_lcc)))
          # Modification to make it work with poly() ---------------------------------
          datageo_fit$covariate <- tibble::as_tibble(datageo_fit$covariate)
          datageo_valid <- geoR::as.geodata(
            Y_data_sample_lcc[crossV, ],
            data.col = which(names(Y_data_sample_lcc) == "dataY"),
            covar.col = 4:(ncol(Y_data_sample_lcc)))
          # Modification to make it work with poly() ---------------------------------
          datageo_valid$covariate <- tibble::as_tibble(datageo_valid$covariate)

          try(trend.sim <- geoR::trend.spatial(as.formula(formula_tested), datageo_fit), silent = TRUE)
          try(trend.pred <- geoR::trend.spatial(as.formula(formula_tested), datageo_valid), silent = TRUE)
          if (grepl("KrigeGLM", modeltype) & !grepl("KrigeGLM.dist", modeltype)) {
            try(predcrossV <- geoR::krige.conv(
              geodata = datageo_fit,
              locations = rbind(datageo_valid$coords),
              krige = geoR::krige.control(trend.d = trend.sim, trend.l = trend.pred,
                                          obj.m = modelX))$predict, silent = TRUE)
          }
        }
      } else {
        try(predcrossV <- mgcv::predict.gam(modelX, newdata = Y_data_crossV, type = "response"),
            silent = TRUE)
      }


      # Correct for the data transformation
      predcrossV <- predcrossV - Sup
      if (grepl("KrigeGLM", modeltype)) {
        predcrossV[predcrossV < 0] <- 0
      }

      #---------- Deviance
      if (grepl("Log", modeltype)) {
        # Correct predictions
        try(predcrossV2 <- exp(predcrossV + 0.5 * var(modelX$residuals)),
            silent = TRUE)
      } else {
        predcrossV2 <- predcrossV
      }

      # Residual deviance of the model
      Zero <- which(Y_data_crossV_dataY == 0)
      Un <- which(Y_data_crossV_dataY != 0)
      Dev <- rep(NA, length(Y_data_crossV_dataY))
      # try(Dev <- Y_data_crossV_dataY * log(Y_data_crossV_dataY/predcrossV) -
      # (Y_data_crossV_dataY- predcrossV))
      try(Dev[Zero] <- -(Y_data_crossV_dataY[Zero] - predcrossV2[Zero]))
      try(Dev[Un] <- Y_data_crossV_dataY[Un] * log(Y_data_crossV_dataY[Un]/predcrossV2[Un]) -
            (Y_data_crossV_dataY[Un] - predcrossV2[Un]))

      try(ResDev_crossV <- 100 * 2 * sum(Dev)/SumDevNull)  #2* sum(Dev)
      #         }

      if (grepl("PA", modeltype)) {
        proba <- predcrossV2
      }
      if (grepl("TweedGLM", modeltype)) {
        # keep errors in 'a' so that no error message is shown in log file
        if (sum(is.na(predcrossV2) == FALSE) == length(predcrossV2)) {
          a <- try(Logl <- sum(apply(t(1:length(Y_data_crossV_dataY)),
                                     2, function(i) {
                                       tweedie::dtweedie.logl(
                                         phi = summary(modelX)$dispersion,
                                         y = Y_data_crossV_dataY[i],
                                         mu = predcrossV2[i],
                                         power = resParam_save[x])
                                     })), silent = TRUE)
          # 1-proba because output p is proba of absence
          a <- try(proba <- 1 - (apply(t(predcrossV2), 2,
                                       function(i) tweedie::dtweedie(
                                         y = 0, xi = resParam_save[x],
                                         mu = i, phi = summary(modelX)$dispersion))))
        }
      }

      if (grepl("PA|TweedGLM", modeltype)) {
        # Minimum predicted probability when presence in data If not all are NA
        if (sum(is.na(proba[Un])) != length(proba[Un])) {
          try(minUn <- min(proba[Un], na.rm = TRUE))
        }
        # Maximum predicted probability when absence in data If not all are NA
        if (sum(is.na(proba[Zero])) != length(proba[Zero])) {
          try(maxZero <- max(proba[Zero], na.rm = TRUE))
        }

        # Amount of errors depending on threshold value
        data01 <- (Y_data_crossV_dataY != 0) * 1
        # for (s in 1:length(seqthd)) {
        #   pred01 <- rep(0, length(proba))
        #   try(pred01[which(proba >= seqthd[s])] <- 1, silent = TRUE)
        #   try(DiffSelSpe[s] <- sum(abs(pred01 - data01)))
        # }
        # Difference Selectivity-Specificity according to threshold value
        try({
          roc1 <- ROCR::prediction(c(proba), c(data01))
          perf1 <- ROCR::performance(roc1, "sens", "spec")
          x <- perf1@alpha.values[[1]]
          y <- abs(perf1@x.values[[1]] - perf1@y.values[[1]])
          dat <- data.frame(cutoff = x, diff = y)
          dat <- dat[is.finite(dat$cutoff),]
        }, silent = TRUE)

        try(DiffSelSpe <- approx(x = dat$cutoff, y = dat$diff,
                                 xout = seqthd)$y)

        try(MeanTHD <- mean(seqthd[which(DiffSelSpe == min(DiffSelSpe))][1]))

        # ROC calculation
        try(ROC_crossV <- ROCR::performance(roc1, "auc")@y.values[[1]], silent = TRUE)
      }


      # ---------- Calculation of Indicators of goodness of fit
      if (!grepl("PA", modeltype)) {
        # Residual deviance of the model
        #           try(ResDev_crossV <- -2 * sum(log(Y_data_crossV_dataY/predcrossV) -
        #           ((Y_data_crossV_dataY - predcrossV)/predcrossV))/SumDevNull)
        # CV of error of prediction
        try(MeanPercentError <- mean(abs(Y_data_crossV_dataY - predcrossV2)/Y_data_crossV_dataY,
                                     na.rm = TRUE))
        # Correlation between prediction and observation
        try(CorPearson <- cor(Y_data_crossV_dataY, predcrossV2, method = "pearson"))
      }
      # RMSE calculation (Recommended by Stanford course)
      # Also work for Lognormal distribution as we are trying to predict the mean
      try(RMSE_crossV <- sqrt(sum((Y_data_crossV_dataY - predcrossV2)^2)/length(Y_data_crossV_dataY)))  # / SSQNull)
    }  # end of is.na(modelX)

    if (grepl("PA|TweedGLM", modeltype)) {
      res <- c(resAIC, resUBRE, resDev, ResDev_crossV, minUn, maxZero, MeanTHD, DiffSelSpe,
               ROC_crossV, RMSE_crossV, Logl)
    } else {
      res <- c(resAIC, resUBRE, resDev, ResDev_crossV, MeanPercentError, RMSE_crossV,
               CorPearson, Logl)
    }
    if (grepl("PA|TweedGLM", modeltype)) {
      names(res) <- c("resAIC", "resUBRE", "resDev", "ResDev_crossV", "minUn",
                      "maxZero", "MeanTHD", paste0("DiffSelSpe", 1:length(DiffSelSpe)),
                      "ROC_crossV", "RMSE_crossV", "Logl")
    } else {
      names(res) <- c("resAIC", "resUBRE", "resDev", "ResDev_crossV", "MeanPercentError",
                      "RMSE_crossV", "CorPearson", "Logl")
    }
    res
  })  # end of apply

  #   resTOT <- as.data.frame(resTOT)
  return(resTOT)
}

#' Mean difference between a distribution and a set of others
#' with or without weights
#'
#' @param x matrix with distribution in rows
#' @param n Number of the column to compare with
#' @param w vector of weights with length = ncol(x)
#' @param comp.n Logical Whether to output x-x[n,] matrix
#' @export

meandiff_distri <- function(x, n, w, comp.n = TRUE)
{
  Comp.n <- t(apply(x, 1,
                    function(y) {y - x[n, ]}))
  if (missing(w)) {
    rM <- rowMeans(Comp.n, na.rm = TRUE)
  } else {
    rM <- apply(Comp.n, 1, function(y) {
      y.na <- which(is.na(y))
      if (length(y.na) > 0) {
        res <- sum(y[-y.na] * w[-y.na]/sum(w[-y.na]))
      } else {
        res <- sum(y * w/sum(w))
      }
      return(res)
      # weighted.mean(y, w, na.rm = TRUE)
    })
  }
  if (comp.n) {
    return(list(comp.n = Comp.n, means = rM))
  } else {
    return(list(means = rM))
  }
}


#' A function to rank distributions (of same length) and statistically compare
#' them to the best one
#'
#' @param x Typically a matrix where rows are different distributions of the
#' same length to be compared while paired
#' @param w vector of weights with the same length a ncol(x) if outputs do not
#' have the same weight. Used for weighted.mean and for p-value calculation.
#' @param test test used to compare distribution as used by
#' \code{\link[survey]{svyranktest}}
#' @param na.max proportion maximum of NA value allowed in one distribution.
#' If proportion of NA is upper na.max, model is ranked at the end
#' and no p-value is calculated
#' @param p.min minimum p-value under which the order of distribution is not
#' important because following distributions will not be kept...
#' If set, when p-value is lower than p.min, distributions are supposed
#' significantly "worse" than the best one. Remaining distributions are ordered
#' according to their mean and p-values are not calculated.
#' @param silent Logical Whether to show \% remained or not
#' @param cl a cluster as made with \code{\link[snow]{makeCluster}}.
#' If empty, nbclust in \code{\link{modelselect_opt}} will be used.
#'
#' @return
#' orderModels: number of columns of x re-ordered from best to worse
#' p.values: p-values of difference between all distributions and the best one
#' ordered like orderModels
#' p.min.test: Logical. FALSE if distribution is ordered after the first
#' distribution occurring with a p.value lower than p.min.
#' Indeed, large distributions with high outliers may be not significantly
#' different than distribution 1.
#'
#' @details This function has been developed to compare indices of goodness of
#' fit calculated after a cross-validation procedure.
#' The best distribution is the one being the best on average for all cross-
#' validation sub-samples.
#' The best average hides extreme values that may be
#' due to particular crossV samples (chosen randomly).
#' Distribution are then compared statistically to the best one with paired test.
#' Because the k-fold may return folds with different lengths, the weight of each
#' fold may be corrected with the w parameter.
#' \itemize{
#'   \item Because values compared do not necessarily follow a normal distribution,
#'   t.test is not the best mean comparison test. Wilcoxon do not require
#'   normality and is thus more appropriate here.
#'   \item Size of validation set is not equal, in particular if there are
#'   factor covariates. wilcoxon.test is thus weighted accordingly
#'   \item the only weighted wilcoxon test is from library(survey).
#'   See \code{\link[survey]{svyranktest}}
#' }
#' @export

best_distri <- function(x, w,
                        test = c("wilcoxon", "vanderWaerden", "median","KruskalWallis"),
                        na.max = 0.5, p.min = 0.01, silent = TRUE, cl = NULL)
{

  if (is.null(cl)) {
    cl_inside <- TRUE
  } else {
    cl_inside <- FALSE
  }
  nbclust <- modelselect_opt$nbclust

  x_n <- x_n.xNA <- 1:nrow(x)
  x.NA.count <- apply(x, 1, function(y) {sum(is.na(y))/length(y)})

  x.NA <- which(x.NA.count > na.max)
  if (length(x.NA) == nrow(x)) {
    stop("na.max is too restrictive, there are no distribution left")
  }
  if (length(x.NA) > 0) {
    x <- x[-x.NA, ]
    x_n <- x_n[-x.NA]
    x_n.xNA <- 1:nrow(x)
  }

  if (missing(w)) {
    x_mean <- rowMeans(x, na.rm = TRUE)
  } else {
    # x_mean <- apply(x, 1, function(y) weighted.mean(y, w, na.rm = TRUE))
    x_mean <- apply(x, 1, function(y) {
      y.na <- which(is.na(y))
      if (length(y.na) > 0) {
        res <- sum(y[-y.na] * w[-y.na]/sum(w[-y.na]))
      } else {
        res <- sum(y * w/sum(w))
      }
      return(res)
    })
  }

  orderModels.xNA <- orderModels <- numeric(0)
  ttest <- logical(0)
  count <- 0
  for (orderN in 1:length(x_n.xNA)) { # From best to less best
    if (!silent & orderN %in% round(seq(1, nrow(x), length.out = 10))) {
      print(paste0(count * 10, "%"))
      count <- count + 1
    }
    if (orderN == 1) {
      meanLineMin.tmp <- meanLineMin <- order(x_mean)[1]
      x_tmp <- x
      x_n_tmp <- x_n.xNA
    } else {
      meanLineMin.tmp <- meanLineMin <- order(x_mean[-orderModels.xNA])[1]
      x_tmp <- x[-orderModels.xNA,]
      x_n_tmp <- x_n.xNA[-orderModels.xNA]
    }
    # The last one is added directly
    if (orderN != nrow(x)) {
      maxit <- 30
      for (it in 1:maxit) {
        # Diff should be upper than zero, otherwise it is not the best
        rM <- meandiff_distri(x = x_tmp, n = meanLineMin, w = w, comp.n = FALSE)

        if (order(rM)[1] != meanLineMin) {
          meanLineMin <- order(rM$means)[1]
          # Because there are NA values, we can loop on the same x best
          # distribution without being able to choose.
          # Thus, if a distribution appears twice in the loop,
          # it is considered the best...
          if (meanLineMin %in% meanLineMin.tmp) {
            break
          }
          meanLineMin.tmp <- c(meanLineMin.tmp, meanLineMin)
        } else {
          break
        }
        if (it == maxit) {
          stop(paste("maxit to find the best model has been reached, increase maxit",
                     "or choose a more restrictive na.max"))
        }
      }
      rm(meanLineMin.tmp)
    }
    orderModels.xNA <- c(orderModels.xNA, x_n_tmp[meanLineMin])

    if (orderN == 1) {
      # Calculate all p-value against best model at step 1.
      if (missing(w)) {w <- NULL}
      if (cl_inside) {
        cl <- parallel::makePSOCKcluster(nbclust)
      }

      # Test for the significance of difference with the first model
      # Need to cheat to compare to zero
      Signif.test  <- function(i, xt, w, orderModels.xNA, test) {
        data.tmp <- data.frame(
        val = c(rnorm(500, 10, 1), rep(0, 500)),
        group = rep(1:2, each = 500),
        w = rep(1, 1000)
        )

        if (i != orderModels.xNA[1]) {
          data.tmp <- data.frame(val = c(xt[i, ] - xt[orderModels.xNA[1], ],
                                         rep(0, ncol(xt))),
                                 group = rep(1:2, each = ncol(xt)),
                                 w = rep(w, 2))

          if (is.null(w)) {
            design <- survey::svydesign(ids = ~0, data = data.tmp[,1:2])
          } else {
            design <- survey::svydesign(ids = ~0, data = data.tmp[,1:2],
                                        weights = c(data.tmp[,3]))
          }
          ttest <- survey::svyranktest(formula = val ~ group, design = design,
                                       test = test)$p.value
        } else {
          # ttest = 0 if distributions are equal
          ttest <- 1
        }
        ttest
      }

      ttest.tmp <- snow::parCapply(
        cl, t(1:nrow(x)),
        function(i, xt = x, orderModels.xNA, w, test, Signif.test)
          Signif.test(i = i, xt = xt, orderModels.xNA = orderModels.xNA, w = w, test = test),
        xt = x, orderModels.xNA = orderModels.xNA, w = w,
        test = test, Signif.test = Signif.test)

      if (cl_inside) {
        parallel::stopCluster(cl)
      }
    }

    # if distributions are not statistically equal, no nead to continue ordering
    if (ttest.tmp[orderModels.xNA[orderN]] < p.min) {
      p.min.test <- c(rep(TRUE, orderN - 1), rep(FALSE, nrow(x) - orderN + 1))
      orderModels.xNA <- c(orderModels.xNA, c(1:nrow(x))[-orderModels.xNA][order(x_mean[-orderModels.xNA])])
      ttest <- ttest.tmp[orderModels.xNA]

      orderModels <- x_n[orderModels.xNA]
      if (!silent) {
        print("100%")
      }
      break
    }
    if (orderN == nrow(x)) {
      p.min.test <- rep(TRUE, nrow(x))
      ttest <- ttest.tmp[orderModels.xNA]
      orderModels <- x_n[orderModels.xNA]
    }
  } # end of orderN
  return(list(orderModels = c(orderModels, x.NA),
              p.values = c(ttest, rep(NA, length(x.NA))),
              p.min.test = c(p.min.test, rep(NA, length(x.NA)))
  ))
}
