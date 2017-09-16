#' @title Order fitted models
#' @description Find the models that outperform the others in terms of cross-validation index
#'
#' @param saveWD directory where outputs of the cross-validation procedure have
#' been saved or zip.file of this directory.
#' All necessary information have been saved in this directory to be
#' able to compile results. This is the only output of findBestModel in the
#' global R environment.
#' @param plot logical. Either to create the output plots or not. Saved
#' directly in saveWD.
#' @param zip.file TRUE (default) to save all outputs in a zipfile with saveWD,
#' FALSE for no zipfile output, or path to a new zip file
#' @param cl a cluster as made with \code{\link[snow]{makeCluster}}.
#' If empty, nbclust in \code{\link{modelselect_opt}} will be used.
#'
#' @importFrom grDevices heat.colors png dev.off
#' @import graphics
#'
#' @details
#' According to the family of model, the index of cross-validation procedure is
#' selected to rank all models fitted for all modeltypes conjointly. Model ranks
#' are then statistically compared. Model with lowest ranks not statistically
#' different than the first ranked model are kept.
#'
#' @return
#' VeryBestModels_crossV: Best models among all
#' BestForModeltypes:  Best models among all in the modeltypes type
#' Figures are produced to visualize ranking and eventually help for a choice
#'
#' @export

ModelOrder <- function(saveWD, plot, zip.file = TRUE, cl = NULL)
{

  saveWD <- normalizePath(saveWD)

  if (utils::file_test("-f", saveWD) & file.exists(saveWD) & grepl("zip", saveWD)) {
    if (zip.file == TRUE) {zip.file <- saveWD}
    utils::unzip(saveWD, exdir = gsub(".zip", "", saveWD))
    saveWD <- gsub(".zip", "", saveWD)
  } else if (utils::file_test("-d", saveWD)) {
    if (zip.file == TRUE) {zip.file <- paste0(saveWD, ".zip")}
  } else {
    stop("saveWD is neither an existing directory nor a zip file")
  }

  imgprefix <- paste0(basename(saveWD), "-")
  saveWD <- normalizePath(saveWD)

  lim_pvalue_final <- modelselect_opt$lim_pvalue_final
  # Load variables of model configuration that need to be equivalent as fit
  opt_saved <- readr::read_rds(paste0(saveWD, "/modelselect_opt.save.rds"))
  seqthd <- opt_saved$seqthd

  # Load information
  load(paste0(saveWD, "/Allinfo_all.RData"))
  datatype <- Allinfo_all$datatype
  modeltypes <- Allinfo_all$modeltypes

  # message("Image outputs will be saved in ", saveWD)

  # Combine results of the different modeltypes ----
  # Verif random sub-dataset for crossV used are the same
  allM <- paste0("M", 1:length(modeltypes))
  for (MX in 1:length(allM)) {
    load(file = paste0(saveWD, "/saveAlea", allM[MX], ".RData"))
    if (MX == 1) {
      saveAleaM1 <- saveAlea
    }
    if (MX >= 2) {
      if (!identical(saveAleaM1, saveAlea)) {
        stop(c("Models of ", modeltypes[1], "and", modeltypes[MX],
               "have not been fitted with the same cross-validation random subdatasets"))
      }
    }
  }
  load(file = paste0(saveWD, "/saveAlea.nbM", 1, ".RData"))

  # Combine all results
  resAIC_unlist <- logical(0)
  test_models_unlist <- logical(0)
  ToutResDev_crossV <- logical(0)
  ToutAUC_crossV <- logical(0)
  ToutRMSE_crossV <- logical(0)
  nb_models_per_nb <- logical(0)
  nb_per_modeltypes <- logical(0)
  BestTHD_crossV_unlist <- logical(0)
  DiffSelSpeTHD_crossV_unlist <- logical(0)
  MinUn_crossV <- logical(0)
  MaxZero_crossV <- logical(0)
  MeanPercentError_crossV <- logical(0)
  l_Max_Nb_Var <- logical(0)
  ToutLogl_crossV <- logical(0)
  resParam_unlist <- logical(0)

  allM <- paste0("M", 1:length(modeltypes))
  for (MX in 1:length(allM)) {
    load(file = paste0(saveWD, "/resAIC_save", allM[MX], ".RData"))
    resAIC_unlist <- c(resAIC_unlist, unlist(resAIC_save))
    nb_per_modeltypes <- c(nb_per_modeltypes, length(unlist(resAIC_save)))
    if (grepl("TweedGLM|KrigeGLM", modeltypes[MX])) {
      load(file = paste0(saveWD, "/resParam_save", allM[MX], ".RData"))
      resParam_unlist <- c(resParam_unlist, unlist(resParam_save))
    } else {
      resParam_unlist <- c(resParam_unlist, rep(NA, length(unlist(resAIC_save))))
    }
    load(file = paste0(saveWD, "/test_models", allM[MX], ".RData"))
    test_models_unlist_A <- unlist(test_models)
    test_models_unlist <- c(test_models_unlist, test_models_unlist_A)
    load(file = paste0(saveWD, "/Output_models", allM[MX], ".RData"))
    ToutResDev_crossV_A <- apply(t(1:length(Output_models)), 2, function(x)
      Output_models[[x]][4,,])
    ToutResDev_crossV_A2 <- do.call(rbind, ToutResDev_crossV_A)
    ToutResDev_crossV <- rbind(ToutResDev_crossV, ToutResDev_crossV_A2)
    l_Max_Nb_Var <- c(l_Max_Nb_Var, length(ToutResDev_crossV_A))
    if (grepl("PA|TweedGLM", modeltypes[MX])) {
      ToutAUC_crossV_A <- apply(t(1:length(Output_models)), 2, function(x)
        Output_models[[x]][7 + length(seqthd) + 1,,])
      ToutAUC_crossV_A2 <- do.call(rbind, ToutAUC_crossV_A)
      ToutAUC_crossV <- rbind(ToutAUC_crossV, ToutAUC_crossV_A2)
    }
    if (grepl("Cont|Count|Gamma", modeltypes[MX])) {
      MeanPercentError_crossV_A <- apply(t(1:length(Output_models)),
                                         2, function(x) Output_models[[x]][5,,])
      MeanPercentError_crossV_A2 <- do.call(rbind, MeanPercentError_crossV_A)
      MeanPercentError_crossV <- rbind(MeanPercentError_crossV,
                                       MeanPercentError_crossV_A2)
    }
    if (grepl("PA|TweedGLM", modeltypes[MX])) {
      N_RMSE <- 7 + length(seqthd) + 2
    } else {
      N_RMSE <- 6
    }
    ToutRMSE_crossV_A <- apply(t(1:length(Output_models)), 2, function(x)
      Output_models[[x]][N_RMSE,,])
    ToutRMSE_crossV_A2 <- do.call(rbind, ToutRMSE_crossV_A)
    ToutRMSE_crossV <- rbind(ToutRMSE_crossV, ToutRMSE_crossV_A2)
    nb_models_per_nb <- c(nb_models_per_nb, unlist(lapply(ToutRMSE_crossV_A, nrow)))

    if (grepl("PA|TweedGLM", modeltypes[MX])) {
      load(file = paste0(saveWD, "/BestTHD_crossV", allM[MX], ".RData"))
      load(file = paste0(saveWD, "/DiffSelSpeTHD_crossV", allM[MX], ".RData"))
      BestTHD_crossV_unlist <- c(BestTHD_crossV_unlist, unlist(BestTHD_crossV))
      DiffSelSpeTHD_crossV_unlist_A <- do.call(rbind, DiffSelSpeTHD_crossV)
      DiffSelSpeTHD_crossV_unlist <- rbind(DiffSelSpeTHD_crossV_unlist,
                                           100 * DiffSelSpeTHD_crossV_unlist_A/dim(saveAlea)[1])

      MinUn_crossV_A <- apply(t(1:length(Output_models)), 2, function(x)
        Output_models[[x]][5,,])
      MinUn_crossV_A2 <- do.call(rbind, MinUn_crossV_A)
      MinUn_crossV <- rbind(MinUn_crossV, MinUn_crossV_A2)
      MaxZero_crossV_A <- apply(t(1:length(Output_models)), 2, function(x)
        Output_models[[x]][6,,])
      MaxZero_crossV_A2 <- do.call(rbind, MaxZero_crossV_A)
      MaxZero_crossV <- rbind(MaxZero_crossV, MaxZero_crossV_A2)
    }

    save(BestTHD_crossV_unlist, file = paste0(saveWD, "/BestTHD_crossV_unlist.RData"))
    save(MinUn_crossV, file = paste0(saveWD, "/MinUn_crossV.RData"))
    save(MaxZero_crossV, file = paste0(saveWD, "/MaxZero_crossV.RData"))

    if (grepl("TweedGLM", modeltypes[MX])) {
      ToutLogl_crossV_A <- apply(t(1:length(Output_models)), 2, function(x)
        Output_models[[x]][7 + length(seqthd) + 3,,])
      ToutLogl_crossV_A2 <- do.call(rbind, ToutLogl_crossV_A)
      ToutLogl_crossV <- rbind(ToutLogl_crossV, ToutLogl_crossV_A2)
    }
  }  # END of allMX

  if (grepl("Krige", datatype)) {
    allmodeltypes <- rep(paste(modeltypes, Model, sep = "_"), nb_per_modeltypes)
  } else {
    allmodeltypes <- rep(modeltypes, nb_per_modeltypes)
  }


  # Fig for the choice of the best thd for 'PA' pred ----
  if (plot & grepl("PA", datatype))
  {
    png(filename = paste0(saveWD, "/", imgprefix, "Thresholds_for_prediction.png"),
        width = 18, height = 8, units = "cm", pointsize = 7, res = 300)
    # x11()
    par(mfrow = c(1, 3))
    MinUn_crossV[MinUn_crossV == "Inf"] <- 1
    for (i in c("BestTHD_crossV_unlist", "MinUn_crossV", "MaxZero_crossV")) {
      if (i == "BestTHD_crossV_unlist") {
        histSeuil <- hist(BestTHD_crossV_unlist, main = i)
      } else {
        histSeuil <- hist(apply(get(i), 1,
                                          function(x) mean(x, na.rm = TRUE)),
                                    main = i)
      }
      abline(v = 0.5, col = "blue", lwd = 2)
      SeuilFromHist_crossV <- histSeuil$mids[which(histSeuil$density ==
                                                     max(histSeuil$density))]
      abline(v = SeuilFromHist_crossV, col = "red", lwd = 2)
      axis(1, line = 1, at = SeuilFromHist_crossV, col = "red",
           col.axis = "red")
    }
    dev.off()
  }  # end of datatype == 'PA'


  # Choose the very best models ----
  # Define indice depending on datatype
  minus_index <- 1
  if (grepl("PA", datatype)) {
    IndexForRank_crossV <- -ToutAUC_crossV
    minus_index <- -1
    nameIndex <- "Mean_AUC_crossV"
  } else if (grepl("Cont", datatype)) {
    # IndexForRank_crossV <- log(MeanPercentError_crossV)
    IndexForRank_crossV <- ToutRMSE_crossV
    nameIndex <- "Mean_RMSE_crossV"
  } else if (grepl("Count|Krige", datatype)) {
    # IndexForRank_crossV <- log(ToutRMSE_crossV)
    IndexForRank_crossV <- ToutRMSE_crossV
    nameIndex <- "Mean_RMSE_crossV"
  } else if (grepl("Tweed", datatype)) {
    # IndexForRank_crossV <- -ToutLogl_crossV
    IndexForRank_crossV <- ToutRMSE_crossV
    nameIndex <- "Mean_RMSE_crossV"
  }  # log likelihood

  if (length(which(IndexForRank_crossV == "Inf")) > 0) {
    IndexForRank_crossV[which(IndexForRank_crossV == "Inf")] <- 1e+15
  }
  if (length(which(IndexForRank_crossV == "-Inf")) > 0) {
    IndexForRank_crossV[which(IndexForRank_crossV == "-Inf")] <- -1e+15
  }


  # Rank models ----
  # according to their quality of prediction for each
  # cross-validation sub-dataset
  IndexForRank_crossVMean <- rowMeans(IndexForRank_crossV, na.rm = T)

  ToutRMSE_crossVMean <- apply(ToutRMSE_crossV, 1, function(x) {
    mean(x, na.rm = T)
  })

  # The best model is the one which is the best on average for all crossV samples
  # meanLineMin <- order(OrderIndexForRank_crossVMean)[1]

  # The model that is always smaller than the other is not necessarily the one
  # with the smallest average rank because of the calculation of the rank and
  # the possible same ranking
  # It is also not necessary the one with the smallest average value
  # as comparisons test will be paired and because there are NA values.
  # Although models with too many NA are removed from the analysis
  # To avoid comparison of all distributions against all distributions,
  # we do an iterative procedure in maxit step maximum.
  # This should be better for memory...
  # 1/ The model with the smallest mean is set as the best model
  # 2/ Mean comparison is done for all against this model
  # 3/ If a model shows a mean difference lower than zero, then
  # the lowest is retained and return to 2/

  best.distri <- best_distri(x = IndexForRank_crossV, w = saveAlea.nb,
                             test  = "wilcoxon", na.max = 0.5,
                             p.min = lim_pvalue_final, silent = TRUE,
                             cl = cl)

  orderModels <- best.distri$orderModels
  Line_keep <- best.distri$orderModels[which(best.distri$p.min.test)]

  DiffTo1 <- meandiff_distri(x = IndexForRank_crossV[orderModels,],
                             n = 1, w = saveAlea.nb, comp.n = TRUE)

  # Keep models in orderModels ----
  # Find models having ranks not significantly different (p>=0.001) than the best
  # model lim_pvalue can be large for this step to keep a little more
  # models than necessary
  #   Line_keep <- numeric(0)
  #   orderModels <- order(OrderIndexForRank_crossVMean)
  # ttest <- t.test(IndexForRank_crossV[meanLineMin, ],
  #  IndexForRank_crossV[i,], paired = TRUE)$p.value
  # Size of validation set is not equal, in particular if there are
  # factor covariates. t.test is thus weighted accordingly
  #   for (i in c(1:nrow(IndexForRank_crossV))[orderModels]) {
  #   ttest <- weights::wtd.t.test(x = IndexForRank_crossV[meanLineMin, ] -
  #     IndexForRank_crossV[i,], y = 0, weight = saveAlea.nb,
  #     mean1 = TRUE)$coefficients["p.value"]
  #     if (is.na(ttest) == TRUE) {
  #       Line_keep <- c(Line_keep, i)
  #     } else {
  #       if (ttest >= lim_pvalue_final) {
  #         Line_keep <- c(Line_keep, i)
  #       } else {
  #         break
  #       }
  #     }
  #   }
  #   if (length(Line_keep) > 30) {
  #     Line_keep <- Line_keep[1:30]
  #   }


  # Figures of ranks on Percent Errors ----
  if (plot)
  {
    # Rank models, Select the best model and
    # the best ones that are not significantly
    # different from the 1st
    if (grepl("PA", datatype)) {
      titre <- "AUC on crossV"
    } else {
      titre <- "Root of Mean Squared Error on crossV"
    }

    png(filename = paste0(saveWD, "/", imgprefix, "Ordered_Models.png"),
        width = 50, height = 8, units = "cm", pointsize = 5, res = 300)
    par(mai = c(0.2, 0.4, 0.25, 0.3))
    plot(-10,-10, main = "Mean difference of index with 1st model",
         ylim = c(quantile(c(DiffTo1$comp.n), probs = 0.005, na.rm = TRUE),
                  quantile(c(DiffTo1$comp.n), probs = 0.995, na.rm = TRUE)),
         pch = 20, xlab = "", ylab = "Mean difference with 1st model",
         xlim = c(0,
                  nrow(IndexForRank_crossV) + 1),
         xaxs = "i")

    res <- apply(t(1:nrow(DiffTo1$comp.n)), 2, function(x) {
      boxplot(DiffTo1$comp.n[x, ], at = x, pch = 20,
              border = c("black", "forestgreen")[1 + best.distri$p.min.test[x]],
              pch = 20, cex = 0.5, yaxt = "n", xaxt = "n",
              outline = TRUE, add = TRUE)
    })
    #
    #         boxplot(t(DiffTo1$comp.n),
    #           border = c("black", "forestgreen")[1 + best.distri$p.min.test], pch = 20, cex = 0.5)

    points(DiffTo1$means, col = "green", pch = 20, cex = 0.5)
    abline(h = 0, col = "grey", lty = "dashed")


    par(new = TRUE)
    plot(best.distri$p.values, pch = 20, col = "red",
         xaxt = "n", yaxt = "n", cex = 0.5, xaxs = "i",
         xlim = c(0, nrow(IndexForRank_crossV) + 1),
         ylim = c(0,1), yaxs = "i", xlab = "", ylab = "")
    axis(4, col.axis = "orangered")
    mtext("p-value", side = 4, line = 2)

    abline(h = c(lim_pvalue_final), col = "orangered", lwd = 2)
    abline(v = max(which(best.distri$p.min.test)) + 0.5, col = "green")
    par(new = FALSE)
    dev.off()


    IndexForRank_crossV.plot <- minus_index * IndexForRank_crossV

    png(filename = paste0(saveWD, "/", imgprefix, "Index_Error_on_crossV.png"),
        width = 50, height = 8, units = "cm", pointsize = 5, res = 300)

    plot(-10,-10, main = titre,
         ylim = c(stats::quantile(c(IndexForRank_crossV.plot), probs = 0.005, na.rm = TRUE),
                  stats::quantile(c(IndexForRank_crossV.plot), probs = 0.995, na.rm = TRUE)),
         pch = 20, xlab = "", ylab = gsub("_", " ", nameIndex),
         xlim = c(0, nrow(IndexForRank_crossV.plot) + 1),
         xaxs = "i")
    axis(1)
    # Using apply accelerate calculation for table of models and
    # reduce use of RAM
    res <- apply(t(1:nrow(IndexForRank_crossV.plot)), 2, function(x) {
      boxplot(IndexForRank_crossV.plot[x, ], at = x, pch = 20,
              col =  c("black", "forestgreen")[1 + x %in% Line_keep],
              border = c("black", heat.colors(length(Line_keep)))[
                c(which(Line_keep %in% x) + 1, 1)[1]], yaxt = "n", xaxt = "n",
              outline = FALSE, add = TRUE)
      if (x %in% Line_keep) {
        axis(1, at = x, labels = "", line = 1.5,
             col = heat.colors(length(Line_keep))[which(Line_keep %in% x)],
             col.axis = heat.colors(length(Line_keep))[which(Line_keep %in% x)])
      }
    })
    rm(res)
    gc()
    axis(1, at = Line_keep[1], labels = Line_keep[1], line = 2, col = "red",
         col.axis = "red")

    # Add a line to separate models with 1, 2, 3 ... parameters
    for (i in 1:length(nb_models_per_nb)) {
      abline(v = c(0, sum(nb_models_per_nb[1:i])) + 0.5, col = "blue",
             lwd = 2)
    }
    abline(v = 0.5, col = "red", lwd = 3)
    mtext(modeltypes[1], side = 1, at = 0.5, cex = 1.5)
    for (i in 1:(length(l_Max_Nb_Var) - 1)) {
      pos <- cumsum(l_Max_Nb_Var)[i]
      abline(v = sum(nb_models_per_nb[1:pos]) + 0.5, col = "red",
             lwd = 3)
      mtext(modeltypes[i+1], side = 1,
            at = sum(nb_models_per_nb[1:pos]) + 0.5, cex = 1.5)
    }

    dev.off()

  }  # end of plot

  # Create table with the models ranked ----
  # and with the non-statistically differantiated best models
  if (grepl("RMSE", nameIndex)) {
    VeryBestModels_crossV_A <- data.frame(
      Num = 1:length(test_models_unlist),
      model = test_models_unlist,
      nameIndex = round(IndexForRank_crossVMean, digits = 2),
      AIC = round(resAIC_unlist)
    )
  } else {
    VeryBestModels_crossV_A <- data.frame(
      Num = 1:length(test_models_unlist),
      model = test_models_unlist,
      nameIndex = minus_index * round(IndexForRank_crossVMean, digits = 2),
      Mean_RMSE_crossV = round(ToutRMSE_crossVMean, digits = 2),
      AIC = round(resAIC_unlist)
    )
  }
  # With tweedie, selection is made on SSQ but we need to get AUC of the
  # binomial part
  if (grepl("Tweed", datatype)) {
    VeryBestModels_crossV_A$Mean_AUC_crossV <- -ToutAUC_crossVMean
  }
  # We keep XI fo tweedie or Lambda for KrigeGLM
  if (length(grep("TweedGLM|KrigeGLM", modeltypes)) > 0) {
    VeryBestModels_crossV_A$Param <- resParam_unlist
  }
  VeryBestModels_crossV_A$modeltype <- allmodeltypes

  names(VeryBestModels_crossV_A)[
    which(names(VeryBestModels_crossV_A) == "nameIndex")] <- nameIndex

  for (i in 3:(ncol(VeryBestModels_crossV_A) - 1)) {
    VeryBestModels_crossV_A[, i] <- as.numeric(as.character(
      VeryBestModels_crossV_A[, i]))
  }


  VeryBestModels_crossV <- VeryBestModels_crossV_A[Line_keep, ]
  # Save all models ordered
  AllModels_crossV_ordered <- VeryBestModels_crossV_A[orderModels, ]
  utils::write.table(AllModels_crossV_ordered,
                     file = paste0(saveWD, "/", imgprefix,
                                   "AllModels_crossV_ordered.csv"),
                     sep = ";", col.names = TRUE, row.names = FALSE)

  LimitNb <- nrow(VeryBestModels_crossV)

  # List of covariates ranked ----
  # with higher score when on the top of the list
  nChar_crossV <- logical(0)
  sumTot_crossV <- 0
  for (i in 1:LimitNb) {
    nChar_crossV[i] <- length(strsplit(gsub(", ", ",", as.character(
      VeryBestModels_crossV[, "model"][i])), split = " ")[[1]])
    sumTot_crossV <- sumTot_crossV + i
  }
  freq_fact_crossV2 <- tapply(rep(c(LimitNb:1), nChar_crossV),
                              unlist(strsplit(gsub(", ", ",", as.character(
                                VeryBestModels_crossV[, "model"])),
                                split = " ")), sum)

  ScoreTab_A <- t(t(round(t(t(freq_fact_crossV2[
    order(freq_fact_crossV2, decreasing = TRUE)]))/sumTot_crossV * 100)))
  ScoreTab <- t(t(ScoreTab_A[-which(rownames(ScoreTab_A) %in%
                                      c("+", "~", "dataY", "log(dataY)")), ]))
  utils::write.table(ScoreTab, file = paste0(saveWD, "/", imgprefix, "ScoreTab_VeryBest.csv"),
              sep = ";", col.names = TRUE, row.names = TRUE)

  freq_fact_crossVlist2 <- strsplit(gsub(", ", ",", as.character(
    VeryBestModels_crossV[, "model"])), split = " ")

  # Compare models expression to the first 20 covariates chosen by the
  # best models
  NbParam <- min(20, nrow(ScoreTab))
  ParamMat <- matrix(NA, nrow = NbParam, ncol = nrow(VeryBestModels_crossV))
  for (i in 1:nrow(VeryBestModels_crossV)) {
    VeryBestModels_crossV[i, "Num"]
    for (j in 1:NbParam) {
      if (sum(freq_fact_crossVlist2[[i]] %in% rownames(ScoreTab)[j]) >= 1) {
        ParamMat[j, i] <- 1
      }
    }
  }

  if (plot) {
    png(filename = paste0(saveWD, "/", imgprefix, "Help_for_model_Choice_NbParam-",
                          NbParam, ".png"),
        width = 18, height = 9, units = "cm",
        pointsize = 7, res = 300)
    # x11(h=6,w=12) par(mai=c(0.5,4,0.8,0.1))
    mai_left <- max(strwidth(c(paste0(nameIndex),
                               paste0(rownames(ScoreTab)[NbParam:1], " ")),
                             cex = 1, units = "inches"))
    par(mai = c(0.3, mai_left + 0.1, 0.7, 0.4))
    # par(mai = c(0.4, 3, 0.6, 0.4))
    image(1:ncol(ParamMat), 1:nrow(ParamMat), t(ParamMat[nrow(ParamMat):1,]),
          xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "forestgreen")
    axis(2, at = 1:nrow(ParamMat), rownames(ScoreTab)[NbParam:1], las = 1)
    axis(3, at = 1:ncol(ParamMat), labels = NA, las = 2)
    axis(3, at = 1:ncol(ParamMat)-0.2, VeryBestModels_crossV[, "Num"],
         las = 2, tick = FALSE, lwd = 0)
    axis(3, at = 1:ncol(ParamMat)+0.2, VeryBestModels_crossV[, "modeltype"],
         las = 2, cex.axis = 0.6, tick = FALSE, lwd = 0)
    #       SumScore <- apply(ParamMat, 2, function(x) sum(x, na.rm = TRUE))
    #       for (i in min(SumScore):max(SumScore)) {
    #         axis(1, at = c(1:ncol(ParamMat))[which(SumScore == i)],
    #           SumScore[which(SumScore == i)], tick = FALSE,
    #           col.axis = rev(heat.colors(max(SumScore) - min(SumScore) + 1))[
    #             i - min(SumScore) + 1])
    #       }
    axis(1, at = 1:ncol(ParamMat), labels = VeryBestModels_crossV[,nameIndex],
         tick = FALSE, las = 2, cex.axis = 0.8)
    abline(v = (1:ncol(ParamMat)) + 0.5)
    abline(h = (1:nrow(ParamMat)) + 0.5)
    arrows(0.5, nrow(ParamMat) + (nrow(ParamMat)/4.5), ncol(ParamMat) + 0.5,
           nrow(ParamMat) + (nrow(ParamMat)/4.5), xpd = TRUE,
           angle = 20, length = 0.15)
    text((ncol(ParamMat) + 0.5)/2, nrow(ParamMat) + (nrow(ParamMat)/4),
         "+++   From Best to less best model   +",
         xpd = TRUE, cex = 1.5)
    arrows(ncol(ParamMat) + 0.5 + 0.03 * ncol(ParamMat), nrow(ParamMat) + 0.5,
           ncol(ParamMat) + 0.5 + 0.03 * ncol(ParamMat),
           0.5, xpd = TRUE, angle = 15, length = 0.1)
    text(ncol(ParamMat) + 0.5 + 0.05 * ncol(ParamMat), 0.5 * nrow(ParamMat) + 0.5,
         "+++ From Best to less best factor +",
         xpd = TRUE, cex = 1.5, srt = 270)

    mtext(paste0(nameIndex, ": "), side = 1, line = 1,
          at = 0.5, adj = 1, cex = 1)

    dev.off()
  }  # end of plot

  # Save the table of the very best models
  utils::write.table(VeryBestModels_crossV,
              file = paste0(saveWD, "/", imgprefix,
                            "VeryBestModels_crossV.csv"),
              sep = ";", col.names = TRUE, row.names = FALSE)
  save(VeryBestModels_crossV,
       file = paste0(saveWD, "/", imgprefix, "VeryBestModels_crossV.RData"))

  # Save the 2 best models for each modeltypes
  if (length(modeltypes) > 1) {
    w.bestmodeltype <- NULL
    for (modeltype in levels(as.factor(modeltypes))) {
      w.bestmodeltype <- c(w.bestmodeltype, orderModels[which(
        VeryBestModels_crossV_A[orderModels, "modeltype"] == modeltype)[1:2]])
    }
    BestForModeltypes <- VeryBestModels_crossV_A[w.bestmodeltype, ]
    if (grepl("TweedGLM|KrigeGLM", datatype)) {
      BestForModeltypes <- cbind(BestForModeltypes, resParam_unlist[w.bestmodeltype])
      names(BestForModeltypes)[ncol(BestForModeltypes)] <- "Param"
    }
    utils::write.table(BestForModeltypes, file = paste0(saveWD, "/", imgprefix,
                                                    "BestForModeltypes.txt"),
                col.names = TRUE, row.names = FALSE)
  }
  # Zip all outputs in a unique file.
  if (zip.file != FALSE) {
    # Options j allows to save directly files without path
    utils::zip(zip.file, files = saveWD, flags = "-urj9X", extras = "",
        zip = Sys.getenv("R_ZIPCMD", "zip"))
    # Return the file path of the zip file
    return(list(VeryBestModels_crossV = VeryBestModels_crossV,
                BestForModeltypes = BestForModeltypes,
                zip.file = zip.file))
  } else {
    return(list(VeryBestModels_crossV = VeryBestModels_crossV,
                BestForModeltypes = BestForModeltypes))
  }
}  # end of ModelOrder function

#' Characteristics of the model on which to produce outputs
#'
#' @param saveWD directory where outputs of model selection are saved, including
#' AllModels_crossV_ordered.
#' @param new.data SpatialPointsDataFrame dataset used to fit the model.
#' @param modeltype The model type as chosen among modeltypes available
#' (see \code{\link{modelselect_opt}})
#' @inheritParams ModelResults
#'
#' @return Num, modeltype, Family, formulaX
#' @importFrom splines ns
#' @importFrom mgcv gam
#'
#' @export

model_select <- function(saveWD, new.data, Num = NULL, model = NULL,
                         modeltype = NULL, powerXI = NULL)
{
  Param <- NULL
  # Read list of models
  AllModels_crossV_ordered <- readr::read_delim(
    file = paste0(saveWD, "/", basename(saveWD), "-AllModels_crossV_ordered.csv"),
    delim = ";")#, col_types = "ncnnncn")

  new.data.orig <- new.data
  if (is(new.data)[1] == "SpatialPointsDataFrame") {
    new.data <- dplyr::as.tbl(new.data@data)
  } else {
    new.data <- dplyr::as.tbl(new.data)
  }

  # If Num is not specified, choose the first best model
  if (is.null(Num)) {
    Num <- AllModels_crossV_ordered$Num[1]
  }
  if (!is.numeric(Num)) {Num <- as.numeric(as.character(Num))}
  # Verify if model defined by user exists in already fitted models
  if (!is.null(model)) {
    if (sum(AllModels_crossV_ordered$model == model) == 1) {
      Num <- AllModels_crossV_ordered$Num[which(
        AllModels_crossV_ordered$model == model)]
      Num <- as.numeric(as.character(Num))
      message(c("Your model corresponds to model Num = ", Num))
    } else {
      Num <- "UserModel"
    }
  }

  if (Num != "UserModel") {
    # Usual case with model already fitted in crossV
    modeltype <- as.character(AllModels_crossV_ordered[
      which(AllModels_crossV_ordered$Num == Num), "modeltype"])
    formulaX <- as.character(AllModels_crossV_ordered[
      which(AllModels_crossV_ordered$Num == Num), "model"])
    if (grepl("PA|Tweed", modeltype)) {
      load(file = paste0(saveWD, "/BestTHD_crossV_unlist.RData"))
      load(file = paste0(saveWD, "/MinUn_crossV.RData"))
      load(file = paste0(saveWD, "/MaxZero_crossV.RData"))

      if (grepl("Tweed", modeltype)) {
        # THD values saved only for Tweed model, not for other Cont ones.
        ord.Num <- order(unlist(AllModels_crossV_ordered[,"Num"]))
        Num.thd <- which(unlist(AllModels_crossV_ordered[ord.Num,][
          which(grepl("Tweed", unlist(AllModels_crossV_ordered[ord.Num,"modeltype"]))), "Num"]) == Num)
      } else {
        Num.thd <- Num
      }
      Seuil <- BestTHD_crossV_unlist[Num.thd]
      # Minimum value where we see presence
      SeuilMinUn <- mean(MinUn_crossV[Num.thd, ][
        is.finite(MinUn_crossV[Num.thd, ])], na.rm = TRUE)
      # Maximum value where we see absence
      SeuilMaxZero <- mean(MaxZero_crossV[Num.thd, ][
        is.finite(MaxZero_crossV[Num.thd, ])], na.rm = TRUE)
    }
  } else {
    # Case with User defined model not already in crossV
    # No guaranty on the utility of this type of operation
    formulaX <- model
    if (grepl("PA|Tweed", modeltype)) {
      ## Write the function to calculate Seuil(s) by crossV if model not
      ## already estimated
      Seuil <- 0.5
      SeuilMinUn <- 0.95  # Minimum value where we see presence
      SeuilMaxZero <- 0.05  # Maximum value where we see absence
      warning(c("Presence-absence 'Seuil(s)' values have been fixed for",
                " compatibility purpose but they have not any sense for the",
                " model chosen"))
    }
  }


  # Modify formulaX for factor covariates
  formulaX.tmp <- as.character(formulaX)
  w.NotNum <- which(unlist(lapply(new.data[,all.vars(as.formula(formulaX.tmp))],
                                  function(x) {!is.numeric(x)})))
  if (length(w.NotNum) > 0) {
    for (var in w.NotNum) {
      if (!grepl(paste0("as.factor(as.character(",
                        all.vars(as.formula(formulaX.tmp))[var], "))"),
                 formulaX.tmp)) {
        formulaX <- gsub(all.vars(as.formula(formulaX.tmp))[var],
                         paste0("as.factor(as.character(", all.vars(as.formula(formulaX.tmp))[var], "))"),
                         formulaX.tmp, fixed = TRUE)
      }
    }
  }



  if (grepl("PA", modeltype)) {
    Family <- "binomial"
  } else if (grepl("Cont", modeltype)) {
    Family <- "gaussian"
  } else if (grepl("Gamma", modeltype)) {
    Family <- "Gamma(link=\"log\")"
  } else if (grepl("Count", modeltype)) {
    Family <- "poisson"
  } else if (grepl("Tweed", modeltype)) {
    if (Num != "UserModel") {
      Param <- unlist(AllModels_crossV_ordered[
        which(AllModels_crossV_ordered$Num == Num), "Param"])
    } else {
      Param <- tweedie::tweedie.profile(as.formula(formulaX), data = new.data,
                                        xi.vec = seq(1.1, 1.9, 0.1), do.plot = FALSE, fit.glm = FALSE,
                                        do.ci = FALSE, method = "series")$xi.max[1]
    }
    Family <- statmod::tweedie(var.power = Param, link.power = 0)
  } else if (grepl("KrigeGLM", modeltype)) {
    ## To do
  }

  ## Add data + 1
  ## Add warning + output param to say it is + 1 data
  # For box-cox or Log transformation, data should be strictly positive Add +1 to
  # data if needed
  Sup <- 0
  ## Add for KrigeGLM with lambda
  if (grepl("Log|Gamma", modeltype)) {
    if (min(new.data$dataY) < 0) {
      stop("For Box-Cox or Log transformation or Gamma model, data should be strictly positive")
    }
    if (min(new.data$dataY) == 0) {
      Sup <- 1
      warning("For Box-Cox or Log transformation or Gamma model, +1 was added to data")
    }
  }
  new.data$dataY <- new.data$dataY + Sup
  # For future KrigeGLM
  #  if (!is.na(Y_data_sample_lcc)) {
  #    Y_data_sample_lcc$dataY <- Y_data_sample_lcc$dataY + Sup
  #  }

  if (grepl("GLM", modeltype)) {
    modelX <- glm(as.formula(formulaX), data = new.data, family = Family,
                  control = list(epsilon = 1e-05, maxit = 1000))
  } else {
    modelX <- mgcv::gam(as.formula(formulaX), method = "REML", data = new.data,
                        family = Family)
  }

  if (grepl("Tweed", modeltype)) {
    ## Verify with TweedGAM
    # 1-proba because output p is proba of absence
    modelX$proba <- 1 - (apply(t(modelX$fitted), 2, function(i) tweedie::dtweedie(y = 0,
                                                                         xi = Param, mu = i, phi = summary(modelX)$dispersion)))
  }

  res <- list(Num = Num, modeltype = as.character(modeltype), Family = Family,
              formulaX = formulaX, modelX = modelX, Param = Param,
              new.data = new.data.orig, Sup = Sup
  )

  if (grepl("PA|Tweed", modeltype)) {
    res$Seuil <- Seuil
    res$SeuilMinUn <- SeuilMinUn
    res$SeuilMaxZero <- SeuilMaxZero
  }
  return(res)
}

#' Study the outputs and predictions of a specific model
#'
#' @param saveWD directory where outputs of the cross-validation procedure have
#' been saved or zip.file of this directory.
#' All necessary information have been saved in this directory to be
#' able to compile results. This is the only output of findBestModel in the
#' global R environment.
#' @param plot logical Whether to do plot outputs or not. Although there are not
#' a lot of other types of outputs...
#' @param zip.file TRUE (default) to save all outputs in a zipfile with saveWD,
#' FALSE for no zipfile output, or path to a new zip file
#' @param Num numeric The 'Num' value of the model retained for further analysis
#' like those saved in AllModels_crossV_ordered. Not useful if model and Family
#' are set. Default to the first model in AllModels_crossV_ordered.
#' @param model character string A user-defined model specification. This requires
#' to set the Family parameter. This overrides effects of 'Num', however, if the
#' model exists with the exact same expression in the list of model fitted,
#' 'Num' will be re-attributed.
#' @param modeltype The type of the model as chosen among modeltypes available
##' (see \code{\link{modelselect_opt}})
#' @param powerXI XI parameter for a tweedie model.
#' @param Marginals Calculate full marginal prediction for each covariate
#' @param cl a cluster as made with \code{\link[snow]{makeCluster}}.
#' If cl is empty, nbclust in \code{\link{modelselect_opt}} will be used
#' to create a cluster.
#'
#' @importFrom magrittr "%>%"
#' @importFrom sp spTransform CRS
#' @import graphics
#' @importFrom grDevices dev.off png heat.colors colorRampPalette
#' @importFrom ggplot2 ggplot geom_pointrange aes facet_wrap ggsave
#' @importFrom yarrr pirateplot
#' @import stats
#'
#' @export

ModelResults <- function(saveWD, plot, Num = NULL, model = NULL,
                             modeltype = NULL, powerXI = NULL, Marginals = FALSE,
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
  }

  if (!is.null(Num) & !is.numeric(Num)) {
    Num <- as.numeric(as.character(Num))
  }

  imgprefix <- paste0(basename(saveWD), "-")
  saveWD <- normalizePath(saveWD)

  message("Image outputs will be saved in ", saveWD)

  # Load information
  load(paste0(saveWD, "/Allinfo_all.RData"))
  datatype <- Allinfo_all$datatype

  nbclust <- modelselect_opt$nbclust

  # Load variables of model configuration that need to be equivalent as fit
  opt_saved <- readr::read_rds(paste0(saveWD, "/modelselect_opt.save.rds"))
  # Options ----
  nbMC <- opt_saved$nbMC
  seqthd <- opt_saved$seqthd
  MaxDist <- opt_saved$MaxDist
  Phi <- opt_saved$Phi
  Model <- opt_saved$Model
  Lambda <- opt_saved$Lambda
  modeltypes <- opt_saved$modeltypes
  lcc_proj <- opt_saved$lcc_proj

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
  # Precision to keep for predictions
  prec.data <- prec_data(Y_data_sample$dataY) - 1

  # Define the model on which to produce outputs
  model_selected <- model_select(saveWD = saveWD,
                                 new.data = Y_data_sample,
                                 Num = Num, model = model,
                                 modeltype = modeltype, powerXI = powerXI)

  Num <- model_selected$Num
  modeltype <- model_selected$modeltype
  Family <- model_selected$Family
  formulaX <- model_selected$formulaX
  modelX <- model_selected$modelX
  Param <- model_selected$Param
  if (grepl("PA|Tweed", modeltype)) {
    Seuil <- model_selected$Seuil
    SeuilMinUn <- model_selected$SeuilMinUn
    SeuilMaxZero <- model_selected$SeuilMaxZero
  }

  # Outputs of the best chosen model ====
  # Find the prediction closest to mean average prediction ----
  # This is to build some marginal distributions of covariates
  datapred.tbl <- model_selected$new.data[,all.vars(stats::as.formula(modelX))] %>%
    dplyr::as.tbl() %>%
    dplyr::mutate_at(dplyr::vars(-dataY), .funs = function(x) {
      if (is.numeric(x)) {
        minmax <- stats::quantile(x, probs = c(0.1, 0.9), na.rm = TRUE)
        x[(x <= minmax[1] | x >= minmax[2])] <- NA
      }
      x
    }) %>%
    dplyr::mutate(pred = stats::predict(modelX),
                  na.row = apply(., 1, function(x) sum(is.na(x)) != 0),
                  mean = ifelse(grepl("PA", modeltype), 0, mean(dataY)),
                  pred.absmean = abs(pred - mean)) %>%
    dplyr::filter(!na.row) %>%
    dplyr::arrange(pred.absmean)

  CovForMeanPred <- datapred.tbl[1,]

  readr::write_rds(CovForMeanPred,
                   path = paste0(saveWD, "/", basename(saveWD),
                                 "-CovForMeanPred_", model_selected$Num, ".rds"))

  # Deviance explained ----
  ## Also add AUC and RMSE
  w.dataY <- "dataY"
  if (model_selected$Sup != 0) {
    w.dataY <- "dataY + 1"
  }

  if (grepl("Log", modeltype)) {
    factlist <- c(paste0("log(", w.dataY, ") ~ 1"), gsub("dataY", w.dataY,
                                                         unlist(strsplit(formulaX, "+", fixed = TRUE))))
  } else {
    factlist <- c(paste(w.dataY, "~ 1"), gsub("dataY", w.dataY, unlist(strsplit(
      formulaX,
      "+", fixed = TRUE))))
  }
  aov_save <- as.data.frame(matrix(ncol = 5, nrow = length(factlist)))
  rownames(aov_save) <- factlist
  # For crossV indices
  load(file = paste0(saveWD, "/saveAleaM1.RData"))
  RMSE_crossV_save <- matrix(0, nrow = length(factlist), ncol = 2)
  if (cl_inside) {
    cl <- parallel::makePSOCKcluster(nbclust)
  }
  for (fl in 1:length(factlist)) {
    if (fl == 1) {
      model <- factlist[1]
      if (grepl("GLM", modeltype)) {
        assign(paste0("modelX", fl), stats::glm(as.formula(model),
                                         data = Y_data_sample, family = Family, control = list(epsilon = 1e-06,
                                                                                               maxit = 1000)))
      } else {
        assign(paste0("modelX", fl), mgcv::gam(as.formula(model),
                                         method = "REML", data = Y_data_sample, family = Family))
      }
    } else {
      model <- paste(factlist[2:fl], collapse = "+")
      if (grepl("GLM", modeltype)) {
        assign(paste0("modelX", fl), stats::glm(as.formula(paste(factlist[2:fl],
                                                          collapse = "+")), data = Y_data_sample, family = Family,
                                         control = list(epsilon = 1e-06, maxit = 1000)))
      } else {
        assign(paste0("modelX", fl), mgcv::gam(as.formula(paste(factlist[2:fl],
                                                          collapse = "+")), method = "REML", data = Y_data_sample,
                                         family = Family))
      }
      # Compare anova
      if (fl == 2) {
        aov1 <- anova(get(paste0("modelX", fl - 1)), get(paste0("modelX",
                                                                fl)), test = "Chisq")
        aov_save[1:2, ] <- aov1
        names(aov_save) <- names(aov1)
      } else {
        aov_save[fl, ] <- anova(get(paste0("modelX", fl - 1)),
                                get(paste0("modelX", fl)), test = "Chisq")[2,]
      }
    }
    print(model)
    # Calculate crossV indices ----
    RMSE_crossV <- snow::parApply(
      cl, t(1:nbMC), 2,
      function(MC,
               formula, modeltype, saveAlea, Y_data_sample,
               seqthd, resParam_save, nb, Y_data_sample_lcc,
               MaxDist, Phi, Model, lambda) {
        res <- crossV_indices(
          MC = MC, formulas = t(t(formula)),
          modeltype = modeltype, saveAlea = saveAlea, Y_data_sample = Y_data_sample,
          seqthd = seqthd,
          resParam_save = resParam_save, Y_data_sample_lcc = Y_data_sample_lcc,
          MaxDist = MaxDist, Phi = Phi, model = Model,
          lambda = lambda)
        if (grepl("PA", modeltype)) {
          res <- res["ROC_crossV",]
          # Test for threshold value
          # res <- res[grep("DiffSelSpe", rownames(res)),]
        } else {
          res <- res["RMSE_crossV",]
        }
        return(unlist(res))
        c(res)
      }, formula = gsub(" + 1", "", model, fixed = TRUE),
      modeltype = modeltype,
      saveAlea = saveAlea, Y_data_sample = Y_data_sample,
      seqthd = seqthd, resParam_save = Param,
      Y_data_sample_lcc = Y_data_sample_lcc, MaxDist = MaxDist, #nb = nb,
      Phi = Phi, Model = Model[which(modeltypes == modeltype)],
      lambda = Lambda[which(modeltypes == modeltype)])

    # Test for threshold value
    # plot(seqthd, RMSE_crossV[,1], type = "l")
    # apply(RMSE_crossV, 2, function(x) lines(seqthd, x))
    # lines(seqthd,
    #       apply(RMSE_crossV, 1, function(x) mean(x, na.rm = TRUE)),
    #       col = "red", lwd = 2)

    RMSE_crossV_save[fl,1] <- mean(RMSE_crossV, na.rm = TRUE)
    if (fl >= 2) {
      RMSE_crossV_save[fl,2] <- RMSE_crossV_save[fl,1] - RMSE_crossV_save[fl - 1,1]
    }
  }
  # Round values of RMSE properly
  diff.rmse <- max(RMSE_crossV_save[,2]) - min(RMSE_crossV_save[fl,2])
  if (grepl("PA", modeltype)) {
    prec <- 2
  } else {
    prec <- 1
  }

  if (diff.rmse > 1) {
    # precision defined by length of the number
    n <- nchar(round(diff.rmse))
  } else {
    # precision defined by number of zero after comma
    n.w0 <- (nchar(diff.rmse) - 2) # nchar with zero after comma
    n.wo0 <- nchar(diff.rmse * 10 ^ n.w0) # nchar without zero
    n <- -n.w0 + n.wo0 # diff of nchar
  }
  RMSE_crossV_save <- round(RMSE_crossV_save / 10 ^ (n - prec)) * 10 ^ (n - prec)
  if (cl_inside) {
    parallel::stopCluster(cl); gc()
  }
  aov_save[,5] <- cut(aov_save[, 5],
                      breaks = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
                      labels = c("p<0.0001", "p<0.001", "p<0.01", "p<0.05", "p<0.1", "NS"))
  names(aov_save)[5] <- "p"
  aov_save <- rbind(aov_save, c(rep(NA, 3),
                                aov_save[nrow(aov_save), 2], NA, NA))

  rownames(aov_save)[nrow(aov_save)] <- "Residuals"
  aov_save <- cbind(aov_save, round(aov_save[, "Deviance"] /
                                      aov_save[1, 2] * 100, digits = 2),
                    rbind(RMSE_crossV_save, NA))
  if (grepl("PA", modeltype)) {
    names(aov_save)[6:8] <- c("%Exp.Dev", "AUC", ">Diff")
  } else {
    names(aov_save)[6:8] <- c("%Exp.Dev", "RMSE", ">Diff")
  }


  save(aov_save, file = paste0(saveWD, "/aov_save_", Num, ".RData"))
  utils::write.csv(aov_save,
            file = paste0(saveWD, "/", imgprefix, "aov_save_", Num, ".csv"))


  if (plot) {
    # Fig Deviance ----
    png(filename =
          paste0(saveWD, "/", imgprefix, "Param_Exp_Deviance_", Num, ".png"),
        width = 11, height = 5, units = "cm", pointsize = 8,
        res = 300)
    # Calculate number of character to define par(mai) on left
    mai_left <- max(strwidth(paste0(rownames(aov_save), " "), cex = 0.5, units = "inches"))
    par(mai = c(0, mai_left + 0.05, 0.3, 0.05), cex = 0.5)
    plot(NA, xlim = c(0.75, ncol(aov_save) + 0.25), ylim = c(1, nrow(aov_save)),
         yaxt = "n", xaxt = "n", xlab = "", ylab = "",
         main = "Contribution of factors")
    axis(2, at = (2:nrow(aov_save)), labels = rev(paste0(rownames(aov_save), " "))[-1],
         las = 1, tick = FALSE, line = -1)
    axis(2, at = 1, labels = paste0(rownames(aov_save), " ")[nrow(aov_save)],
         las = 1, font = 3, tick = FALSE, line = -1)
    axis(3, at = c(1:ncol(aov_save)), labels = names(aov_save), las = 1,
         tick = FALSE, line = -0.5)
    for (i in 1:ncol(aov_save)) {
      for (j in c(1:nrow(aov_save))) {
        font <- 1
        if (j == nrow(aov_save)) {font <- 3} # Residuals
        if (i %in% c(1:4)) {
          text(i, nrow(aov_save) - j + 1,
               labels = round(aov_save[j, i], digits = 2), font = font)
        }
        #         if (i == 5) {
        #           text(i, nrow(aov_save) - j + 1,
        #                labels = round(aov_save[j, i], digits = 4), font = font)
        #         }
        if (i %in% c(5,6)) {
          text(i, nrow(aov_save) - j + 1,
               labels = aov_save[j, i], font = font)
        }
        if (i %in% c(7,8)) {
          text(i, nrow(aov_save) - j + 1,
               labels = aov_save[j, i], font = 3)
        }
      }
    }
    abline(h = c(1.5, nrow(aov_save) - 0.5), lty = "dashed", col = "grey")
    dev.off()

    # Fig Residuals analysis ----
    png(filename = paste0(saveWD, "/", imgprefix, "Residual-Analysis_", Num, ".png"),
        width = 15, height = 15, units = "cm", pointsize = 8,
        res = 300)

    par(mfrow = c(2,2))
    plot(fitted(modelX), residuals(modelX), col = "black",
         xlab = "Fitted", ylab = "Resp. residuals", main = "Resp. residuals vs Fitted")
    panel.smooth(fitted(modelX), residuals(modelX, type = "pearson"), col = NA)
    abline(h = 0, v = 0, lty = 3, col = "gray")

    qqnorm(residuals(modelX), main = "Residuals QQplot")
    qqline(residuals(modelX), lty = 3, col = "gray50")

    plot(fitted(modelX), residuals(modelX, type = "pearson"), col = "black",
         xlab = "Fitted", ylab = "Std. residuals", main = "Std residuals vs Fitted")
    panel.smooth(fitted(modelX), residuals(modelX, type = "pearson"), col = NA)
    abline(h = 0, v = 0, lty = 3, col = "gray")

    plot(cooks.distance(modelX), type = "h", ylab = "Cook's distance",
         main = "Cook's distance")
    abline(h = c(0.5, 1), lty = c(3,1) , col = "gray50")
   dev.off()

  }

  # Covariates and Interactions used ====

  # Find if estimates are for numeric, factors or interactions

  # Get list of covariates in the model
  n_tmp <- unlist(strsplit(gsub(", ", ",", as.character(formulaX)), split = " "))
  cov_used_A <- n_tmp[-which(n_tmp %in% c("~", "+", "dataY", "log(dataY)"))]
  cov_used_B <- apply(t(cov_used_A), 2, function(x) strsplit(x, split = ",")[[1]][1])
  if (length(grep("GLM", modeltype)) == 1) {
    # For GLM
    cov_used_B <- cov_used_A
  }

  # Get factors in interactions
  if (length(grep("te(", cov_used_B, fixed = TRUE)) > 0 |
      length(grep(":", cov_used_B, fixed = TRUE)) > 0) {
    cov_used_C <- t(cov_used_B[-c(grep("te(", cov_used_B, fixed = TRUE),
                                  grep(":", cov_used_B, fixed = TRUE))])
    inter_used_A <- cov_used_A[c(grep("te(", cov_used_A, fixed = TRUE),
                                 grep(":", cov_used_B, fixed = TRUE))]  # interaction of continuous factors
    inter_used_B <- matrix(ncol = 2, nrow = length(inter_used_A))
    if (length(grep("GLM", modeltype)) == 1) {
      # For GLM
      inter_used_B <- t(apply(t(inter_used_A), 2, function(x)
        strsplit(x, split = ":")[[1]][1:2]))
      wpol <- unique(grep("poly(", inter_used_B, fixed = TRUE),
                     grep("ns(", inter_used_B, fixed = TRUE))
      if (length(wpol) > 0) {
        inter_used_B[wpol] <- c(apply(t(inter_used_B[wpol]), 2,
                                      function(x) strsplit(x, split = ",", fixed = TRUE)[[1]][1]))
        inter_used_B[wpol] <- apply(t(inter_used_B[wpol]), 2, function(x)
          strsplit(x, split = "(", fixed = TRUE)[[1]][2])
      }
      inter_used <- inter_used_B
    } else {
      # For GAM
      interWithfact <- grep("by=", inter_used_A, fixed = TRUE)
      if (length(interWithfact) == 0) {
        inter_used_B <- t(apply(t(inter_used_A), 2, function(x)
          strsplit(x, split = ",")[[1]][1:2]))
        inter_used_B[, 1] <- apply(t(inter_used_B[, 1]), 2, function(x)
          strsplit(x, split = "(", fixed = TRUE)[[1]][2])
      }
      if (length(interWithfact) != 0 &
          length(interWithfact) < length(inter_used_A)) {
        inter_used_B[-interWithfact, ] <- t(apply(
          t(inter_used_A[-interWithfact]),
          2, function(x) strsplit(x, split = ",")[[1]][1:2]))
        inter_used_B[-interWithfact, 1] <- apply(
          t(inter_used_B[-interWithfact, 1]), 2, function(x)
            strsplit(x, split = "(", fixed = TRUE)[[1]][2])
        inter_used_B[interWithfact, ] <- t(apply(
          t(inter_used_A[interWithfact]), 2, function(x)
            strsplit(x, split = ",")[[1]][c(1, 3)]))
        inter_used_B[interWithfact, 2] <- apply(
          t(inter_used_B[interWithfact, 2]), 2, function(x)
            strsplit(strsplit(x, split = "=", fixed = TRUE)[[1]][2],
                     split = ")")[[1]][1])
        inter_used_B[interWithfact, 1] <- apply(
          t(inter_used_B[interWithfact, 1]), 2, function(x)
            strsplit(x, split = "(", fixed = TRUE)[[1]][2])
      }
      if (length(interWithfact) == length(inter_used_A)) {
        inter_used_B <- t(apply(t(inter_used_A), 2, function(x)
          strsplit(x, split = ",")[[1]][c(1, 3)]))
        inter_used_B[, 1] <- apply(t(inter_used_B[, 1]), 2, function(x)
          strsplit(x, split = "(", fixed = TRUE)[[1]][2])
        inter_used_B[, 2] <- apply(t(inter_used_B[, 2]), 2, function(x)
          strsplit(strsplit(x,  split = "=", fixed = TRUE)[[1]][2],
                   split = ")")[[1]][1])
      }
      inter_used_C <- inter_used_B
      if (length(inter_used_C) >= 3) {
        if (length(grep("k=", inter_used_C[, 1])) >= 1 &
            length(grep("k=", inter_used_C[, 1])) != nrow(inter_used_C)) {
          inter_used <- inter_used_C[-grep("k=", inter_used_C[, 1]), ]
        }
      }
      if (length(grep("k=", inter_used_C[, 1])) == 0) {
        inter_used <- inter_used_C
      }
    }
  } else {
    cov_used_C <- t(cov_used_B)
  }
  cov_used_D <- apply(cov_used_C, 2, function(x) {
    res <- x
    if (grepl("GLM", modeltype)) {
      if (grepl("ns", modeltype)) {
        if (length(strsplit(x, split = "ns(", fixed = TRUE)[[1]]) == 2) {
          res1 <- strsplit(x, split = "ns(", fixed = TRUE)[[1]][2]
          res <- strsplit(res1, split = ",", fixed = TRUE)[[1]][1]
        }
      } else {
        if (length(strsplit(x, split = "poly(", fixed = TRUE)[[1]]) == 2) {
          res1 <- strsplit(x, split = "poly(", fixed = TRUE)[[1]][2]
          res <- strsplit(res1, split = ",", fixed = TRUE)[[1]][1]
        }
      }
    } else {
      if (length(strsplit(x, split = "s(", fixed = TRUE)[[1]]) == 2) {
        res <- strsplit(x, split = "s(", fixed = TRUE)[[1]][2]
      }
    }
    if (length(strsplit(x, split = "as.character(", fixed = TRUE)[[1]]) == 2) {
      res1 <- strsplit(x, split = "as.character(", fixed = TRUE)[[1]][2]
      res <- strsplit(res1, split = ")", fixed = TRUE)[[1]][1]
    }
    res
  })
  cov_used <- cov_used_D #[order(cov_used_D)]
  # NumFact <- 0

  # discretisation of numeric factors for graph representation
  ## S'assurer que le nombre de levels est suffisant pour qu'il y ait un tableau à travailler
  #   IsDataNumeric <- apply(t(cov_used), 2, function(x) is.numeric(Y_data_sample[, x]))
  IsDataNumeric <- (Y_data_sample %>%
                      dplyr::select(dplyr::one_of(cov_used)) %>%
                      sapply(., class) == "numeric")

  save(cov_used, file = paste0(saveWD, "/cov_used_", Num, ".RData"))
  save(IsDataNumeric, file = paste0(saveWD, "/IsDataNumeric_", Num, ".RData"))

  if (sum(grepl("factor", cov_used) == sum(!IsDataNumeric)) !=
      length(IsDataNumeric)) {
    warnings("Verify if all factors have been accounted",
             "as factors in the model selection procedure !")
    IsDataNumeric[grepl("factor", cov_used)] <- FALSE
  }

  n.MargFig <- length(cov_used)
  if (exists("inter_used")) {n.MargFig <- length(cov_used) + length(inter_used)}
  save(n.MargFig, file = paste0(saveWD, "/n.MargFig_", Num, ".RData"))

  # Marginal predictions for covariates ====
  if (Marginals) {
    nb_lev <- 15
    if (sum(IsDataNumeric) >= 2 & length(IsDataNumeric) >= 6) {
      nb_lev <- 10
    }

    AllFact <- lapply(cov_used, function(f) { # f <- cov_used[4]
      #     y <- Y_data_sample[, which(names(Y_data_sample) == f)]
      y <- Y_data_sample %>% dplyr::select(dplyr::matches(f))

      if (IsDataNumeric[cov_used == f]) {
        seq(min(y), max(y), (max(y) - min(y))/nb_lev)
      } else {
        unique(unlist(y))
      }
    })
    names(AllFact) <- cov_used

    # Create the fractionnar plan for predictions
    frac.key <- planor::planor.designkey(
      factors = names(AllFact),
      nlevels = lengths(AllFact),
      model = as.formula(paste0("~(", paste(names(AllFact), collapse = "+"), ")^2")),
      base = as.formula(paste0("~(", paste(names(AllFact), collapse = "+"), ")"))
    )

    frac.plan <- planor::planor.design(frac.key)

    GridFact <- as.data.frame(apply(t(1:length(AllFact)), 2, function(x) {
      AllFact[[x]][frac.plan@design[,x]]
    }))

    if (sum(IsDataNumeric) != length(!IsDataNumeric)) {
      GridFact <- GridFact %>%
        dplyr::mutate_if(.predicate = !IsDataNumeric, .funs = function(x) as.factor(as.character(x))) %>%
        dplyr::mutate_if(.predicate = IsDataNumeric, .funs = function(x) as.numeric(as.character(x)))
    }
    names(GridFact) <- names(AllFact)
    # head(GridFact)
    rm(frac.key, frac.plan); gc()

    # Calculate predictions while accounting for uncertainty of prediction
    if (grepl("PA", modeltype)) {
      PredGridFact <- cbind(GridFact, predict(modelX, newdata = GridFact,
                                              type = "link", se.fit = TRUE))
      names(PredGridFact) <- c(names(GridFact), "pred", "se")
      PredGridFact_boot2 <- t(apply(PredGridFact[, c("pred", "se")], 1,
                                    function(x) {
                                      1/(1 + exp(-qnorm(c(0.5, 0.05, 0.95), mean = x[1], sd = x[2])))
                                    }))
    } else if (grepl("Cont", modeltype) & !grepl("Log", modeltype)) {
      PredGridFact <- cbind(GridFact, predict(modelX, newdata = GridFact,
                                              type = "response", se.fit = TRUE))
      names(PredGridFact) <- c(names(GridFact), "pred", "se")
      PredGridFact_boot2 <- t(apply(PredGridFact[, c("pred", "se")], 1,
                                    function(x) {
                                      qnorm(c(0.5, 0.05, 0.95), mean = x[1], sd = x[2])
                                    }))
    } else if (grepl("Gamma", modeltype)) {
      PredGridFact <- cbind(GridFact, predict(modelX, newdata = GridFact,
                                              type = "link", se.fit = TRUE))
      names(PredGridFact) <- c(names(GridFact), "pred", "se")
      # Gamma distribution
      # beta = scale = (se*se)/pred
      # alpha= shape = 1/(CV*CV)
      PredGridFact_boot2 <- t(apply(PredGridFact[, c("pred", "se")], 1,
                                    function(x) {
                                      exp(qnorm(c(0.5, 0.05, 0.95), mean = x[1], sd = x[2]))
                                    }))
    } else if (grepl("Log", modeltype)) {
      PredGridFact <- cbind(GridFact, predict(modelX, newdata = GridFact,
                                              type = "response", se.fit = TRUE))
      names(PredGridFact) <- c(names(GridFact), "pred", "se")
      # Log-normal distribution with Laurent transformation
      PredGridFact_boot2 <- t(apply(PredGridFact[, c("pred", "se")], 1,
                                    function(x) {
                                      exp(qnorm(c(0.5, 0.05, 0.95), mean = x[1],
                                                sd = x[2]) + 0.5 * var(modelX$residuals))
                                    }))
    } else if (grepl("Count", modeltype)) {
      PredGridFact <- cbind(GridFact,
                            predict(modelX, newdata = GridFact, type = "link", se.fit = TRUE))
      names(PredGridFact) <- c(names(GridFact), "pred", "se")
      PredGridFact_boot2 <- t(apply(PredGridFact[, c("pred", "se")], 1,
                                    function(x) {
                                      exp(qnorm(c(0.5, 0.05, 0.95), mean = x[1],
                                                sd = x[2]) + 0.5 * var(modelX$residuals))
                                    }))
    } else if (grepl("Tweed", modeltype)) {
      PredGridFact <- cbind(GridFact,
                            predict(modelX, newdata = GridFact, type = "response", se.fit = TRUE))

      maxY <- max(Y_data_sample$dataY) * 1.1
      # Ne pas recalculer pour toutes les valeurs de fit
      # Faire pour 10^-2 la précision des datas.
      PredGridFact$fit2 <- round(PredGridFact$fit * 10^-prec.data) * 10^prec.data
      PredGridFact$fit2[which(PredGridFact$fit2 == min(PredGridFact$fit2, na.rm = TRUE))] <- min(PredGridFact$fit, na.rm = TRUE)
      fit.unique <- data.frame(unique(PredGridFact$fit2))
      names(fit.unique) <- "fit2"

      if (cl_inside) {
        cl <- parallel::makePSOCKcluster(nbclust)
      }
      fit.unique2 <- t(snow::parApply(
        cl, t(fit.unique), 2,
        function(x, powerXI, modelX, maxY) {
          if (x > maxY) {
            ## Calculation is very long for too high values
            a <- rep(maxY, 3)
          } else {

            a <- rep(NA, 3)
            if (x == 0) {
              a <- rep(0, 3)
            } else {
              try(a <- tweedie::qtweedie(p = c(0.5, 0.05, 0.95), xi = powerXI,
                                         mu = x, phi = summary(modelX)$dispersion))
            }
          }
          a
        }, powerXI = Param, modelX = modelX, maxY = maxY))
      if (cl_inside) {
        parallel::stopCluster(cl)
      }

      PredGridFact_boot2 <- dplyr::left_join(PredGridFact, cbind(fit.unique, fit.unique2)) %>%
        dplyr::select(dplyr::one_of(c("1","2","3")))
    }

    # Add +1 when Sup = 1
    if (model_selected$Sup != 0) {
      PredGridFact_boot2 <- PredGridFact_boot2 - model_selected$Sup
    }

    if (grepl("PA", modeltype)) {
      ymax <- 1
      ymin <- 0
    } else {
      ymax <- quantile(PredGridFact_boot2[,1], probs = 0.975, na.rm = TRUE)
      ymin <- quantile(PredGridFact_boot2[,1], probs = 0.025, na.rm = TRUE)
    }
    MaxMed <- ymax

    # Plot Marginals ----
    if (plot) {
      # One figure for each covariate or interaction
      # n.fig <- 0
      if (cl_inside) {
        cl <- parallel::makePSOCKcluster(nbclust)
      }
      if (!exists("Seuil")) {Seuil <- NA}

      tmp <- snow::parLapply(
        cl, 1:length(cov_used),
        function(i,
                 cov_used, saveWD, imgprefix, Num,
                 PredGridFact, PredGridFact_boot2,
                 modeltype, Y_data_sample, IsDataNumeric,
                 ymin, ymax, Seuil)
        {
          f <- cov_used[i]

          grid.pred.2 <- data.frame(PredGridFact[,f], PredGridFact_boot2)
          names(grid.pred.2) <- c("grid", "med", "min", "max")
          if (!IsDataNumeric[cov_used == f]) {
            grid.pred.2[,1] <- as.numeric(grid.pred.2[,1])
          }
          # n.fig <- n.fig + 1
          png(filename =
                paste0(saveWD, "/", imgprefix, "Marginal_Effects_", Num, "_", i, ".png"),
              # width = fig.dim$ncol * 5, height = fig.dim$nrow * 5,
              width = 8, height = 8,
              units = "cm", pointsize = 8, res = 300)
          # x11()
          # par(mfrow = c(fig.dim$nrow, fig.dim$ncol), mar = c(3.5, 3, 0.5, 0.5))
          par(mar = c(3.25, 3, 0.5, 0.5))
          # Figures of marginals of covariates
          # for (f in cov_used) { # f <- cov_used[1]
          wdth <- sort(unique(grid.pred.2[,1]))[2] - sort(unique(grid.pred.2[,1]))[1]
          if (IsDataNumeric[cov_used == f]) {
            plot(NA, xlim = c(min(grid.pred.2[,1]) - 0.5 * wdth,
                              max(grid.pred.2[,1]) + 0.5 * wdth),
                 ylim = c(ymin, ymax), #main = f, cex.main = 1,
                 ylab = NA, xlab = NA, xaxs = "i")
          } else {
            plot(NA,
                 xlim = c(min(grid.pred.2[,1] - 0.5 * wdth),
                          max(grid.pred.2[,1]) + 0.5 * wdth),
                 ylim = c(ymin, ymax), #main = f, cex.main = 1,
                 xaxt = "n", ylab = NA, xlab = NA, xaxs = "i")
            axis(1, at = 1:length(unique(PredGridFact[,f])))
          }
          title(line = 2, xlab = f, ylab = "Y", cex.lab = 1.2)
          # tmp <- lapply(1:length(unique(PredGridFact[,f])), function(nb.lev)
          #   .marge_cov(nb.lev, grid = PredGridFact[,f], grid.pred = PredGridFact_boot2))

          for (z in 1:3) {
            p <- c("med", "min", "max")[z]
            col <- c("grey30", "forestgreen", "orangered")[z]
            m.wdth <- c(0.49, 0.3, 0.2)[z]
            pt.o <- c(0.005, 0, 0)[z]

            for (val in sort(unique(grid.pred.2$grid))) {
              # val <- sort(unique(grid.pred.2$grid))[6]
              pred <- as.data.frame(grid.pred.2[which(grid.pred.2$grid == val), p])
              if (sum(!duplicated(pred)) == 1) {
                boxplot(
                  pred, at = val,
                  boxwex = 2 * wdth * m.wdth, border = col,
                  add = TRUE)
              } else {

                yarrr::pirateplot(formula = as.formula(paste(p, "~ grid")),
                                  data = as.data.frame(
                                    grid.pred.2[which(grid.pred.2$grid == val), ]),
                                  # avg.line.fun = median, inf = "iqr",
                                  avg.line.fun = mean,
                                  point.o = pt.o,# 0.05 * m.wdth,
                                  point.cex = 0.5,
                                  bean.b.o = 1, bean.lwd = 0.8,
                                  avg.line.lwd = 0.8,
                                  inf.f.o = 0.35, inf.b.o = 0.35,
                                  pal = col,
                                  # bw = "bcv",
                                  jitter.val = 0.1 * wdth * m.wdth,
                                  at = val,
                                  width.min = wdth * m.wdth * 0.99,
                                  width.max = wdth * m.wdth,
                                  xaxt = "n", yaxt = "n",
                                  add = TRUE)
              }
            }

          }
          # rm(tmp)
          if (modeltype %in% c("PA", "PAGLM")) {
            abline(h = Seuil, lwd = 2, lty = "dashed")
          }
          # rug(Y_data_sample[,which(names(Y_data_sample) == f)], ticksize = 0.02)
          rug(unlist(dplyr::select(Y_data_sample, dplyr::matches(f))), ticksize = 0.02)
          dev.off()
        },
        cov_used = cov_used, saveWD = saveWD, imgprefix = imgprefix, Num = Num,
        PredGridFact = PredGridFact, PredGridFact_boot2 = PredGridFact_boot2,
        modeltype = modeltype, Y_data_sample = Y_data_sample,
        IsDataNumeric = IsDataNumeric, ymin = ymin, ymax = ymax, Seuil = Seuil
      ) # End of parApply on cov_used

      # Figures of marginales of interactions
      if (exists("inter_used")) {
        if (!is.matrix(inter_used)) {inter_used <- matrix(inter_used, ncol = 2)}
        # if (length(inter_used) >= 3) {
        tmp <- snow::parLapply(cl, (length(cov_used) + 1):n.MargFig,
                         function(i,
                                  inter_used, cov_used, saveWD, imgprefix, Num,
                                  PredGridFact, PredGridFact_boot2,
                                  modeltype, Y_data_sample, IsDataNumeric,
                                  ymin, ymax, MaxMed)
                         {
                           #for (inter in 1:nrow(inter_used)) {
                           inter <- i - length(cov_used)
                           # n.fig <- n.fig + 1
                           png(filename =
                                 paste0(saveWD, "/", imgprefix, "Marginal_Effects_", Num, "_", i, ".png"),
                               # width = fig.dim$ncol * 5, height = fig.dim$nrow * 5,
                               width = 5, height = 5,
                               units = "cm", pointsize = 8, res = 300)
                           # x11()
                           # par(mfrow = c(fig.dim$nrow, fig.dim$ncol), mar = c(3.5, 3, 0.5, 0.5))
                           par(mar = c(3.5, 3, 0.5, 0.5))

                           MEAN <- tapply(PredGridFact_boot2[, 1],
                                          list(PredGridFact[,
                                                            which(names(PredGridFact) == inter_used[inter, 1])],
                                               PredGridFact[,
                                                            which(names(PredGridFact) == inter_used[inter, 2])]),
                                          FUN = median)
                           if (is(PredGridFact[,
                                               which(names(PredGridFact) == inter_used[inter, 2])])[1]
                               == "numeric") {
                             x1 <- as.numeric(as.character(colnames(MEAN)))
                           } else {
                             x1 <- as.numeric(as.factor(colnames(MEAN)))
                           }
                           if (is(PredGridFact[,
                                               which(names(PredGridFact) == inter_used[inter, 1])])[1]
                               == "numeric") {
                             y1 <- as.numeric(as.character(rownames(MEAN)))
                           } else {
                             y1 <- as.numeric(as.factor(rownames(MEAN)))
                           }
                           image(x1, y1, t(MEAN), col = rev(heat.colors(255)),
                                 breaks = c(seq(ymin, ymax, (ymax - ymin)/254),
                                            MaxMed), xlab = NA, ylab = NA)
                           title(line = 2, xlab = inter_used[inter,2],
                                 ylab = inter_used[inter,1])
                           #         rug(Y_data_sample[,
                           #                           which(names(Y_data_sample) == inter_used[inter, 2])])
                           #
                           rug(unlist(dplyr::select(Y_data_sample, dplyr::matches(inter_used[inter, 2]))),
                               ticksize = 0.02)
                           rug(unlist(dplyr::select(Y_data_sample, dplyr::matches(inter_used[inter, 1]))),
                               ticksize = 0.02, side = 2)

                           #         rug(Y_data_sample[,
                           #                           which(names(Y_data_sample) == inter_used[inter, 1])], side = 2)
                           points(
                             unlist(dplyr::select(Y_data_sample, dplyr::matches(inter_used[inter, 2]))),
                             unlist(dplyr::select(Y_data_sample, dplyr::matches(inter_used[inter, 1]))),
                             #           Y_data_sample[, which(names(Y_data_sample) == inter_used[inter,2])],
                             #           Y_data_sample[, which(names(Y_data_sample) == inter_used[inter,1])],
                             pch = 20, cex = 0.1)
                           if (modeltype %in% c("PA", "PAGLM")) {
                             contour(x1, y1, t(MEAN), levels = Seuil, lwd = 2,
                                     add = TRUE)
                           }
                           dev.off()
                         },# End of inter_used
                         cov_used = cov_used, saveWD = saveWD, imgprefix = imgprefix, Num = Num,
                         PredGridFact = PredGridFact, PredGridFact_boot2 = PredGridFact_boot2,
                         modeltype = modeltype, Y_data_sample = Y_data_sample,
                         IsDataNumeric = IsDataNumeric, ymin = ymin, ymax = ymax, Seuil = Seuil,
                         MaxMed = MaxMed
        ) # end of parLapply
      } # End of if(exists("inter_used"))
      if (cl_inside) {
        parallel::stopCluster(cl); gc()
      }
    }  # end of plot
  }

  # Simple Marginals close to mean or 0 (PA) ----
  AllFact.simple <- lapply(cov_used, function(f) { # f <- cov_used[4]
    #     y <- Y_data_sample[, which(names(Y_data_sample) == f)]
    y <- Y_data_sample %>% dplyr::select(dplyr::matches(f))

    if (IsDataNumeric[cov_used == f]) {
      # seq(min(y), max(y), (max(y) - min(y))/)
      res <- data.frame(
        unique(quantile(unlist(y), probs = seq(0, 1, 0.05))),
        CovForMeanPred[,cov_used[-which(cov_used == f)]],
        cov = f,
        cov.val = unique(quantile(unlist(y), probs = seq(0, 1, 0.05))))
      names(res)[1] <- f
    } else {
      res <- data.frame(
        unique(unlist(y)),
        CovForMeanPred[,cov_used[-which(cov_used == f)]],
        cov = f,
        cov.val = as.numeric(as.character(unique(unlist(y)))))
      names(res)[1] <- f
    }
    res %>%
      dplyr::mutate_if(is.factor, as.character)
  })
  AllFact.simple.tbl <- dplyr::as.tbl(do.call(dplyr::bind_rows, AllFact.simple)) %>%
    dplyr::bind_cols(data.frame(predict(modelX, newdata = .,
                                 type = "response", se.fit = TRUE))) %>%
    dplyr::mutate(cov2 = factor(cov, levels = cov_used))

  if (plot) {
    figdim <- Fig_split(length(cov_used))

    g <- ggplot(AllFact.simple.tbl) +
      geom_pointrange(aes(x = cov.val, y = fit,
                          ymin = fit-se.fit,
                          ymax = fit+se.fit)) +
      facet_wrap(~cov2, scales = "free_x",
                 nrow = figdim$nrow,
                 ncol = figdim$ncol)
    rm(AllFact.simple)

    ggsave(plot = g,
           filename = paste0(
             saveWD, "/", basename(saveWD),
             "-CovForMeanPred_Marginals_", model_selected$Num, ".png"),
           width = 16, height = 16*(figdim$nrow / figdim$ncol), units = "cm",
           dpi = 300)
  }

  # Comparison of predictions against observations ====
  if (!grepl("PA", modeltype) & plot) {

    if (grepl("Log", model_selected$modeltype)) {
      fits <- exp(modelX$fitted + 0.5 * var(residuals(modelX)))
    } else {
      fits <- modelX$fitted
    }

    png(filename =
          paste0(saveWD, "/", imgprefix, "Obs-Pred_", Num, ".png"),
        width = 12, height = 6,
        units = "cm", pointsize = 8, res = 300)
    par(mai = c(0.35, 0.35, 0.1, 0.1), mfrow = c(1,2), cex.axis = 0.8)

    smoothScatter(Y_data_sample$dataY, fits, nbin = 256,
                  xlab = "", ylab = "", col = "grey")
    mtext(side = 1, line = 1.75, "Observations")
    mtext(side = 2, line = 1.75, "Predictions in obs. scale")
    abline(b = 1, a = 0, lty = "dashed", col = "grey")

    smoothScatter(Y_data_sample$dataY, Y_data_sample$dataY - fits, nbin = 256,
                  colramp = colorRampPalette(c("white", "red")),
                  xlab = "", ylab = "", col = "grey")
    mtext(side = 1, line = 1.75, "Observations")
    mtext(side = 2, line = 1.75, "Obs-Pred Difference")
    abline(h = 0, lty = "dashed", col = "grey")
    dev.off()

  }
  if (grepl("PA|Tweed", modeltype) & plot)
  {
    png(filename =
          paste0(saveWD, "/", imgprefix, "Thresholds_for_prediction_", Num, ".png"),
        width = 9, height = 8,
        units = "cm", pointsize = 7, res = 300)
    par(mai = c(0.3, 0.5, 0.3, 0.9), yaxs = "i")
    if (grepl("Tweed", modeltype)) {
      dataBin <- unlist(1 * (Y_data_sample$dataY > 0))
      dataFit <- modelX$proba
    } else {
      dataBin <- unlist(Y_data_sample$dataY)
      dataFit <- modelX$fitted
    }
    pirateplot(formula = dataFit ~ dataBin,
               data = Y_data_sample,
               avg.line.fun = mean,
               point.o = .05, point.cex = 0.5,
               bean.b.o = 1, bean.lwd = 1,
               avg.line.lwd = 0.8,
               inf.f.o = 0.35, inf.b.o = 0.35,
               pal = "black",
               bw = "bcv",
               xlab = "",
               ylab = "Model prediction (probability of presence)",
               main = "Comparison of observation and prediction",
               ylim = c(0, 1)
    )
    abline(h = c(SeuilMinUn, SeuilMaxZero), col = c("red", "olivedrab3"),
           lty = "dashed")
    segments(1.5, SeuilMinUn, 3, SeuilMinUn, col = "red", lwd = 2)
    segments(1.5, SeuilMaxZero, 0, SeuilMaxZero, col = "olivedrab3",
             lwd = 2)
    mtext(side = 1, line = 1.5, text = "Field observation")
    axis(4, at = SeuilMinUn, labels = round(SeuilMinUn, digits = 2),
         col.axis = "red", las = 1)
    axis(4, at = SeuilMaxZero, labels = round(SeuilMaxZero, digits = 2),
         col.axis = "olivedrab3", las = 1)
    axis(4, at = Seuil, labels = round(Seuil, digits = 2),
         col.axis = "forestgreen", las = 1)
    abline(h = Seuil, col = "forestgreen", lwd = 2.5)
    # axis(4,at=Seuil,labels=round(Seuil,digits=2),col.axis='forestgreen',las=1)
    text(1.5, Seuil, "Best THD", pos = 3, col = "forestgreen", font = 2)

    # segments(2.9,SeuilMaxZero,3,SeuilMaxZero,col='olivedrab3',lty='dashed',xpd=TRUE)
    segments(2.6, 0, 3, 0, lty = "dashed", xpd = TRUE)
    segments(2.6, 1, 3, 1, lty = "dashed", xpd = TRUE)
    arrows(3, 0, 3, SeuilMaxZero, xpd = TRUE, code = 3, length = 0.05)
    text(3.1, 0, "95% of observation of absence", xpd = TRUE, srt = 270, adj = 1)
    arrows(3, 1, 3, SeuilMaxZero, xpd = TRUE, code = 3, length = 0.05,
           col = "olivedrab3")
    text(3.1, 1, "Sure presence", xpd = TRUE, srt = 270, col = "olivedrab3",
         adj = 0)

    segments(3.2, SeuilMinUn, 3.3, SeuilMinUn, col = "red", lty = "dashed",
             xpd = TRUE)
    segments(3.2, 0, 3.3, 0, lty = "dashed", xpd = TRUE)
    segments(3.2, 1, 3.3, 1, lty = "dashed", xpd = TRUE)
    arrows(3.3, 0, 3.3, SeuilMinUn, xpd = TRUE, code = 3, length = 0.05,
           col = "red")
    text(3.4, 0, "Sure absence", xpd = TRUE, srt = 270, col = "red",
         adj = 1)
    arrows(3.3, 1, 3.3, SeuilMinUn, xpd = TRUE, code = 3, length = 0.05)
    text(3.4, 1, "95% of observation of presence", xpd = TRUE, srt = 270,
         adj = 0)
    box()
    text(0.65, SeuilMaxZero, "MaxZero", pos = 3, col = "olivedrab3", xpd = TRUE)
    text(2.35, SeuilMinUn, "MinUn", pos = 1, col = "red", xpd = TRUE)
    dev.off()

  }

  # Zip all outputs in a unique file.
  if (zip.file != FALSE) {
    # Options j allows to save directly files without path
    utils::zip(zip.file, files = saveWD, flags = "-urj9X", extras = "",
        zip = Sys.getenv("R_ZIPCMD", "zip"))
    # Return the file path of the zip file
    return(zip.file)
  } else {
    # Return the folder path where all results have been saved
    return(saveWD)
  }
}  # end of ModelResults
