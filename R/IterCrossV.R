#' @title Find best model
#' @description Multiple k-fold cross-validation to test for the best combination
#' of covariates, model type and family.
#' @param x dataset (dataframe or SpatialPointsDataFrame) prepared with
#' \code{\link{Prepare_dataset}} function.
#' All column are supposed to be covariates except 'dataY' column.
#' Factor/Character columns should have a name starting with 'factor_'
#' @param datatype The data type to be chosen for modelling among
#' 'PA', 'Density', 'ContPosNull', 'Count', 'TweedGLM', 'KrigeGLM'.
#' (See \code{\link{modelselect_opt}} for more options and information)
#' @param corSpearman dataframe of correlation between covariates as calculated by
#' function \code{\link{Param_corr}}
#' @param saveWD directory where to save all outputs of the cross-validation procedure.
#' Folder is created if not exists. Tmp file if not defined.
#' @param zip.file Logical or file path where to save all outputs in zip. In tmpdir by default.
#' If FALSE, outputs are not zipped and the path to saveWD is returned.
#' @param restart numeric vector. If you stopped the analysis for any reasons, you can
#' restart it at the modeltype step you want. Provide a vector of values,
#' so that all modeltypes with the corresponding positions will be re-calculated.
#' In this procedure, the modelselect_opt.save file saved in saveWD will be loaded.
#' You can unzip your previously saved file and define this unzip folder as saveWD.
#' @param MPIcalc Logical. Whether the function is run within a MPI cluster or locally.
#' @param verbose Numeric. 0: no message, 1:few messages, 2:all messages
#'
#' @export

findBestModel <- function(x, datatype, corSpearman,
                          saveWD = paste0(tempdir(), "/outputs"),
                          zip.file = TRUE,
                          restart = NULL, MPIcalc = FALSE,
                          verbose = 0) {

  saveWD <- normalizePath(saveWD)
  if (!dir.exists(saveWD)) {dir.create(saveWD)}

  if (!missing(datatype)) {
    modelselect_opt$datatype <- datatype
  } else {
    datatype <- modelselect_opt$datatype
  }

  nbclust <- modelselect_opt$nbclust

  if (is.null(restart)) {
    # Save config.file in the output folder
    # file.copy(config.file, paste0(saveWD, "/config.file.save.R"), overwrite = TRUE)
    readr::write_rds(modelselect_opt(), paste0(saveWD, "/modelselect_opt.save.rds"))
  } else {
    modelselect_opt(readr::read_rds(paste0(saveWD, "/modelselect_opt.save.rds")))
    message("modelselect_opt options were loaded from previous start.")
  }

  # Options ----
  Max_nb_Var <- modelselect_opt$Max_nb_Var
  modeltypes <- modelselect_opt$modeltypes
  graph <- modelselect_opt$graph
  Max_K <- modelselect_opt$Max_K
  Max_K_te <- modelselect_opt$Max_K_te
  Max_K_Poly <- modelselect_opt$Max_K_Poly
  fixXI <- modelselect_opt$fixXI

  Interaction <- modelselect_opt$Interaction
  MinNbModel <- modelselect_opt$MinNbModel
  MaxDist <- modelselect_opt$MaxDist
  Phi <- modelselect_opt$Phi
  Model <- modelselect_opt$Model
  lcc_proj <- modelselect_opt$lcc_proj
  Lambda <- modelselect_opt$Lambda
  fix.Lambda <- modelselect_opt$fix.Lambda

  k_fold <- modelselect_opt$k_fold
  nbMC <- modelselect_opt$nbMC

  seqthd <- modelselect_opt$seqthd
  lim_pvalue <- modelselect_opt$lim_pvalue

  seed <- modelselect_opt$seed

  if (grepl("KrigeGLM", datatype)) {
    warning(c("You chose 'datatype = KrigeGLM*', which is still experimental, ",
              "please verify that modelselect_opt is correctly set up: ",
              "Default values are not safe at all !"))
  }

  # Data ----
  Y_data_sample <- x

  # Added for remind ----------------------
  if (is(x)[1] == "SpatialPointsDataFrame") {
    physParam_NewSp <- x@data[, -which(names(x@data) == "dataY")]
  } else {
    physParam_NewSp <- x[, -which(names(x) == "dataY")]
  }

  # MPIcalc if available ----

  if (grepl("KrigeGLM", datatype)) {
    Y_data_sample_lcc <- sp::spTransform(Y_data_sample, CRS(lcc_proj))
  } else {
    Y_data_sample_lcc <- NA
  }

  if (!MPIcalc) {
    cl <- parallel::makePSOCKcluster(nbclust)
  } else {
    cl <- snow::makeCluster()
    options(rasterClusterObject = cl)
    options(rasterClusterCores = length(cl))
    options(rasterCluster = TRUE)
    options(rasterClusterExclude = NULL)
  }

  set.seed(seed)
  # if (!is.null(restart) & file.exists(paste0(saveWD, "/saveAleaM1.RData"))) {
  if (!is.null(restart) & file.exists(paste0(saveWD, "/saveAleaM1.rds"))) {
    # load(paste0(saveWD, "/saveAleaM1.RData"))
    saveAlea <- readr::read_rds(paste0(saveWD, "/saveAleaM1.rds"))
    # load(paste0(saveWD, "/saveAlea.nbM1.RData"))
    saveAlea.nb <- readr::read_rds(paste0(saveWD, "/saveAlea.nbM1.rds"))
    warning("If you changed part of the dataset, you should remove saveAlea*.RData files before running this script")
  } else {
    # Generate 10 samples for a 10-fold cross-validation ----
    # 10 times the 10-fold cross-validation
    for (MC in 1:(nbMC/k_fold)) {
      saveAlea_tmp <- matrix(nrow = trunc(length(Y_data_sample$dataY)/k_fold) +
                               2, ncol = k_fold)
      if (grepl("PA", datatype)) {
        # Sample randomly but with the same 0/1 ratio than the dataset A
        # VERIFIER ##
        samp0 <- sample(rep(1:k_fold, length = sum(Y_data_sample$dataY ==
                                                     0)))
        samp1 <- sample(rep(1:k_fold, length = sum(Y_data_sample$dataY ==
                                                     1)))
        for (i in 1:k_fold) {
          saveAlea_tmp[, i] <- c(which(Y_data_sample$dataY == 0)[which(samp0 ==
                                                                         i)], which(Y_data_sample$dataY == 1)[which(samp1 ==
                                                                                                                      i)], rep(NA, nrow(saveAlea_tmp) - sum(c(samp0, samp1) ==
                                                                                                                                                              i)))
        }
      } else {
        samp10 <- sample(rep(1:k_fold, length = length(Y_data_sample$dataY)))

        for (i in 1:k_fold) {
          # i <- 1
          saveAlea_tmp[, i] <- c(which(samp10 == i), rep(NA, nrow(saveAlea_tmp) -
                                                           sum(samp10 == i)))
          # Add a data to fit dataset if a level is missing
          if (length(grep("factor", names(physParam_NewSp))) !=
              0) {
            # all class factors
            class.in <- names(physParam_NewSp)[grep("factor",
                                                    names(physParam_NewSp))]
            # Test if some levels not in fit dataset. Add one randomly.
            for (class.data in class.in) {
              # lev.not.in <- levels(physParam_NewSp[, class.data])[(!levels(physParam_NewSp[,
              # class.data]) %in% physParam_NewSp[na.omit(saveAlea_tmp[,
              # i]), class.data])]
              # if (length(lev.not.in) != 0) { add.line <- apply(t(1 :
              # length(lev.not.in)), 2, function(x) sample(which(physParam_NewSp[,
              # class.data] == lev.not.in[x]), size = 1))
              # saveAlea_tmp[is.na(saveAlea_tmp[,i]),i][1 : length(add.line)] <-
              # add.line } Remove new levels in validation dataset when factor data
              # This method is prefered. At the end, RMSE t.test will be weighted
              # according to proportion of validation set on total dataset test if
              # new levels in valid
              test.out <- apply(apply(t(class.in), 2, function(x) {
                !(physParam_NewSp[na.omit(saveAlea_tmp[, i]), x] %in%
                    physParam_NewSp[-na.omit(saveAlea_tmp[, i]),
                                    x])
              }), 1, sum)
              saveAlea_tmp[test.out != 0, i] <- NA
            }
          }
        }
      }
      if (MC == 1) {
        saveAlea <- saveAlea_tmp
      } else {
        saveAlea <- cbind(saveAlea, saveAlea_tmp)
      }
    }
    # Proportion of observation in the validation dataset compared to total
    saveAlea.nb <- apply(saveAlea, 2, function(x) {
      sum(!is.na(x))
    })/nrow(Y_data_sample)
  }


  # Save critical information for output analysis ------------------------------
  Allinfo_all <- list(data = x, datatype = datatype, modeltypes = modeltypes,
                      MaxDist = MaxDist, Phi.lim = Phi, fixXI = fixXI, Interaction = Interaction)
  # save(Allinfo_all, file = paste0(saveWD, "/Allinfo_all.RData"))
  readr::write_rds(Allinfo_all, paste0(saveWD, "/Allinfo_all.rds"))

  if (!is.null(restart)) {
    all_datac <- restart
  } else {
    all_datac <- 1:length(modeltypes)
  }

  # Run loops separately for each sub-type of model ----------------------------
  for (n_datac in all_datac) {
    # n_datac <- 1
    modeltype <- modeltypes[n_datac]

    if (verbose >= 1) {
      print(" ----------------------------- ")
      print(paste("modeltype:", modeltype))
    }

    # Save information
    Allinfo <- list(modeltype = modeltypes[n_datac], Model = Model[n_datac],
                    fix.lambda = fix.Lambda[n_datac], Lambda = Lambda[n_datac])


    # For KrigeGLM compatibility
    fix.lambda <- fix.Lambda[n_datac]

    # _List all variables tested remove model with correlated ones ----

    allParamN <- names(physParam_NewSp)
    combN_combN <- as.matrix(corSpearman[, c("Var1_num", "Var2_num")])  #[,4:5]
    corSpearmanValid_combN <- corSpearman[, "valid"]  #corSpearmanValid
    if (grepl("GLM", modeltype)) {
      if (sum(grepl("factor_", allParamN)) == 0) {
        ParamNoClass <- allParamN
      } else {
        ParamNoClass <- allParamN[-grep("factor_", allParamN)]
      }
      if (length(ParamNoClass) != 0) {
        for (Pn in 2:Max_K_Poly) {
          # Add parameters for polynoms to allow different tests on polynom
          # degrees
          newParam <- paste0(ParamNoClass, "_P", Pn)
          allParamN <- c(allParamN, newParam)
        }
        for (Pn in 1:length(ParamNoClass)) {
          # Remove models where same variables exist with different degrees of
          # freedom
          combN_combN_A <- t(utils::combn(grep(ParamNoClass[Pn], allParamN),
                                          2))
          combN_combN <- rbind(combN_combN, combN_combN_A)
          corSpearmanValid_combN <- c(corSpearmanValid_combN, rep(1,
                                                                  nrow(combN_combN_A)))
        }
        for (Pn in 1:length(ParamNoClass)) {
          # Remove models where variables with different degrees of freedom have
          # original variable correlated to another variable Ex: If x1 correlated
          # to x2, then poly(x1) correlated to x2
          CombExist <- which(combN_combN[, 1] == grep(ParamNoClass[Pn],
                                                      allParamN)[1] | combN_combN[, 2] == grep(ParamNoClass[Pn],
                                                                                               allParamN)[1])
          CorExist <- CombExist[which(corSpearmanValid_combN[CombExist] ==
                                        1)]
          if (length(CorExist) != 0) {
            combN_combN_A <- expand.grid(combN_combN[CorExist,
                                                     ][which(combN_combN[CorExist, ] != grep(ParamNoClass[Pn],
                                                                                             allParamN)[1])], grep(ParamNoClass[Pn], allParamN)[-1])
            combN_combN_A2 <- matrix(unlist(combN_combN_A[which(combN_combN_A[,
                                                                              1] != combN_combN_A[, 2]), ]), ncol = 2)
            if (length(combN_combN_A2) != 0) {
              combN_combN <- rbind(combN_combN, combN_combN_A2)
              corSpearmanValid_combN <- c(corSpearmanValid_combN,
                                          rep(1, nrow(combN_combN_A2)))
            }
          }
        }
        dupComb <- duplicated(apply(combN_combN, 1, function(x) {
          paste(min(x), max(x), sep = "-")
        }))
        if (sum(dupComb == TRUE) != 0) {
          combN_combN <- combN_combN[-which(dupComb == TRUE), ]
          corSpearmanValid_combN <- corSpearmanValid_combN[-which(dupComb ==
                                                                    TRUE)]
        }
      }
    }


    # _Prepare for output record ----
    test_name_list <- list()  # Store the combinations of parameters to test in the models (free of correlated variables)
    test_factors <- list()  # Store test_name_list with names of parameters
    test_factor_int <- list()
    test_factors_all <- list()
    test_models <- list()  # Store models tested
    test_models_tmp <- list()  # Store temporary models tested
    Output_models <- list()  # Store AIC, GCV and Deviance explained
    BestTHD_crossV <- list()  # Store Best threshold for PA models
    DiffSelSpeTHD_crossV <- list()  # Store nb of errors for each possible threshold in PA models
    n_keep <- list()  # Store the models kept in the iterative procedure
    resAIC <- list()  # Store AIC
    resAIC_save <- list()  # Store all AIC
    resParam <- list()
    resParam_save <- list()

    # _Run the iterative loop for model selection ----
    # Must be sequential from model with 1 parameter to model with Max_nb_Var parameters
    if (Max_nb_Var > length(allParamN)) {Max_nb_Var <- length(allParamN)}

    for (nb in 1:Max_nb_Var) {
      if (verbose >= 1) {
        print(paste("nb:", nb))
      }
      # __List all models to test ----------------------------
      if (nb == 1) {
        n_keep[[1]] <- t(t(1:length(allParamN)))
        test_name_list[[1]] <- as.data.frame(t(t(n_keep[[1]])))
        names(test_name_list[[1]]) <- "Var1"
      } else {
        # All possible combinations of 'nb' parameters but not those previously
        # choosen 'n_keep'.

        # Add new parameter if there are not better ordered in
        # n_keep with same 1:(nb-2) parameter. Example : n_keep =
        # c(5-32,5-34,6-32,6-36,6-34); new could be : 5-32-34 but not 5-34-32;
        # 6-32-34 but not 6-34-32 nor 6-32-5
        test_name_tmp4 <- logical(0)
        for (ii in 1:nrow(n_keep[[nb - 1]])) {
          if (nb == 2) {
            n_keep_tmp <- as.numeric(as.character(unique(
              c(n_keep[[nb - 1]][c(1:ii), ]))))
          } else {
            # test on 1:(nb-2) == ii
            n_keep_tmp <- as.numeric(as.character(unique(
              c(n_keep[[nb - 1]][c(1:ii)[which(
                apply(t(1:ii), 2, function(x) {
                  paste(n_keep[[nb - 1]][x, 1:(nb - 2)], collapse = "-")}
                ) == paste(n_keep[[nb - 1]][ii, 1:(nb - 2)], collapse = "-"))], ],
                n_keep[[nb - 1]][c(1:ii), 1:(nb - 2)]))))
            # test on 1:(nb-2) != ii but before
          }
          if (length(n_keep_tmp) < length(allParamN)) {
            test_name_tmp4 <- rbind(test_name_tmp4,
                                    cbind(matrix(replicate(length(allParamN[-n_keep_tmp]),
                                                           cbind(n_keep[[nb - 1]][ii, ]), simplify = "matrix"),
                                                 ncol = (nb - 1), byrow = TRUE),
                                          c(1:length(allParamN))[-n_keep_tmp]))
          }
        }

        # Verification
        if (nb >= 3) {
          if (max(table(apply(test_name_tmp4[, 1:(nb - 1)], 1,
                              function(x) {
                                paste(x, collapse = "-")}))) >= length(allParamN))
          {
            print("Error in Combination of parameters")
          }
        }

        if (max(apply(test_name_tmp4, 1, function(x) {
          sum(duplicated(x))})) > 0)
        {
          print("Error : Parameters are asked twice in same model")
          test_name_tmp4 <- test_name_tmp4[which(apply(test_name_tmp4, 1,
                                                       function(x) {
                                                         sum(duplicated(x))
                                                       }) == 0), ]
        }


        # __Remove those for which there are correlated variables ----
        if (length(which(corSpearmanValid_combN == 1)) == 1) {
          tmp <- apply(test_name_tmp4, 1, function(x) {
            sum(
              apply(t(
                cbind(
                  combN_combN[, 1],
                  combN_combN[,2])[
                    which(corSpearmanValid_combN == 1), ]),
                2,
                function(y) {
                  sum(as.numeric(as.character(x)) %in%
                        as.numeric(as.character(y))) == 2
                }))
          })
        } else {
          tmp <- snow::parApply(
            cl, test_name_tmp4, 1,
            function(x, combN_combN, corSpearmanValid_combN) {
              sum(
                apply(
                  cbind(
                    combN_combN[, 1],
                    combN_combN[, 2])[
                      which(corSpearmanValid_combN == 1), ],
                  1,
                  function(y) {
                    sum(as.numeric(as.character(x)) %in%
                          as.numeric(as.character(y))) == 2
                  }))
            }, combN_combN = combN_combN,
            corSpearmanValid_combN = corSpearmanValid_combN)
        }
        # If there is no model left, break the loop
        if (sum(tmp == 0) == 0) {
          print("no more models to test")
          break
        } else {
          test_name_list[[nb]] <- test_name_tmp4[which(tmp == 0), ]
        }
      }  # end of nb>= 2

      # if there is only one model left
      if (is.null(dim(test_name_list[[nb]]))) {test_name_list[[nb]] <- t(test_name_list[[nb]])}

      # __Save list of models as names (not numbers) ----
      test_factors[[nb]] <- t(apply(test_name_list[[nb]], 1, function(x) allParamN[as.numeric(as.character(x))]))

      if (nb == 1) {
        test_factors[[nb]] <- t(test_factors[[nb]])
        if (verbose >= 2) {
          print("All covariates tested")
          print(test_factors[[nb]])
        }
      }
      if (length(nrow(test_factors[[nb]])) == 0) {
        test_factors[[nb]] <- t(test_factors[[nb]])
      }

      # __Add interactions ----
      if (nb >= 2 & Interaction)
      {
        test_factor_int[[nb]] <- t(apply(
          test_factors[[nb]],
          1, function(x) {
            # x <- test_factors[[nb]][1,]
            inter <- t(utils::combn(x, 2))
            c(apply(inter, 1, function(y) {
              # y <- inter[2,]
              res <- paste("te(", paste(y, collapse = ","), ",k=",
                           Max_K_te, ")", sep = "")
              if (length(grep("factor_", y)) == 1) {
                # if one is a factor
                res <- paste("te(", y[-grep("factor_", y)],
                             ",k=", Max_K_te, ",by=", y[grep("factor_",
                                                             y)], ")", sep = "")
              }
              if (length(grep("factor_", y)) == 2) {
                # if both are factors or in GLM
                res <- paste(y, collapse = ":")
              }
              if (grepl("GLM", modeltype)) {
                y2 <- apply(t(y), 2, function(yy) {
                  res2 <- yy
                  res1 <- apply(t(2:Max_K_Poly), 2, function(z) {
                    if (grepl(paste0("_P", z, "$"), yy)) {
                      if (!grepl("ns", modeltype)) {
                        res1 <- paste0("poly(",
                                       gsub(paste0("_P", z, "$"), "", yy),
                                       ",", z, ")")
                      } else {
                        res1 <- paste0("ns(",
                                       gsub(paste0("_P", z, "$"), "", yy),
                                       ", df=", z, ")")
                      }
                    } else {
                      res1 <- NA
                    }
                    res1
                  })
                  if (sum(is.na(res1)) != length(res1)) {
                    res2 <- res1[which(!is.na(res1))]
                  }
                  res2
                })
                res <- paste(y2, collapse = ":")
              }
              res
            }))
          }))

        if (nrow(test_factor_int[[nb]]) == 1 |
            length(nrow(test_factor_int[[nb]])) == 0) {
          test_factor_int[[nb]] <- t(test_factor_int[[nb]])
        }

        if (nb == 2) {
          l2 <- length(test_factors)
          test_factors[[Max_nb_Var + 1]] <- cbind(test_factors[[nb]],
                                                  test_factor_int[[nb]])
        }
        if (nb >= 3)
        {
          assign(paste0("l", nb), length(test_factors) + 1)
          if (length(nrow(test_factor_int[[nb]])) != 0) {
            # 1 interaction
            for (int_n in 1:ncol(test_factor_int[[nb]])) {
              test_factors[[length(test_factors) + 1]] <- cbind(test_factors[[nb]],
                                                                test_factor_int[[nb]][, int_n])
            }
            # 2 interactions
            comb2 <- utils::combn(ncol(test_factor_int[[nb]]), 2)
            for (int_n in 1:ncol(comb2)) {
              test_factors[[length(test_factors) + 1]] <- cbind(test_factors[[nb]],
                                                                test_factor_int[[nb]][, c(comb2[1, int_n])],
                                                                test_factor_int[[nb]][, c(comb2[2, int_n])])
            }
            # 3 interactions
            comb3 <- utils::combn(ncol(test_factor_int[[nb]]), 3)
            for (int_n in 1:ncol(comb3)) {
              test_factors[[length(test_factors) + 1]] <- cbind(test_factors[[nb]],
                                                                test_factor_int[[nb]][, c(comb3[1, int_n])],
                                                                test_factor_int[[nb]][, c(comb3[2, int_n])],
                                                                test_factor_int[[nb]][, c(comb3[3, int_n])])
            }
          } else {
            print("Only one factor combination left before Interactions")
            # 1 interaction
            for (int_n in 1:ncol(test_factor_int[[nb]])) {
              test_factors[[length(test_factors) + 1]] <- t(cbind(test_factors[[nb]],
                                                                  test_factor_int[[nb]][, int_n]))
            }
            # 2 interactions
            comb2 <- utils::combn(ncol(test_factor_int[[nb]]), 2)
            for (int_n in 1:ncol(comb2)) {
              test_factors[[length(test_factors) + 1]] <- t(cbind(test_factors[[nb]],
                                                                  test_factor_int[[nb]][, c(comb2[1, int_n])],
                                                                  test_factor_int[[nb]][, c(comb2[2, int_n])]))
            }
            # 3 interactions
            comb3 <- utils::combn(ncol(test_factor_int[[nb]]), 3)
            for (int_n in 1:ncol(comb3)) {
              test_factors[[length(test_factors) + 1]] <- t(cbind(test_factors[[nb]],
                                                                  test_factor_int[[nb]][, c(comb3[1, int_n])],
                                                                  test_factor_int[[nb]][, c(comb3[2, int_n])],
                                                                  test_factor_int[[nb]][, c(comb3[3, int_n])]))
            }
          }
        }  # end of nb>=4

      }  # end of Interaction


      # Save number of models for each step of the iterative procedure
      nb_all <- nb
      if (Interaction) {
        if (nb == 2) {
          nb_all <- c(nb, length(test_factors))
        }
        if (nb >= 3) {
          nb_all <- c(nb, get(paste("l", nb, sep = "")):length(test_factors))
        }
      }

      # Write models Y ~ param1 + param2 + ...  nb2 <- nb_all[1] x <-
      # test_factors[[nb2]][120,]
      for (nb2 in nb_all) {
        test_models_tmp[[nb2]] <- apply(
          test_factors[[nb2]], 1,
          function(x) {
            if (grepl("GLM", modeltype)) {
              x2 <- apply(t(x), 2, function(y) {
                res <- y
                res1 <- apply(t(2:Max_K_Poly), 2, function(z) {
                  if (grepl(paste0("_P", z, "$"), y)) {
                    if (!grepl("ns", modeltype)) {
                      res1 <- paste0("poly(",
                                     gsub(paste0("_P", z, "$"), "", y),
                                     ",", z, ")")
                    } else {
                      res1 <- paste0("ns(",
                                     gsub(paste0("_P", z, "$"), "", y),
                                     ", df=", z, ")")
                    }
                  } else {
                    res1 <- NA
                  }
                  res1
                })
                if (sum(is.na(res1)) != length(res1)) {
                  res <- res1[which(!is.na(res1))]
                }
                res
              })
            } else {
              x2 <- paste("s(", x, ",k=", Max_K, ")", sep = "")
            }

            # if models are factors, no smooth parameter
            testFact <- unique(c(grep("factor_", x, fixed = TRUE),
                                 grep("te(", x, fixed = TRUE), grep(":", x)))
            if (sum(testFact) >= 1) {
              x2[testFact] <- paste(x[testFact])
            }
            if (grepl("Log", modeltype)) {
              res2 <- paste("log(dataY) ~", paste(x2, collapse = " + "))
            } else {
              res2 <- paste("dataY ~", paste(x2, collapse = " + "))
            }
            if (grepl("KrigeGLM", modeltype)) {
              res2 <- paste("~", paste(x2, collapse = " + "))
            }
            res2
          })

        # Carefull nb2 are included in nb <=> Models with interaction in same
        # list than models without
        if (nb2 == nb_all[1]) {
          test_models_A <- test_models_tmp[[nb2]]
        } else {
          test_models_A <- c(test_models_A, test_models_tmp[[nb2]])
        }
      }  # end of nb2 loop
      test_models_A <- matrix(test_models_A)
      Models_tmp_nb <- test_models_A
      if (verbose >= 1) {
        print(paste("last model tested:", utils::tail(test_models_A, 1)))
      }
      if (verbose >= 2) {
        print(paste("nb_models:", length(Models_tmp_nb)))
      }

      # __First selection on AIC ----
      # to avoid problems in crossV Randomize the
      # order of models tested to optimize parallelisation
      Goodorder <- 1:nrow(Models_tmp_nb)
      Randorder <- sample(Goodorder, length(Goodorder))

      resAIC_tmp <- snow::parApply(
        cl, t(Randorder), 2,
        function(x, Y_data_sample,
                 Models_tmp_nb, modeltype, fixXI, Y_data_sample_lcc, MaxDist,
                 Phi, model, fix.lambda, lambda, libPaths) {
        .libPaths(c(libPaths, .libPaths()))
        res <- AIC_indices(
          x, Y_data_sample = Y_data_sample, Models_tmp_nb = Models_tmp_nb,
          modeltype = modeltype, fixXI = fixXI, Y_data_sample_lcc = Y_data_sample_lcc,
          MaxDist = MaxDist, Phi = Phi, Model = model,
          fix.lambda = fix.lambda, lambda = lambda)
        res
      }, Y_data_sample = Y_data_sample, Models_tmp_nb = Models_tmp_nb,
      modeltype = modeltype, fixXI = fixXI, Y_data_sample_lcc = Y_data_sample_lcc,
      MaxDist = MaxDist, Phi = Phi, model = Model[n_datac],
      fix.lambda = fix.Lambda[n_datac], lambda = Lambda[n_datac],
      libPaths = .libPaths())

      if (grepl("TweedGLM", modeltype) | grepl("KrigeGLM", modeltype)) {
        resAIC[[nb]] <- resAIC_tmp[1, order(Randorder)]
        resParam[[nb]] <- resAIC_tmp[2, order(Randorder)]
      } else {
        resAIC[[nb]] <- resAIC_tmp[order(Randorder)]
      }

      # stop loop if there is no remaining AIC
      if (sum(!is.na(resAIC[[nb]])) == 0) {
        print("no more models to test")
        break
      }

      # Remove models not fitted
      modelkeep <- which(!is.na(resAIC[[nb]]))

      Models_tmp_nb2 <- t(t(Models_tmp_nb[modelkeep, ]))

      # It is assumed that if a model did not converged with one covariate,
      # it cannot converged with additional ones In resAIC, non convergence
      # is parameterized as NaN Here we save covariates that will be removed
      # for further steps if convergence has not been reached for n=1 only
      # if (nb == 1) {
      #   NaN.covar <- which(is.nan(resAIC[[1]]))
      # }


      test_models[[nb]] <- Models_tmp_nb2
      resAIC_save[[nb]] <- resAIC[[nb]][modelkeep]
      if (grepl("TweedGLM|KrigeGLM", modeltype)) {
        resParam_save[[nb]] <- resParam[[nb]][modelkeep]
      } else {
        resParam_save <- NULL
      }
      # Save factors used
      test_factors_all_A <- test_factors[[nb]]
      if (length(nb_all) >= 2) {
        for (nb2 in nb_all[-1]) {
          test_factors_all_A <- rbind(test_factors_all_A, test_factors[[nb2]][,
                                                                              1:nb])
        }
      }
      test_factors_all[[nb]] <- test_factors_all_A[modelkeep, ]
      if (is.null(nrow(test_factors_all[[nb]])) & nb != 1) {
        test_factors_all[[nb]] <- t(test_factors_all[[nb]])
      } else if (nb == 1) {
        test_factors_all[[nb]] <- matrix(test_factors_all[[nb]])
      }


      if (verbose >= 2) {
        print(paste("nb_models w/o AIC=NA:", length(Models_tmp_nb2)))
      }
      if (verbose >= 1) {
        print(paste0("modeltype=", modeltypes[n_datac], ": ", n_datac,
                     "/", length(modeltypes), " ; nb=", nb, "/", Max_nb_Var))
      }
      # End of first selection on AIC

      # __Cross-validation ----
      # Parallelise on MC: One MC is for testing all models on one random sub-sample
      Output_models_A <- snow::parApply(
        cl, t(1:nbMC), 2,
        function(MC,
                 formulas, modeltype, saveAlea, Y_data_sample,
                 seqthd, resParam_save, Y_data_sample_lcc, #nb,
                 MaxDist, Phi, model, lambda, libPaths) {
          .libPaths(c(libPaths, .libPaths()))
          res <- crossV_indices(
            MC = MC, formulas = formulas,
            modeltype = modeltype, saveAlea = saveAlea, Y_data_sample = Y_data_sample,
            seqthd = seqthd,
            resParam_save = resParam_save, Y_data_sample_lcc = Y_data_sample_lcc, #, nb = nb
            MaxDist = MaxDist, Phi = Phi, model = model,
            lambda = lambda)
          c(res)
        }, formulas = Models_tmp_nb2, modeltype = modeltype,
        saveAlea = saveAlea, Y_data_sample = Y_data_sample,
        seqthd = seqthd, resParam_save = resParam_save[[nb]],
        Y_data_sample_lcc = Y_data_sample_lcc, MaxDist = MaxDist, #nb = nb,
        Phi = Phi, model = Model[n_datac],
        lambda = Lambda[n_datac], libPaths = .libPaths())

      # __end of clusterApply ---------------------
      if (verbose >= 2) {
        print("End of parApply")
      }

      # Re-arrange results in a array
      if (grepl("PA|TweedGLM", modeltype)) {
        Output_models[[nb]] <- array(
          Output_models_A,
          dim = c(7 + length(seqthd) + 3, length(Models_tmp_nb2), nbMC))
      } else {
        Output_models[[nb]] <- array(Output_models_A,
                                     dim = c(8, length(Models_tmp_nb2), nbMC))
      }

      # Estimation of the best threshold (thd) for each model.
      # Choose the threshold value engendering the
      # smallest difference Selectivity-Specificity
      # and being the closest to 0.5
      if (grepl("PA|TweedGLM", modeltype)) {
        BestTHD_crossV[[nb]] <- rep(NA, dim(Output_models[[nb]])[2])
        DiffSelSpeTHD_crossV[[nb]] <- matrix(NA,
                                             nrow = dim(Output_models[[nb]])[2],
                                             ncol = nbMC)
        for (Nmodel in 1:(dim(Output_models[[nb]])[2])) {
          DiffSelSpe_crossVMean <- apply(
            Output_models[[nb]][8:(length(seqthd) + 7), Nmodel, ], 1,
            function(x) mean(x, na.rm = TRUE))
          BestTHD_crossV[[nb]][Nmodel] <- mean(seqthd[
            which(DiffSelSpe_crossVMean ==
                    min(DiffSelSpe_crossVMean, na.rm = TRUE))])
          position <- NA
          # All having the least difference
          try(position <- which(
            abs(seqthd - BestTHD_crossV[[nb]][Nmodel]) ==
              min(abs(seqthd - BestTHD_crossV[[nb]][Nmodel]), na.rm = TRUE)
          ), silent = TRUE)
          DiffSelSpeTHD_crossV[[nb]][Nmodel, ] <- NA
          # The one being the closest to 0.5
          if (length(position) > 1) {
            position <- position[which(abs(seqthd[position] - 0.5) ==
                                         min(abs(seqthd[position] - 0.5)))][1]
          }
          # Errors in calculation
          if (MPIcalc) {
            if (length(position) == 1) {
              if (is.na(position) == TRUE) {
                print(paste("NA - Nmodel=", Nmodel, "Errors=",
                            Output_models[[nb]][8:(length(seqthd) + 7), Nmodel,
                                                ]))
              }
            }
            if (length(position) == 0) {
              position <- NA
              print(paste("0 - Nmodel=", Nmodel, "Errors=", Output_models[[nb]][8:(length(seqthd) +
                                                                                     7), Nmodel, ]))
            }
          }
          # Difference Selectivity-Specificity to choose Best THD
          try(DiffSelSpeTHD_crossV[[nb]][Nmodel, ] <- Output_models[[nb]][7 +
                                                                            position, Nmodel, ], silent = TRUE)
        }  # End of Nmodel
      }


      # __Choose the best models in terms of rank ----
      # Best: 5=CV; 6=RMSE; 7=corPearson
      if (grepl("PA", modeltype)) { # Best AUC
        IndexForRank_crossV <- -Output_models[[nb]][7 + length(seqthd) + 1,,]
      } else if (grepl("TweedGLM", modeltype)) { # Best RMSE
        IndexForRank_crossV <- Output_models[[nb]][7 + length(seqthd) + 2,,]
      } else {# Best RMSE
        IndexForRank_crossV <- Output_models[[nb]][6, , ]
      }

      if (is.null(dim(IndexForRank_crossV))) {IndexForRank_crossV <- t(IndexForRank_crossV)}

      # For computing simplifications
      if (length(which(IndexForRank_crossV == "Inf")) > 0) {
        IndexForRank_crossV[which(IndexForRank_crossV == "Inf")] <- 1e+15
      }
      if (length(which(IndexForRank_crossV == "-Inf")) > 0) {
        IndexForRank_crossV[which(IndexForRank_crossV == "-Inf")] <- -1e+15
      }

      IndexForRank_crossVMean <- rowMeans(IndexForRank_crossV, na.rm = T)

      # A weight is attributed to the test depending on the number of data used
      # in the validation dataset. With k-fold, the number of data may be slightly
      # different among folds. Moreover, factor covariates also unbalance the
      # repartition to assure at least one data for each level.
      best.distri <- best_distri(x = IndexForRank_crossV, w = saveAlea.nb,
                                 test  = "wilcoxon",
                                 na.max = 0.5, p.min = lim_pvalue,
                                 silent = TRUE, cl = cl)

      Line_keep <- best.distri$orderModels[which(best.distri$p.min.test)]

      if (graph) {
        graphics::par(xaxs = "i", mfrow = c(1,2))
        graphics::boxplot(t(IndexForRank_crossV),
                border = c("black", "forestgreen")[
                  1 + best.distri$p.min.test[order(best.distri$orderModels)]],
                pch = 20
        )
        graphics::points(IndexForRank_crossVMean, pch=20, col = 'red', cex = 0.5)

        Diff <- apply(IndexForRank_crossV, 1, function(x) {
          x - IndexForRank_crossV[best.distri$orderModels[1],]})
        graphics::boxplot(Diff,
                ylim = c(min(Diff, na.rm = TRUE),
                         quantile(c(Diff), probs = 0.95, na.rm = TRUE)),
                pch = 20,
                border = c("black", "forestgreen")[1 + best.distri$p.min.test[order(best.distri$orderModels)]],
                col = c("grey", "green")[1 + best.distri$p.min.test[order(best.distri$orderModels)]])
        graphics::points(best.distri$orderModels[1], 0, pch = 20, col = "red", cex = 2)
      }


      # __Save models in n_keep ---------------------
      # If there are
      # interactions, only covariates are stored (not interactions). All
      # possible interactions will then be added again in the next iteration
      # step.
      test_factors_all_nb <- test_factors_all[[nb]]
      for (i in 1:nb) {
        test_factors_all_nb[, i] <- apply(t(test_factors_all[[nb]][,i]), 2,
                                          function(x) {which(allParamN %in% x)})
      }
      # Remove duplicates in n_keep. Only possible if Interaction is TRUE
      if (nb == 1) {
        n_keep[[nb]] <- t(t(test_factors_all_nb[Line_keep, ]))
      } else {
        if (length(Line_keep) == 1) {
          n_keep[[nb]] <- t(test_factors_all_nb[Line_keep, ])
        } else {
          n_keep_A <- test_factors_all_nb[Line_keep, ]
          noDup <- which(
            !duplicated(
              apply(n_keep_A, 1,
                    function(x) paste(x, collapse = "-"))
            ))
          Line_keep <- Line_keep[noDup]
          n_keep[[nb]] <- n_keep_A[noDup,]
        }
      }
      if (is.null(nrow(n_keep[[nb]]))) {
        n_keep[[nb]] <- t(n_keep[[nb]])
      }

      # if not enough models, add the other good ones
      if (nrow(n_keep[[nb]]) < MinNbModel & length(best.distri$orderModels) != 1) {
        allnoDup <- which(!duplicated(
          apply(test_factors_all_nb, 1,
                function(x) paste(x, collapse = "-"))[
                  best.distri$orderModels]))[
                    1:min(length(best.distri$orderModels), MinNbModel)]
        Line_keep <- best.distri$orderModels[allnoDup]
        if (nb == 1) {
          n_keep_B <- t(t(test_factors_all_nb[
            best.distri$orderModels[allnoDup],
            ]))
        } else {
          n_keep_B <- test_factors_all_nb[
            best.distri$orderModels[allnoDup], ]
        }
        n_keep[[nb]] <- n_keep_B
      }

      if (verbose >= 1) {
        print(paste("n_keep with ", nb, "covariates (top 10)"))
        if (nrow(n_keep[[nb]]) == 1) {
          allParamN[as.numeric(as.character(n_keep[[nb]]))]
        } else {
          if (nb == 1) {
            print(matrix(apply(
              t(n_keep[[nb]][1:min(nrow(n_keep[[nb]]), 10),
                             ]), 2, function(y)
                               allParamN[as.numeric(as.character(y))])))
          } else {
            print(apply(
              n_keep[[nb]][1:min(nrow(n_keep[[nb]]), 10),
                           ], 2, function(y)
                             allParamN[as.numeric(as.character(y))]))
          }
        }
      }

      # __Save outputs at each iteration ----
      # Maybe usefull if procedure stops in the middle of the loop
      # save(Allinfo, file = paste0(saveWD, "/AllinfoM", n_datac, ".RData"))
      # save(saveAlea, file = paste0(saveWD, "/saveAleaM", n_datac, ".RData"))
      readr::write_rds(saveAlea, paste0(saveWD, "/saveAleaM", n_datac, ".rds"))
      # save(saveAlea.nb, file = paste0(saveWD, "/saveAlea.nbM", n_datac, ".RData"))
      readr::write_rds(saveAlea.nb, paste0(saveWD, "/saveAlea.nbM", n_datac, ".rds"))

      # /!\ Only save results of best models at each iteration /!\
      resAIC_save[[nb]] <- resAIC_save[[nb]][Line_keep]
      if (grepl("TweedGLM|KrigeGLM", modeltype)) {
        resParam_save <- resParam_save[Line_keep]
      }
      test_models[[nb]] <- as.matrix(as.data.frame(test_models[[nb]][Line_keep,]))
      test_factors_all[[nb]] <- as.matrix(data.frame(test_factors_all[[nb]][Line_keep,]))
      Output_models[[nb]] <- Output_models[[nb]][,Line_keep,]
      if (grepl("PA|TweedGLM", modeltype)) {
        BestTHD_crossV[[nb]] <- BestTHD_crossV[[nb]][Line_keep]
        DiffSelSpeTHD_crossV[[nb]] <- DiffSelSpeTHD_crossV[[nb]][Line_keep,]
      }

      readr::write_rds(
        resAIC_save,
        paste0(saveWD, "/resAIC_saveM", n_datac, ".rds"))
      # save(resAIC_save, file = paste0(saveWD, "/resAIC_saveM", n_datac, ".RData"))
      if (grepl("TweedGLM|KrigeGLM", modeltype)) {
        readr::write_rds(resParam_save,
                         paste0(saveWD, "/resParam_saveM", n_datac, ".rds"))
        #save(resParam_save,
        #     file = paste0(saveWD, "/resParam_saveM", n_datac, ".RData"))
      }
      # save(test_models, file = paste0(saveWD, "/test_modelsM", n_datac, ".RData"))
      readr::write_rds(test_models, paste0(saveWD, "/test_modelsM", n_datac, ".rds"))
      # No need to save
      #save(test_factors_all, file = paste0(saveWD, "/test_factors_allM", n_datac,
      #                                     ".RData"))
      #save(Output_models,
      #     file = paste0(saveWD, "/Output_modelsM", n_datac, ".RData"))
      readr::write_rds(Output_models,
                       paste0(saveWD, "/Output_modelsM", n_datac, ".rds"))
      # No need to save
      # save(n_keep, file = paste0(saveWD, "/n_keepM", n_datac, ".RData"))
      # save(test_name_list,
      # file = paste0(saveWD, "/test_name_listM", n_datac, ".RData"))

      if (grepl("PA|TweedGLM", modeltype)) {
        # save(BestTHD_crossV,
        #      file = paste0(saveWD, "/BestTHD_crossVM", n_datac, ".RData"))
        readr::write_rds(BestTHD_crossV,
                         paste0(saveWD, "/BestTHD_crossVM", n_datac, ".rds"))
        # save(DiffSelSpeTHD_crossV,
        #      file = paste0(saveWD, "/DiffSelSpeTHD_crossVM", n_datac, ".RData"))
        readr::write_rds(DiffSelSpeTHD_crossV,
                         paste0(saveWD, "/DiffSelSpeTHD_crossVM", n_datac, ".rds"))
      }



    }  #) # End of 'nb' loop
  }  # end of modeltypes

  if (!MPIcalc) {
    parallel::stopCluster(cl)
  }

  # Zip all outputs in a unique file.
  if (zip.file != FALSE) {
    if (zip.file == TRUE) {zip.file <- paste0(saveWD, ".zip")}
    # Options j allows to save directly files without path
    utils::zip(zip.file, files = saveWD, flags = "-urj9X", extras = "",
        zip = Sys.getenv("R_ZIPCMD", "zip"))
    # Return the file path of the zip file
    return(zip.file)
  } else {
    # Return the folder path where all results have been saved
    return(saveWD)
  }

  ## END of MPIcalc
}  # end of function

