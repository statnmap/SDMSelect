#' Modelselect configuration file
#'
#' @param RESET reset options to default
#' @param LOCAL Logical. set to true to use modified options in a local environment but not in global.
#' @param READ.ONLY Logical. To not modify the options
#' @param ADD Logical. New options can be added after the option function is created by explicitely specifying ADD = TRUE
#' @param ... Other options as listed below
#'
#' @return A list with all options
#'
#' @details
#' The script is not a "everything is working fine by default" script.
#' Different options have to be defined along the script.
#' Users are advised to run the code by little parts to avoid any miss-specifications of parameters.
#' Different types of distribution are available as well as different criteria of selection.
#' Users are advised to verify if options match their dataset
#' \describe{
#'   \item{MPIcalc}{set to TRUE if you are running the functions on a MPI cluster}
#'   \item{nbclust}{number of cores to use. Default to all cores available.}
#'   \item{graph}{Show all figures during calculation process (Recommend FALSE)}
#'   \item{datatype}{data type. Default to ContPosNull. See below}
#'   \item{modeltypes}{All model types to compare depending on datatype.}
#'   \item{MaxDist}{Maximum distance for kriging (For datatype = KrigeGLM)}
#'   \item{Phi}{Vector of length 2 with min and max for fitting Phi (For datatype = KrigeGLM)}
#'   \item{Model}{Vector of length of modeltypes with variogram model as used with
#'  \code{\link[geoR]{cov.spatial}} (For datatype = KrigeGLM)}
#'   \item{lcc_proj}{Projection in meters. Default to Lambert93.(For datatype = KrigeGLM)}
#'   \item{Lambda}{Vector of length of modeltypes with lambda corresponding to
#' box-cox transformation (For datatype = KrigeGLM).}
#'   \item{fix.Lambda}{Vector of length of modeltypes with logical, whether to fix lambda
#' or let the model chose it (For datatype = KrigeGLM).}
#'   \item{Max_nb_Var}{Maximum number of covariates kept in a model}
#'   \item{Max_K Maximum}{degree of freedom for simple variables in GAM models.
#' Default Max_K=5 similar to Max_K_Poly=4}
#'   \item{Max_K_te}{Maximum degree of freedom for tensor interactions in GAM models}
#'   \item{Max_K_Poly}{Maximum degree of polynom for simple variables in GLM models}
#'   \item{fixXI}{0=chosen by model (time consuming), 1=Poisson, 2=Gamma, 1<FixXI<2 compound poisson
#' (For datatype = TweedGLM). See \code{\link[statmod]{tweedie}}}
#'   \item{Interaction}{Logical. Whether to test for covariates interactions in models.
#' A maximum of three interactions will be tested in the same model. This may be
#' highly time consuming. Should be FALSE with any "KrigeGLM*".}
#'   \item{MinNbModel}{minimum number of models kept at each iteration to avoid removing not to bad models}
#'   \item{k_fold Numeric}{Value k of the k-fold cross-validation}
#'   \item{N_k_fold}{Numeric Value N of the N times k-fold cross-validation}
#'   \item{nbMC}{N_k_fold * k_fold (for compatibility)}
#'   \item{seqthd}{Sequence of thresholds tested to cut between 0 and 1 for PA data.}
#'   \item{lim_pvalue}{p-value limit for models retained in each of the "modeltype"
#'   cross-validation stepwise approach.
#'   Procedure seek for models having ranks not
#'   significantly different (p>=lim_pvalue) than the best model .
#'   "lim_pvalue" can be small for each step to keep a little more models than necessary.
#'   If p-value is not significant, models are considered with similar power of prediction,
#'   and thus kept in the following iteration.
#'   The smallest the p-value, the less discriminant the test, thus the highest number
#'   of models retained}
#'   \item{lim_pvalue_final}{Similar to lim_pvalue but used to select best model
#'   among all models at the very end of the procedure. This value may be higher
#'   than lim_pvalue to be more discriminant. This allows to select the best model
#'   among all models that have been fitted after the different cross-validations}
#'   \item{Y.max}{value to which multiplying the maximum value of data observations.
#'   This is used as a maximum predicted value for uncertainty calculation
#'   with Tweedie model as calculation of too high values maybe very long.
#'   This also fixes the maximum prediction for species distribution mapping when
#'   the model has a high uncertainty for sparse positions and lead to impossibly
#'   too high values.}
#'   \item{seed}{Numeric Seed for random number generation (Allow reproducibility
#'    between simulations)}
#'   \item{seqthd}{This allows to choose the best threshold value to predict
#'   presence-absence with probability of presence.
#' Balance between 0 and 1 in data may deviate the best threshold from value 0.5.
#' The best threshold is calculated among all cross-validations.
#' Threshold value chosen is the one the closest to specificity = sensitivity.}
#'
#' datatype options
#' \describe{
#'   \item{PA}{Presence-Absence data (=binomial distribution)}
#'   \item{Cont}{Continuous data. Positive and/or negative. (=gaussian distribution)}
#'   \item{ContPosNull}{Continuous but positive (or null) data like Biomass/Density data.
#'   For data with zero values,
#'   modeltypes allowing only positive values can also be tested (LogNormal, Gamma).
#'   In that case, a box-cox transformation will be applied: model is fitted on
#'   `log(X+1)`, but cross-validation compares `pred(Y)-1`, which is in the scale
#'   of the data. This allows comparison with Gaussian distribution.
#'   (=gaussian distribution; modeltypes are Gaussian, Gamma, Lognormal and Tweedie)}
#'   \item{Count}{model of Count data (=Poisson model)}
#'   \item{PosCount}{Poisson model on positive data (Difference with Count in goodness of fit)}
#'   \item{TweedGLM}{Tweedie distribution (only GLM)}
#'   \item{KrigeGLM}{Co-kriging model with gaussian distribution. Experimental -
#'   Equation of variogram have to be tested separately before any run with this script,
#'   verify that \code{\link{modelselect_opt}} is correctly set up:
#'   Default values are not safe at all !}
#' }
#'
#' modeltypes are types of model tested:
#' GLM, GLM with natural splines (GLMns), GAM.
#' If there is no "GLM" in the name of the modeltype,
#' then GAM is fitted with library mgcv.
#'
#' }
#'
#' @importFrom GlobalOptions set_opt .v
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Reset options
#' modelselect_opt(RESET = TRUE)
#' # List options
#' modelselect_opt()
#' # Show one option
#' modelselect_opt("modeltypes")
#' # Modify datatype
#' modelselect_opt$datatype <- "PA"
#' # Modify modeltypes tested
#' modelselect_opt("modeltypes") <- c("PA", "PAGLM")
#' }

modelselect_opt <- set_opt(
  MPIcalc = list(.value = FALSE),
  nbclust = list(.value = function() parallel::detectCores()),
  graph = list(.value = FALSE),
  datatype = list(.value = "ContPosNull"),
  modeltypes = list(.value = function()
    if (.v$datatype == "PA") {
      c("PAGLM", "PAGLMns", "PA")
    } else if (.v$datatype == "Cont") {
      c("ContGLM", "ContGLMns", "Cont")
    } else if (.v$datatype == "ContPosNull") {
      c("LogContGLM", "GammaGLM", "ContGLM",
        "LogContGLMns", "GammaGLMns", "ContGLMns",
        "LogCont", "Gamma", "TweedGLM", "TweedGLMns",
        "Cont")
    } else if (.v$datatype == "Count") {
      c("CountGLM", "CountGLMns", "Count")
    } else if (.v$datatype == "TweedGLM") {
      c("TweedGLMns", "TweedGLM")
    } else if (.v$datatype == "PosCount") {
      c("PosCountGLM", "PosCountGLMns", "PosCount")
    } else if (.v$datatype == "KrigeGLM") {
      c("KrigeGLM")
    } else if (.v$datatype == "KrigeGLM.dist") {
      c("TweedGLM", "KrigeGLM.dist.lambda", "KrigeGLM.dist")}),
  MaxDist = list(.value = function()
    ifelse(grepl("KrigeGLM", .v$datatype), 6000, NA)),
  Phi = list(.value = function()
    ifelse(grepl("KrigeGLM", .v$datatype), c(1000, .v$MaxDist), NA)),
  Model = list(.value = function()
    if (.v$datatype == "KrigeGLM") {"exponential"
    } else if (.v$datatype == "KrigeGLM.dist") {
      c(NA, "exponential", "exponential")
    } else {rep(NA, length(.v$modeltypes))}),
  Lambda = list(.value = function() if (.v$datatype == "KrigeGLM") {1
  } else if (.v$datatype == "KrigeGLM.dist") {
    c(NA, -1.5, 1)
  } else {rep(NA, length(.v$modeltypes))}),
  fix.Lambda = list(.value = function() if (.v$datatype == "KrigeGLM") {1
  } else if (.v$datatype == "KrigeGLM.dist") {
    c(TRUE, FALSE, TRUE)
  } else {rep(TRUE, length(.v$modeltypes))}),
  lcc_proj = list(.value = "+init=epsg:2154"),
  Max_nb_Var = list(.value = 7),
  Max_K = list(.value = 5),
  Max_K_te = list(.value = 3),
  Max_K_Poly = list(.value = 4),
  fixXI = list(.value = 0),
  Interaction = list(.value = function()
    ifelse(grepl("KrigeGLM", .v$datatype), FALSE, FALSE)),
  MinNbModel = list(.value = 4),
  k_fold = list(.value = 5),
  N_k_fold = list(.value = 20),
  nbMC = list(.value = function() .v$N_k_fold*.v$k_fold),
  seqthd = list(.value = seq(0.1, 0.9, 0.01)),
  lim_pvalue = list(.value = 0.001),
  lim_pvalue_final = list(.value = 0.005),
  Y.max = list(.value = function() ifelse(grepl("PA", .v$datatype), 1, 2)),
  seed = list(.value = 20)
)

#' Test internal modelopt
#'
test_opt <- function() {
  modelselect_opt$Max_nb_Var
}