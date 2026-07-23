## ================= ##
## R Code for MEI    ##
## ================= ##


library(lavaan)  ## load lavaan
library(MASS)  ## load MASS



# ==================== Creating Function "Full_MEI" ==================== #
#' Full Measurement Invariance Test
#'
#' Conduct configural invariance, full metric invariance, and full scalar invariance tests for data from independent groups.
#'
#' Reference for correct.cfi: Widaman, K. F., & Thompson, J. S. (2003). On specifying the null model for incremental fit indices in structural equation modeling. Psychological Methods, 8(1), 16-37.
#'
#' @param model User-specified CFA model.
#' @param data.source A data frame containing the observed variables used in the model.
#' @param Groups Grouping variable for cross-group comparisons.
#' @param Cluster Cluster variable for nested data. The Monte Carlo simulation method should be used for nested data.
#' @param correct.cfi If TRUE (default), baseline model will be corrected for calculation of incremental fit indices based on Widaman & Thompson (2003).
#'
#' @return lavaan outputs, model fit for configural, metric and scalar invariance models, summary of fit statistics for each model.
#' @export
#' @examples
#'
#' ## -- Example A: Measurement Invariance Test Across Groups -- ##
#'
#' # Data file is "Example.A"
#'
#' # Specify the measurement model - Model.A
#' Model.A <- '
#'        WorkLifeConflict =~ R45a + R45b + R45c + R45d + R45e
#'        Engagement =~ R90a + R90b + R90c
#'        Wellbeing =~ R87a + R87b + R87c + R87d + R87e
#' '
#'
#' ## ===== Full Measurement Invariance Test ===== ##
#' Full_MEI(Model.A, Example.A, Groups="Region")
#'
Full_MEI <- function(model, data.source, Groups, Cluster="NULL", correct.cfi=TRUE) {

  if (is.logical(correct.cfi) == FALSE) stop("correct.cfi (corrected baseline model for calculating incremental fit indices) can only be TRUE or FALSE")

  if (Cluster == "NULL") {
    ## Configural invariance model (Model.config) ##
    Model.config <<- lavaan::sem(model, data.source, group=Groups,
                                meanstructure=TRUE,
                                marker.int.zero = TRUE,
                                estimator="MLR",
                                missing = "ML",
                                information="observed")

    ## Metric invariance model (Model.metric) ##
    Model.metric <<- lavaan::sem(model, data.source, group=Groups,
                                meanstructure=TRUE,
                                marker.int.zero = TRUE,
                                estimator="MLR",
                                missing = "ML",
                                information="observed",
                                group.equal = c("loadings"))

    ## Scalar invariance model (Model.scalar) ##
    Model.scalar <<- lavaan::sem(model, data.source, group=Groups,
                                meanstructure=TRUE,
                                marker.int.zero = TRUE,
                                estimator="MLR",
                                missing = "ML",
                                information="observed",
                                group.equal = c("loadings", "intercepts"))
  } else {
    Model.config <<- lavaan::sem(model, data.source, group=Groups, cluster=Cluster,
                                meanstructure=TRUE,
                                marker.int.zero = TRUE,
                                estimator="MLR",
                                missing = "ML",
                                information="observed")

    ## Metric invariance model (Model.metric) ##
    Model.metric <<- lavaan::sem(model, data.source, group=Groups, cluster=Cluster,
                                meanstructure=TRUE,
                                marker.int.zero = TRUE,
                                estimator="MLR",
                                missing = "ML",
                                information="observed",
                                group.equal = c("loadings"))

    ## Scalar invariance model (Model.scalar) ##
    Model.scalar <<- lavaan::sem(model, data.source, group=Groups, cluster=Cluster,
                                meanstructure=TRUE,
                                marker.int.zero = TRUE,
                                estimator="MLR",
                                missing = "ML",
                                information="observed",
                                group.equal = c("loadings","intercepts"))
  } ## end if Cluster == "NULL"


  ## == Print Outputs == ##
  cat(rep("\n",3), "CONFIGURAL INVARIANCE MODEL", rep("\n", 2))
  print(lavaan::summary(Model.config, fit.measures=TRUE, standardized=TRUE, rsq=TRUE))
  cat(rep("\n",3), "METRIC INVARIANCE MODEL", rep("\n", 2))
  print(lavaan::summary(Model.metric, fit.measures=TRUE, standardized=TRUE, rsq=TRUE))
  cat(rep("\n", 3), "SCALAR INVARIANCE MODEL", rep("\n", 2))
  print(lavaan::summary(Model.scalar, fit.measures=TRUE, standardized=TRUE, rsq=TRUE))
  cat(rep("\n",3))

  ## == Compare fit indices across models == ##
#$  FitDiff <- compareFit(Model.config, Model.metric, Model.scalar, nested = TRUE)
#$  PFitDiff <- summary(FitDiff)
#$  print(PFitDiff)


  XX <- sapply(list(Model.config, Model.metric, Model.scalar), lavaan::fitMeasures, c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr", "aic", "bic"))

  if (correct.cfi == TRUE) {
    fitmeasures <- corrected.cfi(PSI.object=Model.scalar, data=data.source, Groups=Groups, cluster=Cluster)
    XX[,3] <- fitmeasures[c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr_bentler", "aic", "bic")]
  } # end (if correct.cfi)

  MODFIT <- matrix(1:27, nrow = 3, dimnames = list(c("Model.config", "Model.metric", "Model.scalar"), c("  chisq.scaled", "  df.scaled", "  pvalue.scaled", "  rmsea.scaled", "  cfi.scaled", "  tli.scaled", "  srmr", "  aic", "  bic")))
  MODFIT[1,1] <- format(round((XX[1,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[1,2] <- format(round((XX[2,1]), digits = 4), nsmall = 0, scientific = FALSE)
  MODFIT[1,3] <- format(round((XX[3,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[1,4] <- format(round((XX[4,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[1,5] <- format(round((XX[5,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[1,6] <- format(round((XX[6,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[1,7] <- format(round((XX[7,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[1,8] <- format(round((XX[8,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[1,9] <- format(round((XX[9,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,1] <- format(round((XX[1,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,2] <- format(round((XX[2,2]), digits = 4), nsmall = 0, scientific = FALSE)
  MODFIT[2,3] <- format(round((XX[3,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,4] <- format(round((XX[4,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,5] <- format(round((XX[5,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,6] <- format(round((XX[6,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,7] <- format(round((XX[7,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,8] <- format(round((XX[8,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,9] <- format(round((XX[9,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[3,1] <- format(round((XX[1,3]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[3,2] <- format(round((XX[2,3]), digits = 4), nsmall = 0, scientific = FALSE)
  MODFIT[3,3] <- format(round((XX[3,3]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[3,4] <- format(round((XX[4,3]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[3,5] <- format(round((XX[5,3]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[3,6] <- format(round((XX[6,3]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[3,7] <- format(round((XX[7,3]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[3,8] <- format(round((XX[8,3]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[3,9] <- format(round((XX[9,3]), digits = 4), nsmall = 4, scientific = FALSE)

  FITDIF <- matrix(1:6, nrow = 2, dimnames = list(c("Model.metric - Model.config", "Model.scalar - Model.metric"), c("  cfi.scaled", "  rmsea.scaled", "    srmr")))
  FITDIF[1,1] <- format(round((XX[5,2] - XX[5,1]), digits = 4), nsmall = 4, scientific = FALSE)
  FITDIF[1,2] <- format(round((XX[4,2] - XX[4,1]), digits = 4), nsmall = 4, scientific = FALSE)
  FITDIF[1,3] <- format(round((XX[7,2] - XX[7,1]), digits = 4), nsmall = 4, scientific = FALSE)
  FITDIF[2,1] <- format(round((XX[5,3] - XX[5,2]), digits = 4), nsmall = 4, scientific = FALSE)
  FITDIF[2,2] <- format(round((XX[4,3] - XX[4,2]), digits = 4), nsmall = 4, scientific = FALSE)
  FITDIF[2,3] <- format(round((XX[7,3] - XX[7,2]), digits = 4), nsmall = 4, scientific = FALSE)

  cat(rep("\n", 3))
  cat("####    MODEL FIT INDICES   ####", rep("\n", 2))
  print(MODFIT, quote=FALSE, right=TRUE)

  ## Compare fit indices across models ##
  print(lavaan::lavTestLRT(Model.config, Model.metric, Model.scalar))

  cat(rep("\n", 3))
  cat("####  DIFFERENCES IN FIT INDICES  ####", rep("\n", 2))
  print(FITDIF, quote=FALSE, right=TRUE)

  no.group <- lavInspect(Model.config, "ngroups")  # Number of groups #
  group.names <<- lavInspect(Model.config, "group.label")
  group.size <- lavaan::lavInspect(Model.config, "nobs")

  ## == Run Measurement Model for Each Group == ##
  Sub.fit.summary <- matrix(0, no.group, 11)
  colnames(Sub.fit.summary) <-
      c("Group", "Sample Size", "  chisq.scaled", "  df.scaled", "  pvalue.scaled", "  rmsea.scaled", "  cfi.scaled", "  tli.scaled", "  srmr", "  aic", "  bic")

  for (group.R in 1:no.group) {
    Sub.data <- subset(data.source, get(Groups) == group.names[group.R])
    if (Cluster == "NULL") {
      Sub.model <- lavaan::sem(model, Sub.data, marker.int.zero = TRUE, meanstructure=TRUE, estimator="MLR", missing = "ML", information="observed")
    } else {
      Sub.model <- lavaan::sem(model, Sub.data, marker.int.zero = TRUE, meanstructure=TRUE, cluster=Cluster, estimator="MLR", missing = "ML", information="observed")
    }
    Sub.fit.summary[group.R, 1] <- group.names[group.R]
    Sub.fit.summary[group.R, 2] <- group.size[group.R]

    Sub.fit.summary[group.R, 3:11] <- format(round(lavaan::fitMeasures(Sub.model, c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled",
        "tli.scaled", "srmr_bentler", "aic", "bic")), digits = 4), nsmall = 4, scientific = FALSE)
    Sub.fit.summary[group.R, 4] <- round(as.numeric(Sub.fit.summary[group.R, 4]), digits = 0)
  }  ## End (for group.R)

  cat(rep("\n",3), "##########  Model Fit Indices for each Group in Configural Invariance Model  ##########", rep("\n", 2))
  print(as.data.frame(Sub.fit.summary), right = TRUE, row.names=F, quote=F)
  cat("\n")

} ## end (Function Full_MEI)

# ==================== Finish Function "Full_MEI" ==================== #





# ==================== Create Function "CompareLoadings" ==================== #
#' Metric Invariance Test
#'
#' Conduct a full metric invariance test and identify a partial metric invariance model for data from independent groups.
#'
#' Reference: Measurement Equivalence/Invariance Test based on "Cheung, G. W. & Lau, R. S. (2012).  A direct comparison approach for testing measurement invariance.  Organizational Research Methods, 15, 167-198."
#'
#'
#' @param model User-specified CFA model.
#' @param data.source A data frame containing the observed variables used in the model.
#' @param Groups Grouping variable for cross-group comparisons.
#' @param Cluster Cluster variable for nested data. The Monte Carlo simulation method should be used for nested data.
#' @param Bootstrap Number of bootstrap samples, must be between 500 and 10,000. If not specified, 1 million Monte Carlo simulated samples (default) will be used.
#' @param Type1 Overall Type I error rate (default is 0.05) for identifying non-invariant items in the List and Delete method.
#' @param Type1Adj Adjustment of Type I error rate for multiple tests for each construct. Default is Possible False Discovery Rate (PFDR = Type I error rate/No. of freely estimated factor loadings/intercepts across all groups), can also use (BON)ferroni adjustment (Type I error rate/No. of pairwise comparisons), or "NULL".
#' @param BMSC Best Model Selection Criterion used to select among models with the same number of estimated parameters. "ChiSquare" = minimum Chi-Square value, "CFI" = maximum CFI value, "RMSEA" = minimum RMSEA value, and "SRMR" = minimum SRMR value (default).
#'
#' @return Partial metric invariance model in the PMI.txt file.
#' @export
#' @examples
#'
#' ## -- Example A: Measurement Invariance Test Across Groups -- ##
#'
#' # Data file is "Example.A"
#'
#' # Specify the measurement model - Model.A
#' Model.A <- '
#'        WorkLifeConflict =~ R45a + R45b + R45c + R45d + R45e
#'        Engagement =~ R90a + R90b + R90c
#'        Wellbeing =~ R87a + R87b + R87c + R87d + R87e
#' '
#'
#' ## Not run:
#' ## ===== Full Measurement Invariance Test ===== ##
#' Full_MEI(Model.A, Example.A, Groups="Region")
#' ## End (Not run)
#'
#' ## ===== Compare Loadings ===== ##
#' CompareLoadings(Model.A, Example.A, Groups="Region", Type1=0.05, Type1Adj="PFDR", BMSC="SRMR")
#'
CompareLoadings <- function(model, data.source, Groups, Cluster="NULL", Bootstrap=0, Type1=0.05, Type1Adj="PFDR", BMSC="SRMR") {

#  Bootstrap = 0 # Number of bootstrap samples
#  model = Model.A
#  data.source = Example.A
#  Groups = "Region"
#  Type1 = 0.05
#  Type1Adj = "PFDR"
#  Cluster="NULL"
#  BMSC="SRMR"

  Type1Adj <- toupper(Type1Adj)
  match.arg(Type1Adj, c("PFDR","BON", "NULL"))
  BMSC <- toupper(BMSC)
  match.arg(BMSC, c("CHISQUARE","CFI", "RMSEA", "SRMR"))

  arg1_char <- deparse(substitute(model))
  arg2_char <- deparse(substitute(data.source))
  arg3_char <- deparse(substitute(Groups))
  arg4_char <- deparse(substitute(Cluster))

  parsed <- lavParseModelString(model, as.data.frame = TRUE)
  parsed <- parsed[, c("lhs", "op", "rhs")]
  parsed <- parsed[parsed$op == "=~",]
  names.lv <- unique(parsed[,"lhs"])
  names.ind <- unique(parsed[, "rhs"])
  no.factor <- length(names.lv)
  no.items.g <- length(names.ind) # number of items per group
  no.items <- matrix(1:no.factor, nrow = 1)  # number of items per factor
  for (factor.no in 1:no.factor) { no.items[factor.no] <- sum(parsed[,"lhs"] == names.lv[factor.no] & parsed[,"op"] == "=~") } ## end loop factor.no

  names.item <- list()
  for (i in 1:no.factor) {
    names.item[[i]] <- parsed[which(parsed[,"lhs"] == names.lv[i] & parsed[,"op"] == "=~"), "rhs"]
  }

  ## Check for bootstrap sample number (Bootstrap) ##
  if (Bootstrap != 0) {
    b.no.integer <- Bootstrap == round(Bootstrap)
    if (b.no.integer == "FALSE") stop("Bootstrap sample number must be an integer")
    if (Bootstrap > 10000) stop("Bootstrap sample number greater than 10,000 is not recommended")
    if (Bootstrap < 500) stop("Bootstrap sample number smaller than 500 is not recommended")
    TYPE = "Bootstrap"
    b.no <- Bootstrap
  } else {
    TYPE = "MonteCarlo"
  } # end (Bootstrap != 0)


  ## ========== Run Configural Model ========== ##

  if (Cluster == "NULL") {
    Model.config <- lavaan::sem(model,
      data.source,
      group = Groups,
      missing = 'fiml',
      auto.fix.first = TRUE,
      marker.int.zero = TRUE,
      meanstructure = T,
      information = 'observed',
      estimator = 'MLR')
  } else {
      Model.config <- lavaan::sem(model,
      data.source,
      group = Groups,
      missing = 'fiml',
      auto.fix.first = TRUE,
      marker.int.zero = TRUE,
      meanstructure = T,
      information = 'observed',
      cluster=Cluster,
      estimator = 'MLR')
   }  # end Cluster

  ## Request summary outputs
  #$   print(lavaan::summary(Model.config, fit.measure = T, standardized = T, rsq = T))

  Model.config <<- Model.config

  ## ===== End (Run Configural Model) ===== ##


  par.est <- lavaan::coef(Model.config)  # sample parameters
  no.group <- lavInspect(Model.config, "ngroups")  # number of groups #
  group.names <<- lavaan::lavInspect(Model.config, "group.label")  # group names

  temp <- lavaan::parameterEstimates(Model.config)
  no.markers <- matrix(1:no.factor, nrow = 1)  # location of marker items
  for (factor.no in 1:no.factor) {
    no.markers[factor.no] <- which(temp[,"lhs"] == names.lv[factor.no] & temp[,"op"] == "=~" & temp[,"est"] == 1 & temp[,"group"] == 1)
    if (factor.no > 1) {no.markers[factor.no] <- no.markers[factor.no] - sum(no.items[1:(factor.no-1)])}
  } ## end for factor.no

  temp <- lavaan::parameterEstimates(Model.config, remove.nonfree=TRUE)
  simvcov <- lavInspect(Model.config, what="vcov")

  ## Extraxt factor loadings ##
  ext <- c(which(temp[,"op"] == "=~" & temp[,"group"] == 1))
  for (i in 2:no.group) {
    ext <- c(ext, which(temp[,"op"] == "=~" & temp[,"group"] == i))
  }
  par.est <- par.est[ext]
  simvcov <- simvcov[ext,ext]
  no.par.g <- length(par.est)/no.group  # number of estimated LX per group #


  cat(rep("\n", 3), "## ======= METRIC INVARIANCE ANALYSIS ======= ##", rep("\n", 2))  ## print heading

  # Calculate alpha for factor[factor.no] #
  ERate <- matrix(1:no.factor, nrow = 1)  # error rate for each factor

  cat("# -- Overall Type I Error Rates for Invariance Tests for Each Factor = ", Type1, " -- #", rep("\n", 2))

  if (Type1Adj == "PFDR") {  ## print adjusted error rate
    cat("# -- Adjusted Type I Error Rates for Invariance Tests (Adjusted by Possible False Discovery Rate) -- #", rep("\n", 2))
    for (factor.no in 1: no.factor) {
      no.adjust <- (no.items[factor.no]-1)*(no.group-1)
      ERate[factor.no] <- Type1 / no.adjust
      cat(paste0("Factor ", factor.no, ": Number of items = ", no.items[factor.no], "; number of free factor loadings = ", no.adjust, "; adjusted error rate = ", sprintf("%.4f", ERate[factor.no]), "\n"))
    } # end (for factor.no)
  } else if (Type1Adj == "BON") {
    cat("# -- Adjusted Type I Error Rates for Invariance Tests (Bonforroni adjustment for number of pairwise comparisons) -- #", rep("\n", 2))
    for (factor.no in 1: no.factor) {
      no.adjust <- (no.items[factor.no]-1)*(no.group*(no.group-1)/2)
      ERate[factor.no] <- Type1 / no.adjust
      cat(paste0("Factor ", factor.no, ": Number of items = ", no.items[factor.no], "; number of pairwise comparisons = ", no.adjust, "; adjusted error rate = ", sprintf("%.4f", ERate[factor.no]), "\n"))
    } # end (for factor.no)
  } else if (Type1Adj == "NULL") {
    cat("# -- Error Rates for Invariance Tests (No adjustment is made) -- #", rep("\n", 2))
    cat(paste0("Type I Error Rate for Invariance Tests = ", Type1, "\n"))
    for (factor.no in 1: no.factor) {
      ERate[factor.no] <- Type1
    } # end (for factor.no)
  } # end (Type1Adj)
  cat(rep("\n",2))

  if (TYPE == "MonteCarlo") {

    mcmc <- MASS::mvrnorm(n=2000000, mu=par.est, Sigma=simvcov, tol = 1e-6)  # Run 2,000,000 simulations
    bootcoef <- mcmc
    bootno <- nrow(mcmc)  # No. of successful simulated samples
    cat(paste0("Number of Successful Simulated Samples = ", bootno, "\n"))

  } else { # TYPE=Bootstrap

    ## == Simplified bootstrapping model == ##
    Config.boot <- lavaan::cfa(model, data = data.source, group=Groups,
                   meanstructure=TRUE,
                   auto.fix.first = TRUE,
                   marker.int.zero = TRUE,
                   ordered=FALSE,
                   missing = "fiml", se = "none", test="none", check.start = FALSE, check.post = FALSE, check.gradient = FALSE)

    ## == Bootstrapping == ##
    bootcoef <- lavaan::bootstrapLavaan(Config.boot, R = b.no, ncpus = max(1L, parallel::detectCores() - 1L), parallel="snow")

    # Remove error and nonadmissible bootstrap samples #
    B.na <- attr(bootcoef,"nonadmissible")
    B.er <- attr(bootcoef,"error.idx")
    B.de <- c(B.na,B.er)
    if (length(B.de) != 0) {
      B.re <- bootcoef[-c(B.de),]
      bootcoef <- B.re
    }
    bootno <- nrow(bootcoef)  # number of successful bootstrap samples
    cat(paste0("Number of Successful Bootstrapped Samples = ", bootno, "\n"))
    bootcoef <- bootcoef[,ext]

  } ## end MonteCarlo or Bootstrap


  ## -- Add group number to items in first group in par.est and bootcoef -- ##
  names(par.est)[1:no.par.g] <- paste0(names(par.est)[1:no.par.g], ".g1")
  colnames(bootcoef)[1:no.par.g] <- paste0(colnames(bootcoef)[1:no.par.g], ".g1")


  ## == Start the factor.no loop for CompareLoadings == ##
  for (factor.no in 1: no.factor) {

    alpha <- ERate[factor.no]
    flY <<- matrix(" ",1, (no.group+2))
    flYY <<- matrix(" ",1, (no.group+2))
    EP <<- 1  # estimated parameter number

    no.k <- no.items[factor.no]*no.group  # number of factor loadings in FL.PMI
    FL.PMI <- matrix(1:no.k, nrow = no.group)  # factor loading matrix
    colnames(FL.PMI) <- paste0(' Item ', 1:no.items[factor.no])
    rownames(FL.PMI) <- lavaan::lavInspect(Model.config, "group.label")

    for (i in 1:no.group) {
      for (j in 1:no.items[factor.no]) {
        if (j == no.markers[factor.no]) {
          FL.PMI[i, j] <- 1
        } else {
          FL.PMI[i, j] <- format(round(par.est[paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[j]],".g",i)], digits=4), scientific = FALSE)
        }  ## end if j
      }  ## end for j
    }  ## end for i
    class(FL.PMI) <- "numeric"

    for (Referent in 1:no.items[factor.no]) {  ## Loop referent item number for comparison

      ## ==  Print Factor Loadings of all groups  == ##
      no.k <- no.items[factor.no]*no.group  # number of factor loadings in FL
      FL <- matrix(1:no.k, nrow = no.group)  # factor loading matrix
      colnames(FL) <- paste0(' Item ', 1:no.items[factor.no])
      rownames(FL) <- lavaan::lavInspect(Model.config, "group.label")

      if (Referent == 1) {  ## if referent is the first item
        for (i in 1:no.group) {
          for (j in 1:no.items[factor.no]) {
            if (j == Referent) {
              FL[i, j] <- 1
            } else {
              FL[i, j] <- FL.PMI[i, j]
            }  ## end if j
          }  ## end loop j
        }  ## end loop i
      } else {  ## referent is not the first item
        for (i in 1:no.group) {
          for (j in 1:no.items[factor.no]) {
            if (j == Referent) {
              FL[i, j] <- 1
            } else {
              FL[i, j] <- FL.PMI[i, j]/FL.PMI[i, Referent]
            }  ## end if j
          }  ## end loop j
        }  ## end loop i
      }  ## end if referent

#$      cat(rep("\n", 3), "Factor Loadings ", rep("\n", 2))
#$      print(round(FL[], digits=4))

    ## == End (Print Factor Loadings)  == ##

      for (nArg in 1: no.items[factor.no]) {  ## Argument item number for comparison
        no.dif <- factorial(no.group) / factorial(no.group - 2) / 2  # No. of pairwise comparisons for each item
        boot.dif.lx <- matrix(0, nrow(bootcoef), no.dif)  # Create bootstrap difference matrix
        samp.dif.lx <- matrix(0, no.group, no.group) # Create sample difference matrix
        if (nArg == Referent) {
          next
        } else {
          Arg <- nArg
        } ## end if nArg


      ## == Calculate bootstrap difference and sample estimate difference == ##
        comp = 0
        for (r in 1:(no.group-1)) {  ## r is referent group
          for (a in (r+1):no.group) {  ## a is argument group
          comp = comp + 1
          if (Referent == 1) {
            boot.dif.lx[,comp] <- bootcoef[, paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Arg]],".g", r)] -
                                  bootcoef[, paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Arg]],".g", a)]
            samp.dif.lx[r,a] <- par.est[paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Arg]],".g", r)] -
                                par.est[paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Arg]],".g", a)]
          } else {
            if (Arg == no.markers[factor.no]) {
              boot.dif.lx[,comp] <- 1/bootcoef[, paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Referent]],".g", r)] -
                                    1/bootcoef[, paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Referent]],".g", a)]
              samp.dif.lx[r,a] <- 1/par.est[paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Referent]],".g", r)] -
                                  1/par.est[paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Referent]],".g", a)]
            } else { # Arg != no.markers[factor.no]
              boot.dif.lx[,comp] <- bootcoef[,paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Arg]],".g", r)]/
                                    bootcoef[,paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Referent]],".g", r)] -
                                    bootcoef[,paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Arg]],".g", a)]/
                                    bootcoef[,paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Referent]],".g", a)]
              samp.dif.lx[r,a] <- par.est[paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Arg]],".g", r)]/
                                  par.est[paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Referent]],".g", r)] -
                                  par.est[paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Arg]],".g", a)]/
                                  par.est[paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Referent]],".g", a)]
            }  ## end if Arg
          }  ## end if Referent
        }  ## end for a
      }  ## end for r

#$ print(samp.dif.lx)

      ## == Calculate Percentile Probability == ##

      F1.comp.pp <- matrix(0, no.group, no.group)  ## Percentile Probability
      pno_nipair <- 0  ## Set number of non-invariant pairs to zero
      comp = 0
      for (r in 1:(no.group-1)) {  ## r is the referent group
        for (a in (r+1):no.group) {  ## a is the argument
          comp = comp + 1
          if (quantile(boot.dif.lx[, comp], probs = 0.5, na.rm = TRUE) > 0) {
            F1.comp.pp[r,a] = 2*(sum(boot.dif.lx[, comp] < 0, na.rm=TRUE)/bootno)
          } else {
            F1.comp.pp[r,a] = 2*(sum(boot.dif.lx[, comp] > 0, na.rm=TRUE)/bootno)
          }  ## end if
        }  ## end loop a
      }  ## end loop r


      ## == Percentile Confidence Intervals of Pairwise Comparisons == ##

      PCI <- matrix(1:(no.dif*10), nrow = no.dif)
      colnames(PCI) <- c("Group 1","Group 2","0.5%","2.5%","5%","Estimate","95%","97.5%","99.5%","p-value")

      comp = 0
      for (r in 1:(no.group-1)) {
        for (a in (r+1):no.group) {
          comp = comp + 1
          PCI[comp, 1] <- r
          PCI[comp, 2] <- a
          PCI[comp, 3] <- quantile(boot.dif.lx[, comp],c(0.005), na.rm = TRUE)
          PCI[comp, 4] <- quantile(boot.dif.lx[, comp],c(0.025), na.rm = TRUE)
          PCI[comp, 5] <- quantile(boot.dif.lx[, comp],c(0.05), na.rm = TRUE)
          PCI[comp, 6] <- samp.dif.lx[r,a]
          PCI[comp, 7] <- quantile(boot.dif.lx[, comp],c(0.95), na.rm = TRUE)
          PCI[comp, 8] <- quantile(boot.dif.lx[, comp],c(0.975), na.rm = TRUE)
          PCI[comp, 9] <- quantile(boot.dif.lx[, comp],c(0.995), na.rm = TRUE)
          PCI[comp,10] <- F1.comp.pp[r,a]
          if (F1.comp.pp[r,a] < alpha) {
            if (pno_nipair == 0)	{
              pnipair <- c(r,a)
              pno_nipair <- pno_nipair + 1
            } else {
              pnipair <- c(pnipair,r,a)
              pno_nipair <- pno_nipair + 1
            }  ## end if pno_nipair
          }  ## end if F1.comp.pp
        }  ## end loop a
      }  ## end loop r

#$      cat("\n")
#$      if (Referent > Arg) {
#$        cat(paste0("Referent Item: ", Referent, "  Argument Item: ", Arg, "\n"))
#$      } else {
#$        cat(paste0("Referent Item: ", Referent, "  Argument Item: ", (Arg+1), "\n"))
#$      }  ## end if Referent
#$      cat("\n")
#$      print(round(PCI[], digits=4))

      ## == Run List and Delete == ##
#$      cat(rep("\n", 2), "Sets of Groups with Invariant Loadings", "\n")
        if (pno_nipair == no.dif) {

          flZ <<- matrix(" ", 1, (no.group+2))
          flX <<- matrix(1:(no.group+2), 1, (no.group+2))
          flX[1, 1:no.group] <- 0
          for (temp.flZ in 1: no.group) {
            flZ[1,temp.flZ] <- paste0("F", factor.no, "L", EP)
            EP <<- EP + 1
          }
          flX[1, (no.group+1)] <- Referent
          flZ[1, (no.group+1)] <- Referent
#?          if (Referent > Arg) {
            flX[1, (no.group+2)] <- Arg
            flZ[1, (no.group+2)] <- Arg
#?          } else {
#?            flX[1, (no.group+2)] <- Arg + 1
#?            flZ[1, (no.group+2)] <- Arg + 1
#?          }  ## end if Referent
#$        print(flX, quote=FALSE)
          flY <<- rbind(flY,flX)
#$        print(flZ, quote=FALSE)
          flYY <<- rbind(flYY,flZ)

        } else if (pno_nipair > 0) {

          flX <<- matrix(" ", 1, (no.group+2))
          flZ <<- matrix(" ", 1, (no.group+2))

          listanddelete.lx(factor.no, no.group, no_nipair=pno_nipair, nipair=pnipair, Referent, Arg) ## Run list and delete lx

        } else {  # (pno_nipair = 0)

          flX <- matrix(1:(no.group+2), 1, (no.group+2))
          flZ <<- matrix(" ", 1, (no.group+2))
          flZ <- matrix(paste0("F", factor.no, "L", EP), 1, (no.group+2))
          EP <<- EP + 1
          flX[1, (no.group+1)] <- Referent
          flZ[1, (no.group+1)] <- Referent
          flX[1, (no.group+2)] <- Arg
          flZ[1, (no.group+2)] <- Arg
          flY <<- rbind(flY,flX)
#$        print(flZ, quote=FALSE)
          flYY <<- rbind(flYY,flZ)

        }  ## end if pno_nipair
      } ## end loop Arg
    } ## end loop Referent


  if (length(unique(flY[, no.group+1])) == 2 & no.items[factor.no] == 2) { # Only 2 groups with more than 2 items
    flX <<- matrix(1:(no.group+2), 1, (no.group+2))
    flZ <<- matrix(" ", 1, (no.group+2))
    flYa <<- flY[nrow(flY), no.group+1]
    flYb <<- flY[nrow(flY), no.group+2]
    flX[1, no.group + 1] <- flYb
    flX[1, no.group + 2] <- flYa
    flY <<- rbind(flY,flX)

    flZ[1, no.group + 1] <- flYb
    flZ[1, no.group + 2] <- flYa
    for (temp.flYY in 1: no.group) {
      flZ[1, temp.flYY] <- paste0("F", factor.no, "I", EP)
      EP <<- EP + 1
    }
    flYY <<- rbind(flYY,flZ)
  }

  flY <<- flY[-c(1),]
#$  cat("\n")
#$  print(flY, quote = FALSE)

  flYY <<- flYY[-c(1),]
#$  cat("\n")
#$  print(flYY, quote = FALSE)

  if (is.matrix(flY) == FALSE) { flY <<- matrix(flY, nrow=1) }
  if (is.matrix(flYY) == FALSE) { flYY <<- matrix(flYY, nrow=1) }

  NIcombine <- table(flY[, no.group+1], flY[, no.group+2])
#  NI.level <- unique(flY[, no.group+1])
#  NIcombine <- NIcombine[NI.level, NI.level]

  tempR <- unique(flYY[, no.group+1])

    for (aR in 1:length(tempR)) {
      R <- tempR[aR]
      for (r in 1:no.items[factor.no]) {
        if (r == 1) {
          assign(paste0("Rset", r), c(length(which(flYY[,no.group+1] < R))+1:(NIcombine[R, r])))
        } else if (r == R) {
          assign(paste0("Rset", r), 0)
        } else {
          assign(paste0("Rset", r), c(length(which(flYY[,no.group+1] < R))+sum(NIcombine[R, 1:(r-1)])+1:(NIcombine[R, r])))
        }
        if (R == 1) { assign("Rset1", 0) }
      }  ## end loop r

      ## expand.grid(Rset1, Rset2, Rset3, Rset4, Rset5 ...)

      temp.row <- eval(parse(text = paste0("expand.grid(", paste0(c("Rset"), 1: no.items[factor.no], collapse=","),")")))
      if (aR == 1) {
        Set.row <- temp.row
      } else {
        Set.row <- rbind(Set.row, temp.row)
      } ## end if R
    }  ## end loop aR

    Model.load <- array(" ", dim = c(no.items[factor.no], no.group, nrow(Set.row)))  # all possible invariance models

    for (Xset in 1: nrow(Set.row)) {
      a <- matrix(" ", no.items[factor.no], no.group)
      for (r in 1:no.items[factor.no]) {
        if (Set.row[Xset, r] == 0) {
          a[r, ] <- 1
        } else {
          a[r, ] <- flYY[Set.row[Xset, r], 1:no.group]
        } ## end if Set.row
#?        if (flYY[r, no.group+1] > 1) { break }
      }  ## end loop r
      Model.load[, , Xset] <- a
    }  ## end loop Xset

#$ print(Model.load)

  ## == Model.summary == ##
  # Column 1 - number of estimated loadings
  # Column 2 - number of measurement invariance items
  # Column 3 - p-value of Model chi-square
  # Column 4 - CFI
  # Column 5 - RMSEA
  # Column 6 - SRMR

    Model.summary <- matrix(0, nrow(Set.row), 7)
    for (R in 1:nrow(Set.row)) {
      Model.summary[R, 1] <- length(unique(c(Model.load[,,R]))) - 1
      no.MI <- no.items[factor.no] - 1
      for (Ra in 1: no.items[factor.no]) {
        if (length(unique(c(Model.load[Ra,,R]))) > 1) {no.MI <- no.MI - 1}
      }  ## end loop Ra
      Model.summary[R, 2] <- no.MI
    }  ## end loop R


    ## == Recommended Models (Rec.Model) - Largest no. of MI items, then smallest number of estimated loadings == ##
    Rec.Model <- which(Model.summary[,2] == max(Model.summary[,2]))
    Rec.Model <- Rec.Model[which(Model.summary[Rec.Model,1] == min(Model.summary[Rec.Model,1]))]

    ## == Find recommended model if no.items > 2 == ##
    if (no.items[factor.no] > 2) {

      Rec.Model.load <- matrix(nrow = length(Rec.Model), ncol = length(lavaan::coef(Model.config)))
      MODFIT <- matrix(0, nrow = length(Rec.Model), ncol = 9)
      colnames(MODFIT) <-
        c("  chisq.scaled", "  df.scaled", "  pvalue.scaled", "  rmsea.scaled", "  cfi.scaled", "  tli.scaled", "  srmr_bentler", "  aic", "  bic")
      Rec.Model.FL <- array(0, dim = c(no.group, no.items[factor.no], length(Rec.Model)))
      colnames(Rec.Model.FL) <- paste0(" Item ", 1:no.items[factor.no])
      rownames(Rec.Model.FL) <- paste0(" Group ", lavInspect(Model.config, "group.label"), "   ")

      for (Model.R in 1: length(Rec.Model)) {

        ## == Specify the partial metric invariance models (PMI.X) == ##

        final.load <- unique(as.vector(Model.load[,, Rec.Model[Model.R]]))
        for (j in 1:no.group) {
          for (i in 1:no.items[factor.no]) {
            if (Model.load[i,j,Rec.Model[Model.R]] != 1) { Model.load[i,j, Rec.Model[Model.R]] <- paste0("F", factor.no, "L", which(Model.load[i,j,Rec.Model[Model.R]] == final.load)) }
          } # end for i
        } # end for j

        PMI.X.1 <- paste0("PMI.X <- '",
                          names.lv[factor.no], " =~ ", "c(", paste0(Model.load[1, , Rec.Model[Model.R]], collapse=","), ")*", names.item[[factor.no]][[1]], " + ")
        for (i in 2: (no.items[factor.no]-1)) {
          PMI.X.1 <- paste0(PMI.X.1, "c(",paste0(Model.load[i, , Rec.Model[Model.R]], collapse=","),")*", names.item[[factor.no]][[i]], " + ")
        } ## end loop i
        PMI.X.1 <- paste0(PMI.X.1, "c(",paste0(Model.load[no.items[factor.no], , Rec.Model[Model.R]], collapse=","),")*", names.item[[factor.no]][[i+1]], "'")

        eval(parse(text = PMI.X.1))

        ## == Run PMI.X == ##
        if (Cluster == "NULL") {
          eval(parse(text=   "PMI.X.fit <- lavaan::sem(PMI.X,
                             data.source,
                             group = Groups,
                             missing = 'fiml',
                             auto.fix.first = FALSE,
                             marker.int.zero = TRUE,
                             meanstructure = T,
                             information = 'observed',
                             estimator = 'MLR')"))
        } else {
          eval(parse(text=   "PMI.X.fit <- lavaan::sem(PMI.X,
                             data.source,
                             group = Groups,
                             missing = 'fiml',
                             auto.fix.first = FALSE,
                             marker.int.zero = TRUE,
                             meanstructure = T,
                             information = 'observed',
                             cluster = Cluster,
                             estimator = 'MLR')"))
        }


        ## == Request summary outputs == ##
#$      eval(parse(text = "print(lavaan::summary(PMI.X.fit, fit.measure = T, standardized = T, rsq = T))"))

        ## == Save fit indices == ##
        XX <- lavaan::fitMeasures(PMI.X.fit,
           c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr_bentler", "aic", "bic"))
        MODFIT[Model.R,] <- format(round((XX), digits = 4), nsmall = 4, scientific = FALSE)

        ## == Save Factor Loadings == ##
        for (G in 1: no.group) {
          Rec.Model.FL[G, , Model.R] <- eval(parse(text = paste0("lavaan::lavInspect(PMI.X.fit, what='est')$'", group.names[G], "'$lambda")))
        } ## end loop G

      } ## end loop Model.R


      ## == Print model fit indices == ##
      rownames(MODFIT) <- c(paste0("Model ", 1: length(Rec.Model)))

#$    cat(rep("\n",2), paste0("Possible Models for Factor Number ", factor.no), rep("\n",2))
#$    cat("####    MODEL FIT INDICES   ####", rep("\n", 2))
#$    print(MODFIT, quote=FALSE, right=TRUE)

    ## == Print factor loadings == ##
#$    for (Model.R in 1: length(Rec.Model)) {
#$      cat(rep("\n", 2), "Factor Loadings of Model", Model.R, rep("\n", 2))
#$      print(round(Rec.Model.FL[,,Model.R], digits=4))
#$    }  ## end loop Model.R

      ## == First model with best model fit (R.Model) == ##
      class(MODFIT) <- "numeric"
      if (BMSC == "Chisquare") {
        Model.R <- which.min(MODFIT[,1])
      } else if (BMSC == "CFI") {
        Model.R <- which.max(MODFIT[,5])
      } else if (BMSC == "RMSEA") {
        Model.R <- which.min(MODFIT[,4])
      } else if (BMSC == "SRMR") {
        Model.R <- which.min(MODFIT[,7])
      }
      Model.R <- Model.R[1]
      S.Model <- Model.load[,,Rec.Model[Model.R]]  # Recommended Model
      un.load <- sort(subset(unique(c(S.Model)), unique(c(S.Model))[] != 1))
      no.un.load <- length(un.load)
      R.Model <- S.Model
      for (j in 1:no.un.load) { R.Model[S.Model == un.load[j]] <- paste0("F", factor.no, "L", j) }  ## end loop j

    ## == Save Recommended Model (Recommend.Model) == ##
      if (factor.no == 1) {Recommend.Model <- matrix("PMI.Model.R <- '", 1)}
      Recommend.Model <-
        rbind(Recommend.Model, paste0("  ", names.lv[factor.no], " =~ ", "c(", paste0(R.Model[1,], collapse=","), ")*", names.item[[factor.no]][[1]], " + "))
      for (i in 2: (no.items[factor.no]-1)) {
        Recommend.Model <- rbind(Recommend.Model, paste0("  c(",paste0(R.Model[i,], collapse=","),")*",names.item[[factor.no]][[i]], " + "))
      }  ## end loop i
      Recommend.Model <-
        rbind(Recommend.Model, paste0("  c(",paste0(R.Model[no.items[factor.no], ], collapse=","),")*", names.item[[factor.no]][[i+1]]))
    }  ## end if no.items[factor.no] > 2

    if (no.items[factor.no] < 3) {   ## Begin Models with 2 items ##
      R.Model <- Model.load[,,1]  # Recommended Model
      un.load <- subset(unique(c(R.Model)), unique(c(R.Model))[] != 1)
      no.un.load <- length(un.load)
      for (j in 1:no.un.load) { R.Model[R.Model == un.load[j]] <- paste0("F", factor.no, "L", j) }  ## end loop j

      ## == Save Recommended Model (Recommend.Model) == ##
      if (factor.no == 1) {Recommend.Model <- matrix("PMI.Model.R <- '", 1)}
      Recommend.Model <-
        rbind(Recommend.Model, paste0("  ", names.lv[factor.no], " =~ ", "c(", paste0(R.Model[1,], collapse=","), ")*", names.item[[factor.no]][[1]], " + "))
      Recommend.Model <-
        rbind(Recommend.Model, paste0("  c(",paste0(R.Model[2, ], collapse=","),")*", names.item[[factor.no]][[2]]))
    }  ## End (Models with 2 items)

  } ## End (loop factor.no)


  Recommend.Model <- rbind(Recommend.Model, "  '")  ## last line of Recommended Model

  cat(rep("\n", 2), "## =====  Recommended Model  ===== ##", rep("\n", 2))
  for (i in 1: nrow(Recommend.Model)) { cat(Recommend.Model[i], "\n") } ## end loop i

  cat(rep("\n", 2))
  cat("## ===== lavaan Outputs of PMI.Model.R ===== ##")
  cat(rep("\n", 2))

  ## == Start print recommended model to PMI.txt (Partial Metric Invariance Model) == ##
  sink('PMI.txt')

  cat("## =====  Recommended Model  ===== ##")
  cat("\n")
  for (i in 1: nrow(Recommend.Model)) { cat(Recommend.Model[i], "\n") } ## end loop i

  cat(rep("\n",2), "## Run PMI.Model.R")
  cat(rep("\n",2), paste0("  PMI.Model.fit <- lavaan::sem(PMI.Model.R,"))
  cat("\n", paste0("    ", arg2_char, ","))
  cat("\n", paste0("    group = ", arg3_char,","))
  cat("\n", "    missing = 'fiml',")
  cat("\n", "    auto.fix.first = FALSE,")
  cat("\n", "    marker.int.zero = TRUE,")
  cat("\n", "    meanstructure = T,")
  cat("\n", "    information = 'observed',")
  if (Cluster != "NULL") { cat("\n", paste0("    cluster = ", arg4_char,",")) }
  cat("\n", "    estimator = 'MLR')")

  cat(rep("\n", 3), "## Request summary outputs")
  cat(rep("\n", 2), paste0("  print(lavaan::summary(PMI.Model.fit, fit.measure = T, standardized = T, rsq = T))"), "\n")

  sink()  ## Stop writing to file


  eval(parse(text = Recommend.Model))
  PMI.Model.R <<- PMI.Model.R
  ## == Run PMI.Model.R == ##
  if (Cluster == "NULL") {
    eval(parse(text=   "PMI.Model.fit <- lavaan::sem(PMI.Model.R,
                       data.source,
                       group = Groups,
                       missing = 'fiml',
                       auto.fix.first = FALSE,
                       marker.int.zero = TRUE,
                       meanstructure = T,
                       information = 'observed',
                       estimator = 'MLR')"))
  } else {
    eval(parse(text=   "PMI.Model.fit <- lavaan::sem(PMI.Model.R,
                       data.source,
                       group = Groups,
                       missing = 'fiml',
                       auto.fix.first = FALSE,
                       marker.int.zero = TRUE,
                       meanstructure = T,
                       information = 'observed',
                       cluster = Cluster,
                       estimator = 'MLR')"))
  }

  print(lavaan::summary(PMI.Model.fit, fit.measure = T, standardized = T, rsq = T))

  XX <- sapply(list(Model.config, PMI.Model.fit), lavaan::fitMeasures, c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr", "aic", "bic"))

  MODFIT <- matrix(1:18, nrow = 2, dimnames = list(c("Model.config", "PMI.Model.fit"), c("  chisq.scaled", "  df.scaled", "  pvalue.scaled", "  rmsea.scaled", "  cfi.scaled", "  tli.scaled", "  srmr", "  aic", "  bic")))
  MODFIT[1,1] <- format(round((XX[1,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[1,2] <- format(round((XX[2,1]), digits = 4), nsmall = 0, scientific = FALSE)
  MODFIT[1,3] <- format(round((XX[3,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[1,4] <- format(round((XX[4,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[1,5] <- format(round((XX[5,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[1,6] <- format(round((XX[6,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[1,7] <- format(round((XX[7,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[1,8] <- format(round((XX[8,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[1,9] <- format(round((XX[9,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,1] <- format(round((XX[1,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,2] <- format(round((XX[2,2]), digits = 4), nsmall = 0, scientific = FALSE)
  MODFIT[2,3] <- format(round((XX[3,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,4] <- format(round((XX[4,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,5] <- format(round((XX[5,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,6] <- format(round((XX[6,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,7] <- format(round((XX[7,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,8] <- format(round((XX[8,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,9] <- format(round((XX[9,2]), digits = 4), nsmall = 4, scientific = FALSE)

  FITDIF <- matrix(1:3, nrow = 1, dimnames = list(c("PMI.Model.fit - Model.config"), c("  cfi.scaled", "  rmsea.scaled", "    srmr")))
  FITDIF[1,1] <- format(round((XX[5,2] - XX[5,1]), digits = 4), nsmall = 4, scientific = FALSE)
  FITDIF[1,2] <- format(round((XX[4,2] - XX[4,1]), digits = 4), nsmall = 4, scientific = FALSE)
  FITDIF[1,3] <- format(round((XX[7,2] - XX[7,1]), digits = 4), nsmall = 4, scientific = FALSE)

  cat(rep("\n", 3))
  cat("####    MODEL FIT INDICES   ####", rep("\n", 2))
  print(MODFIT, quote=FALSE, right=TRUE)

  cat(rep("\n",2), "## Compare fit indices across models ##", rep("\n",2))
  print(lavaan::lavTestLRT(Model.config, PMI.Model.fit))

  cat(rep("\n", 2))
  cat("####  DIFFERENCES IN FIT INDICES  ####", rep("\n", 2))
  print(FITDIF, quote=FALSE, right=TRUE)
  cat(rep("\n", 3))

  ## == Print Factor Loadings == ##
  Final.Model.FL <- matrix(0, sum(no.items), no.group)
  rownames(Final.Model.FL) <- paste0(" Item ", names.ind, "    ")
  colnames(Final.Model.FL) <- paste0(" Group ", lavInspect(PMI.Model.fit, "group.label"), "   ")

  Final.par.est <- lavaan::parameterEstimates(PMI.Model.fit)

  ## Extract factor loadings ##
  for (G in 1: no.group) {
    ext <- c(which(Final.par.est[,"op"] == "=~" & Final.par.est[,"group"] == G))
    Final.Model.FL[, G] <- Final.par.est[ext, "est"]
  }
  cat("\n")
  cat("## === Factor Loadings in Final Model === ##")
  cat("\n")
  print(round(Final.Model.FL, digits=4))
  ## == End (Print Intercepts) == ##

  cat(rep("\n",2),"The recommended model PMI.Model.R is saved in the file 'PMI.txt'", "\n")

} ## end (Function CompareLoadings)

# ==================== Finish Function "CompareLoadings" ==================== #





# ==================== Create Function "CompareMeans" ==================== #
#' Scalar Invariance Test and Compare Latent Means
#'
#' Conduct scalar invariance test, identify a partial scalar invariance model, and compare latent means for data from independent groups.
#'
#' Reference: Measurement Equivalence/Invariance Test based on "Cheung, G. W. & Lau, R. S. (2012).  A direct comparison approach for testing measurement invariance.  Organizational Research Methods, 15, 167-198."
#'
#' Reference for correct.cfi: Widaman, K. F., & Thompson, J. S. (2003). On specifying the null model for incremental fit indices in structural equation modeling. Psychological Methods, 8(1), 16-37.
#'
#' @param model.PMI Partial metric invariance model (PMI.Model.R from CompareLoadings() or user-specified).
#' @param data.source A data frame containing the observed variables used in the model.
#' @param Groups Grouping variable for cross-group comparisons.
#' @param Cluster Cluster variable for nested data. The Monte Carlo simulation method should be used for nested data.
#' @param Bootstrap Number of bootstrap samples, must be between 500 and 10,000. If not specified, 1 million Monte Carlo simulated samples (default) will be used.
#' @param Type1 Overall Type I error rate (default is 0.05) for identifying non-invariant items in the List and Delete method.
#' @param Type1Adj Adjustment of Type I error rate for multiple tests for each construct. Default is Possible False Discovery Rate (PFDR = Type I error rate/No. of freely estimated factor loadings/intercepts across all groups), can also use (BON)ferroni adjustment (Type I error rate/No. of pairwise comparisons), or "NULL".
#' @param BMSC Best Model Selection Criterion used to select among models with the same number of estimated parameters. "ChiSquare" = minimum Chi-Square value, "CFI" = maximum CFI value, "RMSEA" = minimum RMSEA value, and "SRMR" = minimum SRMR value (default).
#' @param correct.cfi If TRUE (default), baseline model will be corrected for calculation of incremental fit indices based on Widaman & Thompson (2003).
#'
#' @return Partial scalar invariance model in PSI.txt file and results of latent means comparisons.
#' @export
#' @examples
#'
#' # Data file is "Example.A"
#'
#' # Specify the measurement model - Model.A
#' Model.A <- '
#'        WorkLifeConflict =~ R45a + R45b + R45c + R45d + R45e
#'        Engagement =~ R90a + R90b + R90c
#'        Wellbeing =~ R87a + R87b + R87c + R87d + R87e
#' '
#'
#' ## Not run:
#' ## ===== Full Measurement Invariance Test ===== ##
#' Full_MEI(Model.A, Example.A, Groups="Region")
#'
#' ## ===== Compare Loadings ===== ##
#' CompareLoadings(Model.A, Example.A, Groups="Region", Type1=0.05, Type1Adj="PFDR", BMSC="SRMR")
#' ## End (Not run)
#'
#' ## ===== Compare Intercepts and Latent Means ===== ##
#' CompareMeans(PMI.Model.R, Example.A, Groups="Region", Type1=0.05, Type1Adj="PFDR", BMSC="SRMR")
#'
CompareMeans <- function(model.PMI, data.source, Groups, Cluster="NULL", Bootstrap=0, Type1 = 0.05, Type1Adj = "PFDR", BMSC="SRMR", correct.cfi=TRUE) {

#  Bootstrap = 0 # Number of bootstrap samples
#  model.PMI = PMI.Model.R
#  data.source = Example.A
#  Groups = "Region"
#  Type1 = 0.05
#  Type1Adj = "PFDR"
#  Cluster="NULL"
#  BMSC="SRMR"
#  correct.cfi=TRUE

  Type1Adj <- toupper(Type1Adj)
  match.arg(Type1Adj, c("PFDR","BON", "NULL"))
  BMSC <- toupper(BMSC)
  match.arg(BMSC, c("CHISQUARE","CFI", "RMSEA", "SRMR"))
  if (is.logical(correct.cfi) == FALSE) stop("correct.cfi (corrected baseline model for calculating incremental fit indices) can only be TRUE or FALSE")

  arg1_char <- deparse(substitute(model.PMI))
  arg2_char <- deparse(substitute(data.source))
  arg3_char <- deparse(substitute(Groups))
  arg4_char <- deparse(substitute(Cluster))

  parsed <- lavParseModelString(model.PMI, as.data.frame = TRUE)
  parsed <- parsed[, c("lhs", "op", "rhs")]
  parsed <- parsed[parsed$op == "=~",]
  names.lv <- unique(parsed[,"lhs"])
  names.ind <- unique(parsed[, "rhs"])
  no.factor <- length(names.lv)
  no.items.g <- length(names.ind) # number of items per group
  no.items <- matrix(1:no.factor, nrow = 1)  # number of items per factor
  for (factor.no in 1:no.factor) { no.items[factor.no] <- sum(parsed[,"lhs"] == names.lv[factor.no] & parsed[,"op"] == "=~") } ## end loop factor.no

  names.item <- list()
  for (i in 1:no.factor) {
    names.item[[i]] <- parsed[which(parsed[,"lhs"] == names.lv[i] & parsed[,"op"] == "=~"), "rhs"]
  }

  ## Check for bootstrap sample number (Bootstrap) ##
  if (Bootstrap !=0) {
    b.no.integer <- Bootstrap == round(Bootstrap)
    if (b.no.integer == "FALSE") stop("Bootstrap sample number must be an integer")
    if (Bootstrap > 10000) stop("Bootstrap sample number greater than 10,000 is not recommended")
    if (Bootstrap < 500) stop("Bootstrap sample number smaller than 500 is not recommended")
    TYPE = "Bootstrap"
    b.no <- Bootstrap
  } else {
    TYPE = "MonteCarlo"
  } # end (Bootstrap != 0)

  count.tx <- 0


  ## ========== Run model.PMI ========== ##

  if (Cluster == "NULL") {
    PMI.Model.fit <- lavaan::sem(model.PMI,
      data.source,
      group = Groups,
      missing = 'fiml',
      auto.fix.first = FALSE,
      marker.int.zero = TRUE,
      meanstructure = T,
      information = 'observed',
      estimator = 'MLR')
  } else {
      PMI.Model.fit <- lavaan::sem(model.PMI,
      data.source,
      group = Groups,
      missing = 'fiml',
      auto.fix.first = FALSE,
      marker.int.zero = TRUE,
      meanstructure = T,
      information = 'observed',
      cluster = Cluster,
      estimator = 'MLR')
   }  # end Cluster

  ## Request summary outputs
  #$  print(lavaan::summary(PMI.Model.fit, fit.measure = T, standardized = T, rsq = T))

  ## ===== End (Run model.PMI) ===== ##


  par.est <- lavaan::coef(PMI.Model.fit)  # sample parameters
  no.group <- lavInspect(PMI.Model.fit, "ngroups")  # number of groups #
  group.names <<- lavaan::lavInspect(PMI.Model.fit, "group.label")  # group names

  # Find out the number of factors and the number of items per factor #
  parsed <- lavParseModelString(model.PMI, as.data.frame = TRUE)
  parsed <- parsed[, c("lhs", "op", "rhs")]
  parsed <- parsed[parsed$op == "=~",]

  temp <- lavaan::parameterEstimates(PMI.Model.fit)
  no.markers <- matrix(1:no.factor, nrow = 1)  # location of marker items
  for (factor.no in 1:no.factor) {
    no.markers[factor.no] <- which(temp[,"lhs"] == names.lv[factor.no] & temp[,"op"] == "=~" & temp[,"est"] == 1 & temp[,"group"] == 1)
    if (factor.no > 1) {no.markers[factor.no] <- no.markers[factor.no] - sum(no.items[1:(factor.no-1)])}
  } ## end loop factor.no

  temp <- lavaan::parameterEstimates(PMI.Model.fit, remove.nonfree=TRUE)
  simvcov <- lavInspect(PMI.Model.fit, what="vcov")

  ## -- Change labels to names -- ##
  names(par.est) <- paste0(temp[,"lhs"], temp[,"op"], temp[,"rhs"],".g", temp[,"group"])
  colnames(simvcov) <- paste0(temp[,"lhs"], temp[,"op"], temp[,"rhs"],".g", temp[,"group"])
  rownames(simvcov) <- paste0(temp[,"lhs"], temp[,"op"], temp[,"rhs"],".g", temp[,"group"])

  ## Extract factor loadings and indicator intercepts ##
  ext <- c(which(temp[,"op"] == "=~" & temp[,"group"] == 1), which(temp[,"op"] == "~1" & temp[,"group"] == 1 & temp[,"lhs"] %in% names.ind))
  for (i in 2: no.group) {
    ext <- c(ext, which(temp[,"op"] == "=~" & temp[,"group"] == i), which(temp[,"op"] == "~1" & temp[,"group"] == i & temp[,"lhs"] %in% names.ind))
  }

  ext <- which(temp[,"op"] == "=~" | temp[, "op"] == "~1")
  par.est <- par.est[ext]
  simvcov <- simvcov[ext,ext]
  no.par.g <- length(par.est)/no.group  # number of estimated parameters per group #
  no.lx.g <- sum(temp[, "op"] == "=~")/no.group # number of estimated LX per group #


  cat(rep("\n", 3), "## ======= SCALAR INVARIANCE ANALYSIS ======= ##", rep("\n", 2))  ## print heading

  # Calculate alpha for factor[factor.no] #
  ERate <- matrix(1:no.factor, nrow = 1)  # error rate for each factor

  cat("# -- Overall Type I Error Rates for Invariance Tests for Each Factor = ", Type1, " -- #", rep("\n", 2))

  if (Type1Adj == "PFDR") {  ## print adjusted error rate
    cat("# -- Adjusted Type I Error Rates for Invariance Tests (Adjusted by Possible False Discovery Rate) -- #", rep("\n", 2))
    for (factor.no in 1: no.factor) {
      no.adjust <- (no.items[factor.no]-1)*(no.group-1)
      ERate[factor.no] <- Type1 / no.adjust
      cat(paste0("Factor ", factor.no, ": Number of items = ", no.items[factor.no], "; number of free intercepts = ", no.adjust, "; adjusted error rate = ", sprintf("%.4f", ERate[factor.no]), "\n"))
    } # end (for factor.no)
  } else if (Type1Adj == "BON") {
    cat("# -- Adjusted Type I Error Rates for Invariance Tests (Bonforroni adjustment for number of pairwise comparisons) -- #", rep("\n", 2))
    for (factor.no in 1: no.factor) {
      no.adjust <- (no.items[factor.no]-1)*(no.group*(no.group-1)/2)
      ERate[factor.no] <- Type1 / no.adjust
      cat(paste0("Factor ", factor.no, ": Number of items = ", no.items[factor.no], "; number of pairwise comparisons = ", no.adjust, "; adjusted error rate = ", sprintf("%.4f", ERate[factor.no]), "\n"))
    } # end (for factor.no)
  } else if (Type1Adj == "NULL") {
    cat("# -- Error Rates for Invariance Tests (No adjustment is made) -- #", rep("\n", 2))
    cat(paste0("Type I Error Rate for Invariance Tests = ", Type1, "\n"))
    for (factor.no in 1: no.factor) {
      ERate[factor.no] <- Type1
    } # end (for factor.no)
  } # end (Type1Adj)
  cat(rep("\n",2))

  if (TYPE == "MonteCarlo") {

    mcmc <- MASS::mvrnorm(n=2000000, mu=par.est, Sigma=simvcov, tol = 1e-6)  # Run 2,000,000 simulations
    bootcoef <- mcmc
    bootno <- nrow(mcmc)  # No. of successful simulated samples
    cat(paste0("Number of Successful Simulated Samples = ", bootno, "\n"))

  } else { # TYPE=Bootstrap

    ## == Simplified bootstrapping model == ##
    PMI.boot <- lavaan::cfa(model.PMI, data = data.source, group=Groups,
                meanstructure = TRUE,
                auto.fix.first = FALSE,
                marker.int.zero = TRUE,
                ordered = FALSE,
                missing = "fiml", se = "none", test = "none", check.start = FALSE, check.post = FALSE, check.gradient = FALSE)

    ## == Bootstrapping == ##
    bootcoef <- lavaan::bootstrapLavaan(PMI.boot, R = b.no, ncpus = max(1L, parallel::detectCores() - 1L), parallel="snow")

    # Remove error and nonadmissible bootstrap samples #
    B.na <- attr(bootcoef,"nonadmissible")
    B.er <- attr(bootcoef,"error.idx")
    B.de <- c(B.na,B.er)
    if (length(B.de) != 0) {
      B.re <- bootcoef[-c(B.de),]
      bootcoef <- B.re
    }
    bootno <- nrow(bootcoef)  # number of successful bootstrap samples
    cat(paste0("Number of Successful Bootstrapped Samples = ", bootno, "\n"))
    bootcoef <- bootcoef[,ext]

  } ## end MonteCarlo or Bootstrap


  ## == Start the factor.no loop for CompareIntercepts == ##
  for (factor.no in 1: no.factor) {

    alpha <- ERate[factor.no]

    txY <<- matrix(" ",1, (no.group+2))
    txYY <<- matrix(" ",1, (no.group+2))
    EP <<- 1  # estimated parameter number

    no.k <- no.items[factor.no]*no.group  # number of factor loadings in FL.PMI
    FL.PMI <- matrix(1:no.k, nrow = no.group)  # factor loading matrix
    colnames(FL.PMI) <- paste0(' Item ', 1:no.items[factor.no])
    rownames(FL.PMI) <- lavaan::lavInspect(PMI.Model.fit, "group.label")

    TX.PMI <- matrix(1:no.k, nrow = no.group)  # Intercept matrix
    colnames(TX.PMI) <- paste0(' Item ', 1:no.items[factor.no])
    rownames(TX.PMI) <- lavaan::lavInspect(PMI.Model.fit, "group.label")


    for (i in 1:no.group) {
      for (j in 1:no.items[factor.no]) {
        if (j == no.markers[factor.no]) {
          FL.PMI[i, j] <- 1
          TX.PMI[i, j] <- 0
        } else {
          FL.PMI[i, j] <- format(round(par.est[paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[j]],".g",i)], digits=4), scientific = FALSE)
          TX.PMI[i, j] <- format(round(par.est[paste0(names.item[[factor.no]][[j]],"~1", ".g", i)], digits=4), scientific = FALSE)
        }  ## end if j
      }  ## end loop j
    }  ## end loop i

    class(FL.PMI) <- "numeric"
    class(TX.PMI) <- "numeric"

    for (Referent in 1:no.items[factor.no]) {  ## Loop referent item number for comparison

      if(length(unique(FL.PMI[,Referent])) != 1) {
        txZ <<- matrix(" ", 1, (no.group+2))
        txX <<- matrix(1:(no.group+2), 1, (no.group+2))
        txX[1, 1:no.group] <- 0
        uni.fl <- unique(FL.PMI[, Referent])
        Luni.fl <- length(uni.fl)
        for (uni in 1: Luni.fl) {
          flfl <- uni.fl[uni]
          for (temp.txZ in 1: no.group) {
            if (txZ[1, temp.txZ] == flfl) {
              txZ[1,temp.txZ] <- paste0("F", factor.no, "I", EP)
            }
          }
          EP <<- EP + 1
        }
        txX[1, (no.group+1)] <- Referent
        txZ[1, (no.group+1)] <- Referent

        Arg <- 1

        if (Referent > Arg) {
          txX[1, (no.group+2)] <- Arg
          txZ[1, (no.group+2)] <- Arg
        } else {
          txX[1, (no.group+2)] <- Arg + 1
          txZ[1, (no.group+2)] <- Arg + 1
        }  ## end if Referent
#$      print(txX, quote=FALSE)
#        txY <<- rbind(txY,txX)
#$      print(txZ, quote=FALSE)
#        txYY <<- rbind(txYY,txZ)

        next
       } # skip if unequal factor loadings


      ## ==  Print Factor Loadings FL of all groups  == ##
      FL <- matrix(1:no.k, nrow = no.group)  # factor loading matrix
      colnames(FL) <- paste0(' Item ', 1:no.items[factor.no])
      rownames(FL) <- lavaan::lavInspect(PMI.Model.fit, "group.label")


      if (Referent == no.markers[factor.no]) {  ## if referent is the marker item
        for (i in 1:no.group) {
          for (j in 1:no.items[factor.no]) {
            if (j == Referent) {
              FL[i, j] <- 1
            } else {
              FL[i, j] <- FL.PMI[i, j]
            }  ## end if j
          }  ## end for j
        }  ## end for i
      } else {  ## referent is not the marker item
        for (i in 1:no.group) {
          for (j in 1:no.items[factor.no]) {
            if (j == Referent) {
              FL[i, j] <- 1
            } else {
              FL[i, j] <- FL.PMI[i, j]/FL.PMI[i, Referent]
            }  ## end if j
          }  ## end for j
        }  ## end for i
      }  ## end if referent

#$      cat(rep("\n", 3), "Factor Loadings ", rep("\n", 2))
#$      print(round(FL[], digits=4))

      ## == End (Print Factor Loadings)  == ##

      ## ==  Print Intercepts of all groups  == ##
      no.k <- no.items[factor.no]*no.group  # number of intercepts in TX
      TX <- matrix(1:no.k, nrow = no.group)  # intercept matrix
      colnames(TX) <- paste0(' Item ', 1:no.items[factor.no])
      rownames(TX) <- lavaan::lavInspect(PMI.Model.fit, "group.label")

      if (Referent == no.markers[factor.no]) {  ## if referent is the marker item
        for (i in 1:no.group) {
          for (j in 1:no.items[factor.no]) {
            if (j == Referent) {
              TX[i, j] <- 0
            } else {
              TX[i, j] <- TX.PMI[i, j]
            }  ## end if j
          }  ## end for j
        }  ## end for i
      } else {  ## referent is not the marker item
        for (i in 1:no.group) {
          for (j in 1:no.items[factor.no]) {
            if (j == Referent) {
              TX[i, j] <- 0
            } else {
              TX[i, j] <- TX.PMI[i, j] - FL.PMI[i, j]/FL.PMI[i, Referent]*TX.PMI[i, Referent]
            }  ## end if j
          }  ## end for j
        }  ## end for i
      }  ## end if referent

#$      cat(rep("\n", 3), "Intercepts ", rep("\n", 2))
#$      print(round(TX[], digits=4))

      ## == End (Print Intercepts)  == ##


      for (nArg in 1: no.items[factor.no]) {  ## Argument item number for comparison
        no.dif <- factorial(no.group) / factorial(no.group - 2) / 2  # No. of pairwise comparisons for each item
        boot.dif.tx <- matrix(0, nrow(bootcoef), no.dif)  # Create bootstrap difference matrix
        samp.dif.tx <- matrix(0, no.group, no.group) # Create sample difference matrix
        if (nArg == Referent) {
          next
        } else {
          Arg <- nArg
        } ## end if nArg

        ## == Calculate bootstrap difference and sample estimate difference == ##
        comp = 0
        for (r in 1:(no.group-1)) {  ## r is referent group
          for (a in (r+1):no.group) {  ## a is argument group
            comp = comp + 1
            if (Referent == no.markers[factor.no]) {
              boot.dif.tx[,comp] <- bootcoef[, paste0(names.item[[factor.no]][[Arg]],"~1.g", r)] -
                                    bootcoef[, paste0(names.item[[factor.no]][[Arg]],"~1.g", a)]
              samp.dif.tx[r,a] <- par.est[paste0(names.item[[factor.no]][[Arg]], "~1.g",r)] -
                                  par.est[paste0(names.item[[factor.no]][[Arg]], "~1.g",a)]

            } else if (Referent != no.markers[factor.no]) {
              if (Arg == no.markers[factor.no]) {
                boot.dif.tx[, comp] <- (-1*bootcoef[, paste0(names.item[[factor.no]][[Referent]], "~1.g", r)]/
                                           bootcoef[, paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Referent]],".g", r)]) -
                                       (-1*bootcoef[, paste0(names.item[[factor.no]][[Referent]], "~1.g", a)]/
                                           bootcoef[, paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Referent]],".g", a)])
                samp.dif.tx[r,a] <- (-1*par.est[paste0(names.item[[factor.no]][[Referent]], "~1.g", r)]/
                                        par.est[paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Referent]],".g", r)]) -
                                    (-1*par.est[paste0(names.item[[factor.no]][[Referent]], "~1.g",a)]/
                                        par.est[paste0(names.lv[factor.no],"=~", names.item[[factor.no]][[Referent]],".g", a)])

              } else if (Arg != no.markers[factor.no]) {
                boot.dif.tx[,comp] <- (bootcoef[, paste0(names.item[[factor.no]][[Arg]], "~1.g",r)]-
                                       bootcoef[, paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Arg]],".g", r)]*
                                       bootcoef[, paste0(names.item[[factor.no]][[Referent]], "~1.g",r)]/
                                       bootcoef[, paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Referent]],".g", r)]) -
                                      (bootcoef[, paste0(names.item[[factor.no]][[Arg]], "~1.g",a)]-
                                       bootcoef[, paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Arg]],".g", a)]*
                                       bootcoef[, paste0(names.item[[factor.no]][[Referent]], "~1.g",a)]/
                                       bootcoef[, paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Referent]],".g", a)])
                samp.dif.tx[r,a] <- (par.est[paste0(names.item[[factor.no]][[Arg]],"~1.g",r)]-
                                     par.est[paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Arg]],".g", r)]*
                                     par.est[paste0(names.item[[factor.no]][[Referent]], "~1.g",r)]/
                                     par.est[paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Referent]],".g", r)]) -
                                    (par.est[paste0(names.item[[factor.no]][[Arg]],"~1.g",a)]-
                                     par.est[paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Arg]],".g", a)]*
                                     par.est[paste0(names.item[[factor.no]][[Referent]], "~1.g", a)]/
                                     par.est[paste0(names.lv[factor.no], "=~", names.item[[factor.no]][[Referent]],".g", a)])

              }  ## end if nArg
            }  ## end if Referent
          }  ## end loop a
        }  ## end loop r

#$ print(paste0("Referent = ", Referent, " Argument = ", nArg))
#$ print(samp.dif.tx)


        ## == Calculate Percentile Probability == ##

        F1.comp.pp <- matrix(0, no.group, no.group)  ## Percentile Probability
        pno_nipair <- 0  ## Set number of non-invariant pairs to zero

        comp = 0
        for (r in 1:(no.group-1)) {  ## r is the referent group
          for (a in (r+1):no.group) {  ## a is the argument
            comp = comp + 1
            if (quantile(boot.dif.tx[, comp], probs = 0.5, na.rm = TRUE) > 0) {
              F1.comp.pp[r,a] = 2*(sum(boot.dif.tx[, comp] < 0, na.rm=TRUE)/bootno)
            } else {
              F1.comp.pp[r,a] = 2*(sum(boot.dif.tx[, comp] > 0, na.rm=TRUE)/bootno)
            }  ## end if
          }  ## end loop a
        }  ## end loop r


        ## == Percentile Confidence Intervals of Pairwise Comparisons == ##

        PCI <- matrix(1:(no.dif*10), nrow = no.dif)
        colnames(PCI) <- c("Group 1","Group 2","0.5%","2.5%","5%","Estimate","95%","97.5%","99.5%","p-value")

        comp = 0
        for (r in 1:(no.group-1)) {
          for (a in (r+1):no.group) {
            comp = comp + 1
            PCI[comp, 1] <- r
            PCI[comp, 2] <- a
            PCI[comp, 3] <- quantile(boot.dif.tx[, comp],c(0.005), na.rm = TRUE)
            PCI[comp, 4] <- quantile(boot.dif.tx[, comp],c(0.025), na.rm = TRUE)
            PCI[comp, 5] <- quantile(boot.dif.tx[, comp],c(0.05), na.rm = TRUE)
            PCI[comp, 6] <- samp.dif.tx[r,a]
            PCI[comp, 7] <- quantile(boot.dif.tx[, comp],c(0.95), na.rm = TRUE)
            PCI[comp, 8] <- quantile(boot.dif.tx[, comp],c(0.975), na.rm = TRUE)
            PCI[comp, 9] <- quantile(boot.dif.tx[, comp],c(0.995), na.rm = TRUE)
            PCI[comp,10] <- F1.comp.pp[r,a]
            if (FL[r,Arg] != FL[a,Arg]) {
              if (pno_nipair == 0)	{
                pnipair <- c(r,a)
                pno_nipair <- pno_nipair + 1
              } else {
                pnipair <- c(pnipair,r,a)
                pno_nipair <- pno_nipair + 1
              }  ## end if pno_nipair
            } else {
              if (F1.comp.pp[r,a] < alpha) {
                if (pno_nipair == 0)	{
                  pnipair <- c(r,a)
                  pno_nipair <- pno_nipair + 1
                } else {
                  pnipair <- c(pnipair,r,a)
                  pno_nipair <- pno_nipair + 1
                }  ## end if pno_nipair
              }  ## end if F1.comp.pp
            }  ## end if FL[, Arg]
          }  ## end loop a
        }  ## end loop r


#$      cat("\n")
#$      if (Referent > Arg) {
#$        cat(paste0("Referent Item: ", Referent, "  Argument Item: ", Arg, "\n"))
#$      } else {
#$        cat(paste0("Referent Item: ", Referent, "  Argument Item: ", (Arg+1), "\n"))
#$      }  ## end if Referent
#$      cat("\n")
#$      print(round(PCI[], digits=4))

        count.tx <- count.tx + 1

        ## == Run List and Delete == ##
#$      cat(rep("\n", 2), "Sets of Groups with Invariant Intercepts", "\n")

        if (pno_nipair == no.dif) {

          txZ <<- matrix(" ", 1, (no.group+2))
          txX <<- matrix(1:(no.group+2), 1, (no.group+2))
          txX[1, 1:no.group] <- 0

          for (temp.txZ in 1: no.group) {
            txZ[1,temp.txZ] <- paste0("F", factor.no, "I", EP)
            EP <<- EP + 1
          }

          txX[1, (no.group+1)] <- Referent
          txZ[1, (no.group+1)] <- Referent
###          if (Referent > Arg) {
            txX[1, (no.group+2)] <- Arg
            txZ[1, (no.group+2)] <- Arg
###          } else {
###            txX[1, (no.group+2)] <- Arg + 1
###            txZ[1, (no.group+2)] <- Arg + 1
###          }  ## end if Referent
#$        print(txX, quote=FALSE)
          txY <<- rbind(txY,txX)
#$        print(txZ, quote=FALSE)
          txYY <<- rbind(txYY,txZ)

        } else if (pno_nipair > 0) {
#          txX <<- matrix(" ", 1, (no.group+2))
#          txZ <<- matrix(" ", 1, (no.group+2))

          listanddelete.tx(factor.no, no.group, no_nipair=pno_nipair, nipair=pnipair, Referent, Arg, txY, txYY, EP) ## Run list and delete tx

#$        print(txX, quote=FALSE)
#$        print(txZ, quote=FALSE)

        } else {  # (pno_nipair = 0)

          txX <<- matrix(1:(no.group+2), 1, (no.group+2))
          txZ <<- matrix(" ", 1, (no.group+2))
          txZ <- matrix(paste0("F", factor.no, "I", EP), 1, (no.group+2))
          EP <<- EP + 1
          txX[1, (no.group+1)] <- Referent
          txZ[1, (no.group+1)] <- Referent
###          if (Referent > Arg) {
            txX[1, (no.group+2)] <- Arg
            txZ[1, (no.group+2)] <- Arg
###          } else {
###            txX[1, (no.group+2)] <- Arg + 1
###            txZ[1, (no.group+2)] <- Arg + 1
###          }  ## end if Referent
#$        print(txX, quote=FALSE)
          txY <<- rbind(txY,txX)
#$        print(txZ, quote=FALSE)
          txYY <<- rbind(txYY,txZ)

        }  ## end if pno_nipair
      } ## end loop Arg
    } ## end loop Referent

    if (length(unique(txY[, no.group+1])) == 2 & no.items[factor.no] == 2) { # Only 2 groups with more than 2 items
      txX <<- matrix(1:(no.group+2), 1, (no.group+2))
      txZ <<- matrix(" ", 1, (no.group+2))
      txYa <<- txY[nrow(txY), no.group+1]
      txYb <<- txY[nrow(txY), no.group+2]
      txX[1, no.group + 1] <- txYb
      txX[1, no.group + 2] <- txYa
      txY <<- rbind(txY,txX)

      txZ[1, no.group + 1] <- txYb
      txZ[1, no.group + 2] <- txYa
      for (temp.txYY in 1: no.group) {
        txZ[1, temp.txYY] <- paste0("F", factor.no, "I", EP)
        EP <<- EP + 1
      }
      txYY <<- rbind(txYY,txZ)
    }

  txY <<- txY[-c(1),]
#$  cat("\n")
#$  print(txY, quote = FALSE)

  txYY <<- txYY[-c(1),]
#$  cat("\n")
#$  print(txYY, quote = FALSE)

  if (is.matrix(txY) == FALSE) { txY <<- matrix(txY, nrow=1) }
  if (is.matrix(txYY) == FALSE) { txYY <<- matrix(txYY, nrow=1) }

    NIcombine <- table(txY[, no.group+1], txY[, no.group+2])

    tempR <- unique(txYY[, no.group+1])

    for (aR in 1:length(tempR)) {
      R <- tempR[aR]
      for (r in 1:(no.items[factor.no])) {
        if (r == 1) {
          assign(paste0("Rset", r), c(length(which(txYY[,no.group+1] < R))+1:(NIcombine[R, r])))
        } else if (r == R) {
          assign(paste0("Rset", r), 0)
        } else {
          assign(paste0("Rset", r), c(length(which(txYY[,no.group+1] < R))+sum(NIcombine[R, 1:(r-1)])+1:(NIcombine[R, r])))
        }
        if (R == 1) { assign("Rset1", 0) }
      }  ## end loop r

      ## expand.grid(Rset1, Rset2, Rset3, Rset4, Rset5 ...)

      temp.row <- eval(parse(text = paste0("expand.grid(", paste0(c("Rset"), 1:no.items[factor.no], collapse=","),")")))
      if (R == 1) {
        Set.row <- temp.row
      } else {
        Set.row <- rbind(Set.row, temp.row)
      } ## end if R
    }  ## end loop aR

    Model.load <- array(" ", dim = c(no.items[factor.no], no.group, nrow(Set.row)))  # all possible invariance models

    for (Xset in 1: nrow(Set.row)) {
      a <- matrix(" ", no.items[factor.no], no.group)
      for (r in 1:no.items[factor.no]) {
        if (Set.row[Xset, r] == 0) {
          a[r, ] <- 0
        } else {
          a[r, ] <- txYY[Set.row[Xset, r], 1:no.group]
        } ## end if Set.row
###@@@        if (txYY[r, no.group+1] > 1) { break }
      }  ## end loop r
      Model.load[, , Xset] <- a

    }  ## end loop Xset

#$ print(Model.load)

##########################

  ## == Model.summary == ##
  # Column 1 - number of estimated intercepts
  # Column 2 - number of measurement invariance items
  # Column 3 - p-value of Model chi-square
  # Column 4 - CFI
  # Column 5 - RMSEA
  # Column 6 - SRMR

    Model.summary <- matrix(0, nrow(Set.row), 7)
    for (R in 1:nrow(Set.row)) {
      Model.summary[R, 1] <- length(unique(c(Model.load[,,R]))) - 1
      no.MI <- no.items[factor.no] - 1
      for (Ra in 1: no.items[factor.no]) {
        if (length(unique(c(Model.load[Ra,,R]))) > 1) {no.MI <- no.MI - 1}
      }  ## end loop Ra
      Model.summary[R, 2] <- no.MI
    }  ## end loop R


    ## == Recommended Models (Rec.Model) - Largest no. of MI items, then smallest number of estimated loadings == ##
    Rec.Model <- which(Model.summary[,2] == max(Model.summary[,2]))
    Rec.Model <- Rec.Model[which(Model.summary[Rec.Model,1] == min(Model.summary[Rec.Model,1]))]


    ## == Find recommended model if no.items > 2 == ##
    if (no.items[factor.no] > 2) {

      Rec.Model.load <- matrix(nrow = length(Rec.Model), ncol = length(lavaan::coef(PMI.Model.fit)))
      MODFIT <- matrix(0, nrow = length(Rec.Model), ncol = 9)
      colnames(MODFIT) <-
        c("  chisq.scaled", "  df.scaled", "  pvalue.scaled", "  rmsea.scaled", "  cfi.scaled", "  tli.scaled", "  srmr_bentler", "  aic", "  bic")
      Rec.Model.FL <- array(0, dim = c(no.group, no.items[factor.no], length(Rec.Model)))
      colnames(Rec.Model.FL) <- paste0(" Item ", 1:no.items[factor.no])
      rownames(Rec.Model.FL) <- paste0(" Group ", lavInspect(PMI.Model.fit, "group.label"), "   ")

      Rec.Model.TX <- array(0, dim = c(no.group, no.items[factor.no], length(Rec.Model)))
      colnames(Rec.Model.TX) <- paste0(" Item ", 1:no.items[factor.no])
      rownames(Rec.Model.TX) <- paste0(" Group ", lavInspect(PMI.Model.fit, "group.label"), "   ")

      PSI.XX <- array(" ", length(Rec.Model))
      for (Model.R in 1: length(Rec.Model)) {

        ## == Specify the factor loadings for partial scalar invariance models (PSI.X) == ##

        PSI.marker <- which(Model.load[,1,Rec.Model[Model.R]] == 0)
        if (length(PSI.marker) > 1) {PSI.marker <<- PSI.marker[1]}
        PSI.FL <- matrix(1:no.k, nrow = no.group)  # factor loading matrix
        PSI.load <- matrix(" ", no.items[factor.no], no.group)  # factor loading matrix

        for (FL.item.g in 1:no.group) {
          for (FL.item in 1:no.items[factor.no]) {
            PSI.FL[FL.item.g, FL.item] <- FL.PMI[FL.item.g,FL.item]/FL.PMI[FL.item.g,PSI.marker]
          }  ## end loop FL.item
        }  ## end loop FL.item.g

        FL.count <- 1
        for (FL.item in 1:no.items[factor.no]) {
          if (PSI.FL[FL.item.g, FL.item] == 1) {
            PSI.load[FL.item,] <- 1
          } else {
            temp.fl <- unique(PSI.FL[,FL.item])
            for (temp.x in 1:length(temp.fl)) {
              for (FL.item.g in 1:no.group) {
                if (PSI.FL[FL.item.g,FL.item] == temp.fl[temp.x]) {
                  PSI.load[FL.item,FL.item.g] <- paste0("F",factor.no,"L",FL.count)
                }  ## end if
              }  ## end loop FL.item.g
              FL.count <- FL.count + 1
            } ## end loop temp.x
          } ## end if PSI.FL
        }  ## end loop FL.item


        ## == Specify the partial metric invariance models (PMI.X) == ##

        PMI.X.1 <- paste0("PSI.X <- '", names.lv[factor.no], " =~ ", "c(", paste0(PSI.load[1, ], collapse=","), ")*", names.item[[factor.no]][[1]], " + ","\n")

        for (i in 2: (no.items[factor.no]-1)) {
          PMI.X.1 <- paste0(PMI.X.1, "c(",paste0(PSI.load[i, ], collapse=","),")*", names.item[[factor.no]][[i]], " + ","\n")
        } ## end loop i
        PMI.X.1 <- paste0(PMI.X.1, "c(",paste0(PSI.load[no.items[factor.no], ], collapse=","),")*", names.item[[factor.no]][[i+1]],"\n")

        ## == Specify the partial scalar invariance models (PSI.X) == ##

        final.intercept <- unique(as.vector(Model.load[,, Rec.Model[Model.R]]))
        for (j in 1:no.group) {
          for (i in 1:no.items[factor.no]) {
            if (Model.load[i,j,Rec.Model[Model.R]] != 0) { Model.load[i,j, Rec.Model[Model.R]] <- paste0("F", factor.no, "I", which(Model.load[i,j,Rec.Model[Model.R]] == final.intercept)) }
          } # end for i
        } # end for j

        PSI.XX[Model.R] <- paste0(PMI.X.1, "\n")
        for (i in 1: no.items[factor.no]) {
          PSI.XX[Model.R] <- paste0(PSI.XX[Model.R], "   ", names.item[[factor.no]][[i]], " ~ c(", paste0(Model.load[i, , Rec.Model[Model.R]], collapse=","), ")*1","\n")
        } ## end loop i
        PSI.XX[Model.R] <- paste0(PSI.XX[Model.R], "'","\n")
        eval(parse(text = PSI.XX[Model.R]))

        ## == Run PSI.X == ##
        if (Cluster == "NULL") {
          eval(parse(text=   "PSI.X.fit <- lavaan::sem(PSI.X,
                           data.source,
                           group = Groups,
                           missing = 'fiml',
                           auto.fix.first = FALSE,
                           marker.int.zero = TRUE,
                           meanstructure = T,
                           information = 'observed',
                           estimator = 'MLR')"))
        } else {
          eval(parse(text=   "PSI.X.fit <- lavaan::sem(PSI.X,
                           data.source,
                           group = Groups,
                           missing = 'fiml',
                           auto.fix.first = FALSE,
                           marker.int.zero = TRUE,
                           meanstructure = T,
                           information = 'observed',
                           cluster = Cluster,
                           estimator = 'MLR')"))
        }

        ## == Request summary outputs == ##
#$      eval(parse(text = "print(lavaan::summary(PSI.X.fit, fit.measure = T, standardized = T, rsq = T))"))


        ## == Save fit indices == ##
        if (correct.cfi == TRUE) {
          fitmeasures <- corrected.cfi(PSI.X.fit, data=data.source, Groups=Groups, cluster=Cluster)
          XX <- fitmeasures[c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr_bentler", "aic", "bic")]
        } else {
          XX <- lavaan::fitMeasures(PSI.X.fit,
                                    c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr_bentler", "aic", "bic"))
        } # end (if correct.cfi)


        MODFIT[Model.R,] <- format(round((XX), digits = 4), nsmall = 4, scientific = FALSE)

        ## == Save Intercepts == ##
        for (G in 1: no.group) {
          temp.tx <- eval(parse(text = paste0("lavaan::lavInspect(PSI.X.fit, what='est')$'", group.names[G], "'$nu")))
          Rec.Model.TX[G, , Model.R] <- temp.tx[1: no.items[factor.no]]
        }  ## end (for G)
      } ## end (for Model.R)

      ## == Print model fit indices == ##
      rownames(MODFIT) <- c(paste0("Model ", 1: length(Rec.Model)))
#$    cat(rep("\n", 3))
#$    cat("####    MODEL FIT INDICES   ####", rep("\n", 2))
#$    print(MODFIT, quote=FALSE, right=TRUE)

      ## == Print factor loadings == ##
#$    for (Model.R in 1: length(Rec.Model)) {
#$      cat(rep("\n", 2), "Factor Loadings of Model", Model.R, rep("\n", 2))
#$      print(round(Rec.Model.FL[,,Model.R], digits=4))
#$    }  ## end (for Model.R)


      ## == Save Recommended Model (Recommend.Model) -- First model with best model fit (R.Model) == ##
      class(MODFIT) <- "numeric"
      if (BMSC == "Chisquare") {
        Model.R <- which.min(MODFIT[,1])
      } else if (BMSC == "CFI") {
        Model.R <- which.max(MODFIT[,5])
      } else if (BMSC == "RMSEA") {
        Model.R <- which.min(MODFIT[,4])
      } else if (BMSC == "SRMR") {
        Model.R <- which.min(MODFIT[,7])
      }
      R.Model <- Model.R[1] # Select the first model if 2 or more models have the same fit
      if (factor.no == 1) {Recommend.Model <- matrix("PSI.Model.R <- '", 1)}  ## Reset Recommend.Model

      Recommend.Model <- rbind(Recommend.Model, PSI.XX[R.Model])

    }  ## end if no.items[factor.no] > 2


    if (no.items[factor.no] < 3) {   ## Begin Models with 2 items ##
        Rec.Model <- Rec.Model[1] # Select the first model if 2 or more models have the same fit
        PSI.marker <- which(Model.load[,1,Rec.Model] == 0)
        if (length(PSI.marker) > 1) {PSI.marker <<- PSI.marker[1]}
        PSI.FL <- matrix(1:no.k, nrow = no.group)  # factor loading matrix
        PSI.load <- matrix(" ", no.items[factor.no], no.group)  # factor loading matrix

        for (FL.item.g in 1:no.group) {
          for (FL.item in 1:no.items[factor.no]) {
            PSI.FL[FL.item.g, FL.item] <- FL.PMI[FL.item.g,FL.item]/FL.PMI[FL.item.g,PSI.marker]
          }  ## end loop FL.item
        }  ## end loop FL.item.g

        FL.count <- 1
        for (FL.item in 1:no.items[factor.no]) {
          if (PSI.FL[FL.item.g, FL.item] == 1) {
            PSI.load[FL.item,] <- 1
          } else {
            temp.fl <- unique(PSI.FL[,FL.item])
            for (temp.x in 1:length(temp.fl)) {
              for (FL.item.g in 1:no.group) {
                if (PSI.FL[FL.item.g, FL.item] == temp.fl[temp.x]) {
                  PSI.load[FL.item,FL.item.g] <- paste0("F",factor.no,"L",FL.count)
                }  ## end if
              }  ## end loop FL.item.g
              FL.count <- FL.count + 1
            } ## end loop temp.x
          } ## end if PSI.FL
        }  ## end loop FL.item

      R.Model <- Model.load[,,Rec.Model]  # Recommended Model

      Recommend.Model <- rbind(Recommend.Model, "\n")

      ## == Save Recommended Model (Recommend.Model) == ##
      if (factor.no == 1) {Recommend.Model <- matrix("PSI.Model.R <- '", 1)}
      Recommend.Model <- rbind(Recommend.Model, paste0("  ", names.lv[factor.no], " =~ ", "c(", paste0(PSI.load[1,], collapse=","), ")*", names.item[[factor.no]][[1]], " + ","\n"))
      Recommend.Model[nrow(Recommend.Model)] <- paste0(Recommend.Model[nrow(Recommend.Model)],"  c(",paste0(PSI.load[2, ], collapse=","),")*", names.item[[factor.no]][[2]],"\n")
      Recommend.Model[nrow(Recommend.Model)] <- paste0(Recommend.Model[nrow(Recommend.Model)],"   ", names.item[[factor.no]][[1]], " ~ c(",paste0(R.Model[1, ], collapse=","), ")*1","\n")
      Recommend.Model[nrow(Recommend.Model)] <- paste0(Recommend.Model[nrow(Recommend.Model)],"   ", names.item[[factor.no]][[2]], " ~ c(",paste0(R.Model[2, ], collapse=","), ")*1","\n")
    }  ## End (Models with 2 items)
  } ## End (loop factor.no)


  Recommend.Model[1] <- sub("''","'", Recommend.Model[1])
  for (i in 1:no.factor) {
    Recommend.Model[i+1] <-  sub("PSI.X <- '", " ", Recommend.Model[i+1])
    Recommend.Model[i+1] <- gsub("\n'", " ", Recommend.Model[i+1])
  } # end loop i
  Recommend.Model <- rbind(Recommend.Model, "  '")  ## last line of Recommended Model

  cat(rep("\n", 2), "## =====  Recommended Model  ===== ##", rep("\n", 2))
  for (i in 1: nrow(Recommend.Model)) {
    if (i == 1) {
      cat(Recommend.Model[i], "\n")
    } else {
      cat(Recommend.Model[i])
    }
  } ## end loop i

  cat(rep("\n", 2))
  cat("## ===== lavaan Outputs of PSI.Model.R ===== ##")
  cat(rep("\n", 2))

  ## == Start print recommended model to PSI.txt (Partial Metric Invariance Model) == ##
  sink('PSI.txt')

  cat("## =====  Recommended Model  ===== ##")
  cat("\n")
    for (i in 1: nrow(Recommend.Model)) {
      if (i == 1) {
        cat(Recommend.Model[i], "\n")
      } else {
        cat(Recommend.Model[i])
      }
    } ## end loop i


  cat(rep("\n",2), "## Run PSI.Model.R")
  cat(rep("\n",2), paste0("  PSI.Model.fit <- lavaan::sem(PSI.Model.R,"))
  cat("\n", paste0("    ", arg2_char, ","))
  cat("\n", paste0("    group = ", arg3_char,","))
  cat("\n", "    missing = 'fiml',")
  cat("\n", "    auto.fix.first = FALSE,")
  cat("\n", "    marker.int.zero = TRUE,")
  cat("\n", "    meanstructure = T,")
  cat("\n", "    information = 'observed',")
  if (Cluster != "NULL") { cat("\n", paste0("    cluster = ", arg4_char,",")) }
  cat("\n", "    estimator = 'MLR')")

  sink()  ## Stop writing to file

  eval(parse(text = Recommend.Model))
  ## == Run PSI.Model.R == ##
  if (Cluster == "NULL") {
    eval(parse(text=   "PSI.Model.fit <- lavaan::sem(PSI.Model.R,
                       data.source,
                       group = Groups,
                       missing = 'fiml',
                       auto.fix.first = FALSE,
                       marker.int.zero = TRUE,
                       meanstructure = T,
                       information = 'observed',
                       estimator = 'MLR')"))
  } else {
    eval(parse(text=   "PSI.Model.fit <- lavaan::sem(PSI.Model.R,
                       data.source,
                       group = Groups,
                       missing = 'fiml',
                       auto.fix.first = FALSE,
                       marker.int.zero = TRUE,
                       meanstructure = T,
                       information = 'observed',
                       cluster = Cluster,
                       estimator = 'MLR')"))
  }

  print(lavaan::summary(PSI.Model.fit, fit.measure = T, standardized = T, rsq = T))

  XX <- sapply(list(PMI.Model.fit, PSI.Model.fit), lavaan::fitMeasures, c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr", "aic", "bic"))
  if (correct.cfi == TRUE) {
    fitmeasures <- corrected.cfi(PSI.object=PSI.Model.fit, data=data.source, Groups=Groups, cluster=Cluster)
    XX[, 2] <- fitmeasures[c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr_bentler", "aic", "bic")]
  } # end (if correct.cfi)

  MODFIT <- matrix(1:18, nrow = 2, dimnames = list(c("PMI.Model.fit", "PSI.Model.fit"), c("  chisq.scaled", "  df.scaled", "  pvalue.scaled", "  rmsea.scaled", "  cfi.scaled", "  tli.scaled", "  srmr", "  aic", "  bic")))
  MODFIT[1,1] <- format(round((XX[1,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[1,2] <- format(round((XX[2,1]), digits = 4), nsmall = 0, scientific = FALSE)
  MODFIT[1,3] <- format(round((XX[3,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[1,4] <- format(round((XX[4,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[1,5] <- format(round((XX[5,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[1,6] <- format(round((XX[6,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[1,7] <- format(round((XX[7,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[1,8] <- format(round((XX[8,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[1,9] <- format(round((XX[9,1]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,1] <- format(round((XX[1,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,2] <- format(round((XX[2,2]), digits = 4), nsmall = 0, scientific = FALSE)
  MODFIT[2,3] <- format(round((XX[3,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,4] <- format(round((XX[4,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,5] <- format(round((XX[5,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,6] <- format(round((XX[6,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,7] <- format(round((XX[7,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,8] <- format(round((XX[8,2]), digits = 4), nsmall = 4, scientific = FALSE)
  MODFIT[2,9] <- format(round((XX[9,2]), digits = 4), nsmall = 4, scientific = FALSE)

  FITDIF <- matrix(1:3, nrow = 1, dimnames = list(c("PSI.Model.fit - PMI.Model.fit"), c("  cfi.scaled", "  rmsea.scaled", "    srmr")))
  FITDIF[1,1] <- format(round((XX[5,2] - XX[5,1]), digits = 4), nsmall = 4, scientific = FALSE)
  FITDIF[1,2] <- format(round((XX[4,2] - XX[4,1]), digits = 4), nsmall = 4, scientific = FALSE)
  FITDIF[1,3] <- format(round((XX[7,2] - XX[7,1]), digits = 4), nsmall = 4, scientific = FALSE)

  cat(rep("\n", 3))
  cat("####    MODEL FIT INDICES   ####", rep("\n", 2))
  print(MODFIT, quote=FALSE, right=TRUE)

  cat(rep("\n",2), "## Compare fit indices across models ##", rep("\n",2))
  print(lavaan::lavTestLRT(PMI.Model.fit, PSI.Model.fit))

  cat(rep("\n", 2))
  cat("####  DIFFERENCES IN FIT INDICES  ####", rep("\n", 2))
  print(FITDIF, quote=FALSE, right=TRUE)
  cat(rep("\n", 3))

  ## == Print Factor Loadings == ##
  Final.Model.FL <- matrix(0, sum(no.items), no.group)
  rownames(Final.Model.FL) <- paste0(" Item ", names.ind, "    ")
  colnames(Final.Model.FL) <- paste0(" Group ", lavInspect(PSI.Model.fit, "group.label"), "   ")

  Final.par.est <- lavaan::parameterEstimates(PSI.Model.fit)

  ## Extraxt factor loadings ##
  for (G in 1: no.group) {
    ext <- c(which(Final.par.est[,"op"] == "=~" & Final.par.est[,"group"] == G))
    Final.Model.FL[, G] <- Final.par.est[ext, "est"]
  }
  cat("\n")
  cat("## === Factor Loadings in Final Model === ##")
  cat("\n")
  print(round(Final.Model.FL, digits=4))
  ## == End (Print Intercepts) == ##



  ## == Print Intercepts == ##
  Final.Model.TX <- matrix(0, sum(no.items), no.group)
  rownames(Final.Model.TX) <- paste0(" Item ", names.ind, "    ")
  colnames(Final.Model.TX) <- paste0(" Group ", lavInspect(PSI.Model.fit, "group.label"), "   ")
  for (G in 1: no.group) { Final.Model.TX[,G] <- eval(parse(text = paste0("lavaan::lavInspect(PSI.Model.fit, what='est')$'", group.names[G], "'$nu"))) }
  cat("\n")
  cat("## === Intercepts in Final Model === ##")
  cat("\n")
  print(round(Final.Model.TX, digits=4))
  ## == End (Print Intercepts) == ##

  ## == Print Latent Means == ##
  Final.Model.LM <- matrix(0, no.factor, no.group)
  rownames(Final.Model.LM) <- paste0(" ", names.lv)
  colnames(Final.Model.LM) <- paste0("    Group ", lavInspect(PSI.Model.fit, "group.label"), "   ")
  for (G in 1: no.group) { Final.Model.LM[,G] <- eval(parse(text = paste0("lavaan::lavInspect(PSI.Model.fit, what='est')$'", group.names[G], "'$alpha"))) }
  cat("\n")
  cat("## === Latent Means in Final Model === ##")
  cat("\n")
  print(round(Final.Model.LM, digits=4))
  ## == End (Print Latent Means) == ##


  ## ========== Compare Latent Means ========== ##

  ## Extraxt latent means ##
  par.est <- lavaan::coef(PSI.Model.fit)  # sample parameters
  temp <- lavaan::parameterEstimates(PSI.Model.fit, remove.nonfree=TRUE)
  simvcov <- lavInspect(PSI.Model.fit, what="vcov")

  ## -- Change labels to names -- ##
  names(par.est) <- paste0(temp[,"lhs"], temp[,"op"], temp[,"rhs"],".g", temp[,"group"])
  colnames(simvcov) <- paste0(temp[,"lhs"], temp[,"op"], temp[,"rhs"],".g", temp[,"group"])
  rownames(simvcov) <- paste0(temp[,"lhs"], temp[,"op"], temp[,"rhs"],".g", temp[,"group"])

  parsed <- temp[temp$op == "=~",]
  names.lm <- unique(parsed[,"lhs"])

  ext <- c(which(temp$op %in% "~1" & temp$lhs %in% names.lm))
  par.est <- par.est[ext]
  simvcov <- simvcov[ext,ext]
  no.lm.g <- length(names.lv)  # number of estimated latent means per group #


  cat(rep("\n", 3), "## ======= COMPARISON OF LATENT MEANS ======= ##", rep("\n", 2))  ## print heading

  if (TYPE == "MonteCarlo") {

    mcmc <- MASS::mvrnorm(n=2000000, mu=par.est, Sigma=simvcov, tol = 1e-6)  # Run 2,000,000 simulations
    bootcoef <- mcmc
    bootno <- nrow(mcmc)  # No. of successful simulated samples
    cat(paste0("Number of Successful Simulated Samples = ", bootno, "\n"))

  } else { # TYPE=Bootstrap

    ## == Simplified bootstrapping model == ##
    PSI.boot <- lavaan::cfa(PSI.Model.R, data = data.source, group=Groups,
                  meanstructure = TRUE,
                  auto.fix.first = FALSE,
                  marker.int.zero = TRUE,
                  ordered = FALSE,
                  missing = "fiml", se = "none", test = "none", check.start = FALSE, check.post = FALSE, check.gradient = FALSE)

    ## == Bootstrapping == ##
    bootcoef <- lavaan::bootstrapLavaan(PSI.boot, R = b.no, ncpus = max(1L, parallel::detectCores() - 1L), parallel="snow")

    # Remove error and nonadmissible bootstrap samples #
    B.na <- attr(bootcoef,"nonadmissible")
    B.er <- attr(bootcoef,"error.idx")
    B.de <- c(B.na,B.er)
    if (length(B.de) != 0) {
      B.re <- bootcoef[-c(B.de),]
      bootcoef <- B.re
    }
    bootno <- nrow(bootcoef)  # number of successful bootstrap samples
    cat(paste0("Number of Successful Bootstrapped Samples = ", bootno, "\n"))
    bootcoef <- bootcoef[,ext]  ## Extract the bootstrapped intercepts

  } ## end MonteCarlo or Bootstrap


  ## == Start the factor.no loop for CompareMeans == ##
  for (factor.no in 1: no.factor) {

    boot.dif.lm <- matrix(0, nrow(bootcoef), no.dif)  # Create bootstrap difference matrix
    samp.dif.lm <- matrix(0, no.group, no.group) # Create sample difference matrix

    ## == Calculate bootstrap difference and sample estimate difference == ##
    comp = 0
    for (r in 1:(no.group-1)) {  ## r is referent group
      for (a in (r+1):no.group) {  ## a is argument group
        comp = comp + 1
        boot.dif.lm[,comp] <- bootcoef[, paste0(names.lv[factor.no], "~1.g", r)] - bootcoef[,paste0(names.lv[factor.no], "~1.g", a)]
        samp.dif.lm[r,a] <- par.est[paste0(names.lv[factor.no], "~1.g", a)] - par.est[paste0(names.lv[factor.no], "~1.g", r)]
        samp.dif.lm[a,r] <- par.est[paste0(names.lv[factor.no], "~1.g", r)] - par.est[paste0(names.lv[factor.no], "~1.g", a)]
      }  ## end loop a
    }  ## end loop r
    colnames(samp.dif.lm) <- paste0(" Group ", lavInspect(PMI.Model.fit, "group.label"), "   ")
    rownames(samp.dif.lm) <- paste0(" Group ", lavInspect(PMI.Model.fit, "group.label"), "   ")


    ## == Calculate Number of Invariant Intercepts == ##
    LM.inv.tau <- matrix(0, no.group, no.group)
    for (r in 1:(no.group-1)) {  ## r is the referent group
      for (a in (r+1):no.group) {  ## a is the argument
        if (factor.no == 1) {
          LM.inv.tau[r,a] <- sum(round(Final.Model.TX[1:no.items[factor.no],r],6) == round(Final.Model.TX[1:no.items[factor.no],a],6)) 
        } else {
          LM.inv.tau[r,a] <- sum(round(Final.Model.TX[(sum(no.items[1:(factor.no-1)])+1):sum(no.items[1:factor.no]),r], 6) == 
                                 round(Final.Model.TX[(sum(no.items[1:(factor.no-1)])+1):sum(no.items[1:factor.no]),a], 6)) 
        }
        LM.inv.tau[a,r] <- LM.inv.tau[r,a]
      } # end for a
    } # end for r 


    ## == Calculate Percentile Probability == ##
    LM.comp.pp <- matrix(0, no.group, no.group)  ## Percentile Probability
    colnames(LM.comp.pp) <- paste0(" Group ", lavInspect(PMI.Model.fit, "group.label"), "   ")
    rownames(LM.comp.pp) <- paste0(" Group ", lavInspect(PMI.Model.fit, "group.label"), "   ")

    comp = 0
    for (r in 1:(no.group-1)) {  ## r is the referent group
      for (a in (r+1):no.group) {  ## a is the argument
        comp = comp + 1
        if (quantile(boot.dif.lm[, comp], probs = 0.5, na.rm = TRUE) > 0) {
          LM.comp.pp[r,a] = 2*(sum(boot.dif.lm[, comp] < 0, na.rm=TRUE)/bootno)
          LM.comp.pp[a,r] = 2*(sum(boot.dif.lm[, comp] < 0, na.rm=TRUE)/bootno)
        } else {
          LM.comp.pp[r,a] = 2*(sum(boot.dif.lm[, comp] > 0, na.rm=TRUE)/bootno)
          LM.comp.pp[a,r] = 2*(sum(boot.dif.lm[, comp] > 0, na.rm=TRUE)/bootno)
        }  ## end if
      }  ## end loop a
    }  ## end loop r

    ## == Percentile Confidence Intervals of Pairwise Comparisons == ##

    PCI <- matrix(1:(no.dif*10), nrow = no.dif)
    colnames(PCI) <- c("Group 1","Group 2","0.5%","2.5%","5%","Estimate","95%","97.5%","99.5%","p-value")

    comp = 0
    for (r in 1:(no.group-1)) {
      for (a in (r+1):no.group) {
        comp = comp + 1
        PCI[comp, 1] <- r
        PCI[comp, 2] <- a
        PCI[comp, 3] <- round(quantile(boot.dif.lm[, comp],c(0.005), na.rm = TRUE), digits=4)
        PCI[comp, 4] <- round(quantile(boot.dif.lm[, comp],c(0.025), na.rm = TRUE), digits=4)
        PCI[comp, 5] <- round(quantile(boot.dif.lm[, comp],c(0.05), na.rm = TRUE), digits=4)
        PCI[comp, 6] <- round(samp.dif.lm[r,a], digits=4)
        PCI[comp, 7] <- round(quantile(boot.dif.lm[, comp],c(0.95), na.rm = TRUE), digits=4)
        PCI[comp, 8] <- round(quantile(boot.dif.lm[, comp],c(0.975), na.rm = TRUE), digits=4)
        PCI[comp, 9] <- round(quantile(boot.dif.lm[, comp],c(0.995), na.rm = TRUE), digits=4)
        PCI[comp,10] <- round(LM.comp.pp[r,a], digits=4)
      }  ## end loop a
    }  ## end loop r

    cat("\n")
    #$      print(PCI[],quote=F, nsmall=4, scientific=FALSE)


    # == Combine LM.comp.pp with LM.inv.tau == #
    LM.comp.pp <- matrix(sprintf("%.4f", LM.comp.pp), nrow = nrow(LM.comp.pp))
    LM.inv.tau <- matrix(sprintf("%.0f", LM.inv.tau), nrow = nrow(LM.inv.tau))
    LM.comp.pp <- matrix(paste0(LM.comp.pp, " (", LM.inv.tau, ")  "), nrow = nrow(LM.comp.pp), ncol = ncol(LM.comp.pp))
    colnames(LM.comp.pp) <- paste0(" Group ", lavInspect(PMI.Model.fit, "group.label"), "   ")
    rownames(LM.comp.pp) <- paste0(" Group ", lavInspect(PMI.Model.fit, "group.label"), "   ")
    diag(LM.comp.pp) <- ""

    samp.dif.lm <- matrix(sprintf("%.4f", samp.dif.lm), nrow = nrow(samp.dif.lm))
    diag(samp.dif.lm) <- ""
    colnames(samp.dif.lm) <- paste0(" Group ", lavInspect(PMI.Model.fit, "group.label"), "   ")
    rownames(samp.dif.lm) <- paste0(" Group ", lavInspect(PMI.Model.fit, "group.label"), "   ")

    cat("\n")
    cat("## ===== Latent Variable: ", names.lv[factor.no], " ===== ##")
    cat("\n")
    cat("== Factor Loadings ==")
    cat("\n")
    if (factor.no == 1) {
      print(format(round(Final.Model.FL[1:no.items[factor.no],],4)), quote=F)
    } else {
      print(format(round(Final.Model.FL[(sum(no.items[1:(factor.no-1)])+1):sum(no.items[1:factor.no]),],4)), quote=F)
    }
    cat("\n")
    cat("== Intercepts ==")
    cat("\n")
    if (factor.no == 1) {
      print(format(round(Final.Model.TX[1:no.items[factor.no],],4)), quote=F)
    } else {
      print(format(round(Final.Model.TX[(sum(no.items[1:(factor.no-1)])+1):sum(no.items[1:factor.no]),],4)), quote=F)
    }
    cat("\n")
    cat("== Latent Means ==")
    cat("\n")
    print(Final.Model.LM[factor.no,], quote=F, digits=5, nsmall=4, scientific=FALSE)
    cat("\n")
    cat("== Pairwise Difference in Latent Means ==")
    cat("\n")
    print(format(samp.dif.lm, justify="right"), quote=F)
    cat("\n")
    cat("== p-values for pairwise comparisons ==")
    cat("\n")
    print(format(LM.comp.pp, justify="right"), quote=F)
    cat("\n")
    cat("Note:", "\n")
    cat("Numbers in parentheses are numbers of items with invariant intercepts.", "\n")
    cat("There are ", no.items[factor.no], " items for this factor. Latent means should only be compared with at least ", ceiling((no.items[factor.no]+1)/2), " items with invariant intercepts.", "\n")
    cat("There are ", choose(no.group, 2), "pairwise comparisons.", "\n")
    cat("A critical p-value of", sprintf("%.4f", 0.05/choose(no.group, 2)), "is recommended for an overall Type I error rate of 0.05 with Bonferroni adjustment.", "\n")
    cat("A critical p-value of", sprintf("%.4f", 0.01/choose(no.group, 2)), "is recommended for an overall Type I error rate of 0.01 with Bonferroni adjustment.", "\n")
  } ## End (factor.no loop)

  cat(rep("\n",2),"The recommended model PSI.Model.R is saved in the file 'PSI.txt'", "\n")

} ## End (Function CompareMeans)

# ==================== Finish Function "CompareMeans" ==================== #





# ==================== Create Function "CompareParameters" ==================== #
#' Compare Defined Parameters Across Groups
#'
#' Compare defined parameters across groups, e.g., direct, indirect and total effects.
#'
#' Reference: Lau, R. S. & Cheung, G. W. (2012). Estimating and comparing specific mediation effects in complex latent variable models. Organizational Research Methods, 15, 3-16.
#'
#'
#' @param model.PMI Partial metric invariance model (PMI.Model.R from CompareLoadings() or user-specified).
#' @param model.PATH model with defined parameters.
#' @param data.source A data frame containing the observed variables used in the model.
#' @param Groups Grouping variable for cross-group comparisons.
#' @param Cluster Cluster variable for nested data. The Monte Carlo simulation method should be used for nested data.
#' @param Bootstrap Number of bootstrap samples, must be between 500 and 5,000. If not specified, 1 million Monte Carlo simulated samples (default) will be used.
#'

#' @return Estimates and confidence intervals for defined parameters in each group and comparisons of defined parameters across groups.
#' @keywords internal
#' @export
#' @examples
#'
#' ## Specify the measurement model - Model.A
#' Model.A <- '
#'        WorkLifeConflict =~ R45a + R45b + R45c + R45d + R45e
#'        Engagement =~ R90a + R90b + R90c
#'        Wellbeing =~ R87a + R87b + R87c + R87d + R87e
#' '
#'
#' ## Not run:
#' ## ===== Full Measurement Invariance Test ===== ##
#' Full_MEI(Model.A, Example.A, Groups = "Region")
#'
#' ## ===== Compare Loadings ===== ##
#' CompareLoadings(Model.A, Example.A, Groups = "Region", Type1 = 0.05, Type1Adj = "PFDR")
#' ## End (Not run)
#'
#' ## ===== Compare Parameters ===== ##
#' # -- Specify Path model - model.PATH (model.DP) [OrgSize and Tenure are control variables] -- #
#' model.DP <- '
#'   Wellbeing ~ Xb1*Engagement + Xc1*WorkLifeConflict + Xd1*OrgSize + Xe1*Tenure
#'   Engagement ~ Xa1*WorkLifeConflict
#'
#'  # Defined Parameters #
#'  IndirectP := Xa1*Xb1  # Indirect effect
#'  DirectP := Xc1  # Direct effect
#'  Total := Xa1*Xb1 + Xc1  # Total effect
#' '
#'
#' # -- Run function CompareParameters using  Monte Carlo simulation -- #
#' CompareParameters(PMI.Model.R, model.DP, Example.A, Groups = "Region")
#'
CompareParameters <- function(model.PMI, model.PATH, data.source, Groups, Cluster="NULL", Bootstrap=0) {

  #  model.PMI <- PMI.Model.R
  #  model.PATH <- model.DP
  #  Bootstrap = 2000 # Number of bootstrap samples
  #  data.source = Example.A
  #  Groups = "country"
  #  Cluster="NULL"

  model.SEM <- rbind(model.PMI, model.PATH)

  arg1_char <- deparse(substitute(model.SEM))
  arg2_char <- deparse(substitute(data.source))
  arg3_char <- deparse(substitute(Groups))

  ## Check for bootstrap sample number (Bootstrap) ##
  if (Bootstrap !=0) {
    b.no.integer <- Bootstrap == round(Bootstrap)
    if (b.no.integer == "FALSE") stop("Bootstrap sample number must be an integer")
    if (Bootstrap > 10000) stop("Bootstrap sample number greater than 10,000 is not recommended")
    if (Bootstrap < 500) stop("Bootstrap sample number smaller than 500 is not recommended")
    TYPE = "Bootstrap"
    b.no <- Bootstrap
  } else {
    TYPE = "MonteCarlo"
  }

  count.tx <- 0


  ## ========== Run model.SEM ========== ##
  if (Cluster == "NULL") {
    SEM.Model.fit <- suppressWarnings(lavaan::sem(model.SEM,
                                                  data.source,
                                                  group = Groups,
                                                  missing = 'fiml',
                                                  #     auto.fix.first = FALSE,
                                                  #     marker.int.zero = TRUE,
                                                  #     meanstructure = T,
                                                  information = 'observed',
                                                  estimator = 'MLR'))
  } else {
    SEM.Model.fit <- suppressWarnings(lavaan::sem(model.SEM,
                                                  data.source,
                                                  group = Groups,
                                                  missing = 'fiml',
                                                  #     auto.fix.first = FALSE,
                                                  #     marker.int.zero = TRUE,
                                                  #     meanstructure = T,
                                                  information = 'observed',
                                                  cluster=Cluster,
                                                  estimator = 'MLR'))
  }  # end Cluster

  ## Request summary outputs
  #$  print(lavaan::summary(SEM.Model.fit, fit.measure = T, standardized = T, rsq = T))

  ## ===== End (Run model.SEM) ===== ##

  no.group <- lavInspect(SEM.Model.fit, "ngroups")  # Number of groups #
  group.names <<- lavInspect(SEM.Model.fit, "group.label")


  ## ========= Create New Model with Defined Parameters for Each Group ======== ##

  est <- parameterEstimates(SEM.Model.fit)
  est <- est[est[, "op"] == "~",]  ## Subset of Data
  par.label <- unique(est[, "label"])
  #$ length(par.label)

  DP <- parameterEstimates(SEM.Model.fit)
  DP <- DP[DP[, "op"] == ":=",]  ## Subset of Data
  DP.name <- DP[, "lhs"]
  no.DP <- length(DP.name) # Number of defined parameters

  sc <- capture.output(cat(model.SEM)) # Change model to vector with row numbers
  LDP <- grep(":=", sc) # Rows with Defined Parameters
  ODP <- sc[LDP] # Extract the Defined Parameters

  Model.wo.DP <- sc # Model without defined parameters
  Model.wo.DP <- Model.wo.DP[-c(LDP)]
  #$ cat(Model.wo.DP)

  for (m in 1: length(par.label)) {
    A <- par.label[m]
    B <- ""
    for (n in 1: no.group) {
      if (n == 1) {
        B <- paste0(B, "c(GP", n, A,",")
      } else if (n == no.group) {
        B <- paste0(B, "GP", n, A,")")
      } else {
        B <- paste0(B, "GP", n, A, ",")
      }
    }
    Model.wo.DP <- gsub(A, B, Model.wo.DP)
  }
  for (i in 1:length(Model.wo.DP)) { Model.wo.DP[i] <- paste0(Model.wo.DP[i], " \n") }
  #$ cat(Model.wo.DP)

  sc <- capture.output(cat(model.SEM))
  LDP <- grep(":=", sc)
  ODP <- sc[LDP]
  NDP <- ODP
  for (n in 1: (no.group-1)) {NDP <- rbind(NDP, ODP)}

  for (n in 1: no.group) {
    for (m in 1: length(par.label)) {
      A <- par.label[m]
      B <- ""
      B <- paste0(B, "GP", n, A)
      NDP[n,] <- gsub(A, B, NDP[n,])
    } # end loop m
    for (p in 1:length(DP.name)) {
      A <- DP.name[p]
      B <- ""
      B <- paste0(B, "GP", n, A)
      NDP[n,] <- gsub(A, B, NDP[n,])
    } # end loop p
  } # end loop n

  k <- length(Model.wo.DP)
  for (i in 1:no.group) {
    for (j in 1:length(LDP)) {
      Model.wo.DP[k+1] <- paste0(NDP[i,j], "\n")
      k <- k + 1
    } # end loop j
  } # end loop i

  Model.DP.1 <- paste0("Model.DP <- ' \n")  # Model.DP is the revised Model
  for (i in 1:length(Model.wo.DP)) {Model.DP.1 <- paste0(Model.DP.1, Model.wo.DP[i])}
  Model.DP.1 <- paste0(Model.DP.1, "'", "\n")

  eval(parse(text = Model.DP.1))
  model.SEM <- Model.DP

  ## ================ ##

  ## ========== Run model.SEM ========== ##
  if (Cluster == "NULL") {
    SEM.Model.fit <- lavaan::sem(model.SEM,
                                 data.source,
                                 group = Groups,
                                 missing = 'fiml',
                                 auto.fix.first = FALSE,
                                 #     marker.int.zero = TRUE,
                                 #     meanstructure = T,
                                 information = 'observed',
                                 estimator = 'MLR')
  } else {
    SEM.Model.fit <- lavaan::sem(model.SEM,
                                 data.source,
                                 group = Groups,
                                 missing = 'fiml',
                                 auto.fix.first = FALSE,
                                 #     marker.int.zero = TRUE,
                                 #     meanstructure = T,
                                 information = 'observed',
                                 cluster=Cluster,
                                 estimator = 'MLR')
  }  # end Cluster


  ## Request summary outputs
  print(lavaan::summary(SEM.Model.fit, fit.measure = T, standardized = T, rsq = T))

  SEM.Model.fit <<- SEM.Model.fit # Save SEM.Model.fit to Global Environment

  ## ===== End (Run model.SEM) ===== ##

  ## ========== Compare Parameters ========== ##

  no.dif <- factorial(no.group) / factorial(no.group - 2) / 2  # No. of pairwise comparisons for defined parameter
  par.est <- lavaan::coef(SEM.Model.fit)  # sample parameters
  simvcov <- lavInspect(SEM.Model.fit, what="vcov")

  if (TYPE == "MonteCarlo") {

    mcmc <- MASS::mvrnorm(n=2000000, mu=par.est, Sigma=simvcov, tol = 1e-6)  # Run 2,000,000 simulations
    bootcoef <- mcmc
    bootno <- nrow(mcmc)  # No. of successful simulated samples
    cat(paste0("Number of Successful Simulated Samples = ", bootno, "\n"))

  } else { # TYPE=Bootstrap

    ## == Simplified bootstrapping model == ##
    SEM.boot <- lavaan::sem(model.SEM, data = data.source, group=Groups,
                            #               meanstructure = TRUE,
                            auto.fix.first = FALSE,
                            #               marker.int.zero = TRUE,
                            ordered = FALSE,
                            missing = "fiml", se = "none", test="none", check.start = FALSE, check.post = FALSE, check.gradient = FALSE)

    ## == Bootstrapping == ##
    bootcoef <- lavaan::bootstrapLavaan(SEM.boot, R = b.no, ncpus = max(1L, parallel::detectCores() - 1L), parallel="snow")
    bootno <- nrow(bootcoef)  # number of successful bootstrap samples

    # Remove error and nonadmissible bootstrap samples #
    B.na <- attr(bootcoef,"nonadmissible")
    B.er <- attr(bootcoef,"error.idx")
    B.de <- c(B.na,B.er)
    if (length(B.de) != 0) {
      B.re <- bootcoef[-c(B.de),]
      bootcoef <- B.re
    }
    bootno <- nrow(bootcoef)  # number of successful bootstrap samples

    ## Print Number of Successful Bootstrap Samples ##

    cat("\n")
    cat(paste0(" ===== Number of requested bootstrap samples, ", b.no, " ==== "))
    cat("\n")
    cat(paste0(" ===== Number of successful bootstrap samples, ", bootno, " ==== "))
    cat("\n")

  } ## end MonteCarlo or Bootstrap



  ## === Calculate Estimated (dfE) and Bootstrapped (dfB) Defined Parameters === ##

  sc <- capture.output(cat(model.SEM))
  LDP <- grep(":=", sc)
  ODP <- sc[LDP]
  NDP <- ODP  # New defined parameters
  BDP <- ODP  # For calculation of bootstrapped defined parameters
  EDP <- ODP  # For calculation of estimated defined parameters

  df1 <- "dfE <- data.frame("  # data.frame for estimated parameters
  df2 <- "dfB <- data.frame("  # data.frame for bootstrapped parameters
  for (n in 1: no.group) {
    for (m in 1: length(par.label)) {
      if (n == 1 & m == 1) {
        df1 <- paste0(df1, "par.est['GP", n, par.label[m],"']")
        df2 <- paste0(df2, "bootcoef[,'GP", n, par.label[m],"']")
      } else {
        df1 <- paste0(df1, ",", "par.est['GP", n, par.label[m],"']")
        df2 <- paste0(df2, ",", "bootcoef[,'GP", n, par.label[m],"']")
      }
    }
  }
  df1 <- paste0(df1,")")
  df2 <- paste0(df2,")")

  eval(parse(text = df1))
  eval(parse(text = df2))

  # Rename columns in dfE #
  k <- 1
  for (n in 1: no.group) {
    for (m in 1: length(par.label)) {
      names(dfE)[k] <- paste0("GP", n, par.label[m])
      colnames(dfB)[k] <- paste0("GP", n, par.label[m])
      k <- k + 1
    } # end loop m
  } # end loop n

  k <- 1
  for (n in 1: no.group) {
    for (p in 1:length(DP.name)) {
      for (m in 1: length(par.label)) {
        A <- paste0("GP", n, par.label[m])
        E <- paste0("dfE['",A,"']")
        B <- paste0("dfB['",A,"']")
        EDP[k] <- gsub(A, E, EDP[k])
        BDP[k] <- gsub(A, B, BDP[k])
      } # end loop m
      A <- paste0("GP",n,DP.name[p])
      E <- paste0("dfE['",A,"']")
      B <- paste0("dfB['",A,"']")
      EDP[k] <- gsub(A, E, EDP[k])
      BDP[k] <- gsub(A, B, BDP[k])
      EDP[k] <- gsub(":=", "<-", EDP[k])
      BDP[k] <- gsub(":=", "<-", BDP[k])
      k <- k + 1
    } # end loop p
  } # end loop n

  # Calculate Defined Parameters #
  for (n in 1: (length(DP.name)*no.group)) {
    eval(parse(text = EDP[n])) # Estimated
    eval(parse(text = BDP[n])) # Bootstrapped
  } # end loop n


  ## === Finish Calculation === ##


  Final.DP <- matrix(0, length(DP.name), no.group) # Create sample defined parameter matrix
  LL99.DP <- matrix(0, length(DP.name), no.group) # Create 99% CI Lower Limit
  LL95.DP <- matrix(0, length(DP.name), no.group) # Create 95% CI Lower Limit
  LL90.DP <- matrix(0, length(DP.name), no.group) # Create 90% CI Lower Limit
  UL99.DP <- matrix(0, length(DP.name), no.group) # Create 99% CI Upper Limit
  UL95.DP <- matrix(0, length(DP.name), no.group) # Create 95% CI Upper Limit
  UL90.DP <- matrix(0, length(DP.name), no.group) # Create 90% CI Upper Limit

  colnames(Final.DP) <- paste0(" Group ", lavInspect(SEM.Model.fit, "group.label"), "   ")
  rownames(Final.DP) <- paste0(" Defined Parameter ", DP.name, "   ")
  for (i in 1: no.group) {
    for (j in 1: length(DP.name)) {
      Final.DP[j,i] <- as.numeric(dfE[no.group*length(par.label)+(i-1)*length(DP.name)+j])
      LL99.DP[j,i] <- round(quantile(dfB[, no.group*length(par.label)+(i-1)*length(DP.name)+j],c(0.005), na.rm = TRUE), digits=4)
      LL95.DP[j,i] <- round(quantile(dfB[, no.group*length(par.label)+(i-1)*length(DP.name)+j],c(0.025), na.rm = TRUE), digits=4)
      LL90.DP[j,i] <- round(quantile(dfB[, no.group*length(par.label)+(i-1)*length(DP.name)+j],c(0.05), na.rm = TRUE), digits=4)
      UL90.DP[j,i] <- round(quantile(dfB[, no.group*length(par.label)+(i-1)*length(DP.name)+j],c(0.95), na.rm = TRUE), digits=4)
      UL95.DP[j,i] <- round(quantile(dfB[, no.group*length(par.label)+(i-1)*length(DP.name)+j],c(0.975), na.rm = TRUE), digits=4)
      UL99.DP[j,i] <- round(quantile(dfB[, no.group*length(par.label)+(i-1)*length(DP.name)+j],c(0.995), na.rm = TRUE), digits=4)
    } # end loop j
  } # end loop i

  ## == Start the DP.no loop for CompareParameters == ##
  for (DP.no in 1: no.DP) {

    boot.dif.DP <- matrix(0, bootno, no.dif)  # Create bootstrap difference matrix
    samp.dif.DP <- matrix(0, no.group, no.group) # Create sample difference matrix

    ## == Calculate bootstrap difference and sample estimate difference == ##
    comp = 0
    for (r in 1:(no.group-1)) {  ## r is referent group
      for (a in (r+1):no.group) {  ## a is argument group
        kr.DP <<- no.group*length(par.label)+(r-1)*length(DP.name)  ## location of DP before r group
        ka.DP <<- no.group*length(par.label)+(a-1)*length(DP.name)  ## location of DP before a group
        comp = comp + 1
        boot.dif.DP[,comp] <- dfB[,(kr.DP+DP.no)] - dfB[,(ka.DP+DP.no)]
        samp.dif.DP[r,a] <- as.numeric(dfE[kr.DP+DP.no] - dfE[ka.DP+DP.no])
      }  ## end loop a
    }  ## end loop r


    ## == Calculate Percentile Probability == ##

    DP.comp.pp <- matrix(0, no.group, no.group)  ## Percentile Probability
    colnames(DP.comp.pp) <- paste0(" Group ", lavInspect(SEM.Model.fit, "group.label"), "   ")
    rownames(DP.comp.pp) <- paste0(" Group ", lavInspect(SEM.Model.fit, "group.label"), "   ")

    comp = 0
    for (r in 1:(no.group-1)) {  ## r is the referent group
      for (a in (r+1):no.group) {  ## a is the argument
        comp = comp + 1
        if (quantile(boot.dif.DP[, comp], probs = 0.5, na.rm = TRUE) > 0) {
          DP.comp.pp[r,a] = 2*(sum(boot.dif.DP[, comp] < 0, na.rm = TRUE)/bootno)
          DP.comp.pp[a,r] = 2*(sum(boot.dif.DP[, comp] < 0, na.rm = TRUE)/bootno)
        } else {
          DP.comp.pp[r,a] = 2*(sum(boot.dif.DP[, comp] > 0, na.rm = TRUE)/bootno)
          DP.comp.pp[a,r] = 2*(sum(boot.dif.DP[, comp] > 0, na.rm = TRUE)/bootno)
        }  ## end if
      }  ## end loop a
    }  ## end loop r

    ## == Percentile Confidence Intervals of Pairwise Comparisons == ##

    no.k <- no.dif*10
    PCI <- matrix(1:no.k, nrow = no.dif)
    colnames(PCI) <- c("Group 1","Group 2","0.5%","2.5%","5%","Estimate","95%","97.5%","99.5%","p-value")

    comp = 0
    for (r in 1:(no.group-1)) {
      for (a in (r+1):no.group) {
        comp = comp + 1
        PCI[comp, 1] <- r
        PCI[comp, 2] <- a
        PCI[comp, 3] <- round(quantile(boot.dif.DP[, comp],c(0.005), na.rm = TRUE), digits=4)
        PCI[comp, 4] <- round(quantile(boot.dif.DP[, comp],c(0.025), na.rm = TRUE), digits=4)
        PCI[comp, 5] <- round(quantile(boot.dif.DP[, comp],c(0.05), na.rm = TRUE), digits=4)
        PCI[comp, 6] <- round(samp.dif.DP[r,a], digits=4)
        PCI[comp, 7] <- round(quantile(boot.dif.DP[, comp],c(0.95), na.rm = TRUE), digits=4)
        PCI[comp, 8] <- round(quantile(boot.dif.DP[, comp],c(0.975), na.rm = TRUE), digits=4)
        PCI[comp, 9] <- round(quantile(boot.dif.DP[, comp],c(0.995), na.rm = TRUE), digits=4)
        PCI[comp,10] <- round(DP.comp.pp[r,a], digits=4)
      }  ## end loop a
    }  ## end loop r

    #$    cat("\n")
    #$      print(PCI[],quote=F, nsmall=4, scientific=FALSE)

    ## Print Defined Parameters ##
    P.DP <- matrix(0, 7, no.group)
    colnames(P.DP) <- paste0(" Group ", lavInspect(SEM.Model.fit, "group.label"), "   ")
    rownames(P.DP) <- c("Estimated Defined Parameter", "Lower Limit 99% CI", "Lower Limit 95% CI", "Lower Limit 90% CI",
                        "Upper Limit 90% CI", "Upper Limit 95% CI", "Upper Limit 99% CI")
    P.DP[1,] <- Final.DP[DP.no,]
    P.DP[2,] <- LL99.DP[DP.no,]
    P.DP[3,] <- LL95.DP[DP.no,]
    P.DP[4,] <- LL90.DP[DP.no,]
    P.DP[5,] <- UL90.DP[DP.no,]
    P.DP[6,] <- UL95.DP[DP.no,]
    P.DP[7,] <- UL99.DP[DP.no,]

    cat("\n")
    cat("## ===== Defined Parameter: ", DP.name[DP.no], " ===== ##")
    cat("\n")
    print(P.DP, quote=F, row.names=TRUE, digits=4, nsmall=4, scientific=FALSE)
    cat("\n")
    cat("== p-values for pairwise comparisons ==")
    cat("\n")
    print(formatC(DP.comp.pp, digits=4, format="f"), quote=F)
    cat("\n")
    cat("Note:", "\n")
    cat("There are ", choose(no.group, 2), "pairwise comparisons.", "\n")
    cat("A critical p-value of", sprintf("%.4f", 0.05/choose(no.group, 2)), "is recommended for an overall Type I error rate of 0.05 with Bonferroni adjustment.", "\n")
    cat("A critical p-value of", sprintf("%.4f", 0.01/choose(no.group, 2)), "is recommended for an overall Type I error rate of 0.01 with Bonferroni adjustment.", "\n")
  } # end loop DP.no

} ## End (Function CompareParameters)

# ==================== Finish Function "CompareParameters" ==================== #





# ==================== Create Function "MLCompareLoadings" ==================== #
#' Metric Invariance Test Between Level 1 and Level 2
#'
#' Conduct configural invariance and full metric invariance test, and identify a partial metric invariance model for multilevel data.
#'
#' Reference: Measurement Equivalence/Invariance Test based on "Cheung, G. W. & Lau, R. S. (2012).  A direct comparison approach for testing measurement invariance.  Organizational Research Methods, 15, 167-198."
#'
#' Since all Level 1 (within-group) variables are group-mean centered, intercepts and latent means of all variables at Level 1 are set to zero. Hence, only configural and metric invariance across levels are tested.
#' Bootstrapping is not available for a multi-level model.
#'
#' Define the measurement model only once without specifying the level, and the function will compare the model across levels.
#'
#'
#' @param model User-specified measurement model.
#' @param data.source A data frame containing the observed variables used in the model.
#' @param Cluster Cluster variable for nested data. The Monte Carlo simulation method should be used for nested data.
#' @param Type1 Overall Type I error rate (default is 0.05) for identifying non-invariant items in the List and Delete method.
#' @param Type1Adj Adjustment of Type I error rate for multiple tests for each construct. Default is Possible False Discovery Rate (PFDR = Type I error rate/No. of freely estimated factor loadings/intercepts across all groups), can also use (BON)ferroni adjustment (Type I error rate/No. of pairwise comparisons), or "NULL".
#' @param BMSC Best Model Selection Criterion used to select among models with the same number of estimated parameters. "ChiSquare" = minimum Chi-Square value, "CFI" = maximum CFI value, "RMSEA" = minimum RMSEA value, and "SRMR" = minimum sum of SRMR-within and SRMR-between values (default).
#' @return Partial metric invariance model in PMI.txt file.
#' @export
#' @examples
#'
#' ## == Example D - Multilevel Confirmatory Factor Analysis == ##
#' # Data file is "Example.D"; cluster variable is "ID"
#'
#' ## Specify the measurement model - Model.D ##
#' Model.D <- 'OLBI =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8'
#'
#' ## ===== Compare Loadings ===== ##
#' MLCompareLoadings(Model.D, Example.D, Cluster="ID", Type1=0.05, Type1Adj="NULL", BMSC="SRMR")
#'
MLCompareLoadings <- function(model, data.source, Cluster="NULL", Type1 = 0.05, Type1Adj = "PFDR", BMSC="SRMR") {

  Type1Adj <- toupper(Type1Adj)
  match.arg(Type1Adj, c("PFDR","BON", "NULL"))
  BMSC <- toupper(BMSC)
  match.arg(BMSC, c("CHISQUARE","CFI", "RMSEA", "SRMR"))

  arg2_char <<- deparse(substitute(data.source))
  arg4_char <<- deparse(substitute(Cluster))

  no.group <- 2
  parsed <- lavParseModelString(model, as.data.frame = TRUE)
  parsed <- parsed[, c("lhs", "op", "rhs")]
  parsed <- parsed[parsed$op == "=~",]
  names.lv <- unique(parsed[,"lhs"])
  names.ind <- unique(parsed[, "rhs"])
  no.factor <- length(names.lv)
  no.items.g <- length(names.ind) # number of items per level
  no.items <- matrix(1:no.factor, nrow = 1)  # number of items per factor
  for (factor.no in 1:no.factor) { no.items[factor.no] <- sum(parsed[,"lhs"] == names.lv[factor.no] & parsed[,"op"] == "=~") } ## end loop factor.no

  names.item <- list()
  for (i in 1:no.factor) {
    names.item[[i]] <- parsed[which(parsed[,"lhs"] == names.lv[i] & parsed[,"op"] == "=~"), "rhs"]
  }

  ## ===== Create Multilevel Model Model.ML.X ===== ##
  ML.X <- paste0("Model.ML.X <- '", "\n", "level: 1", "\n")
  for (i in 1:factor.no) {
    ML.X <- paste0(ML.X, names.lv[i], "w =~ ")
    for (j in 1:no.items[i]) { ML.X <- paste0(ML.X, paste0(names.item[[i]][[j]], " + ")) } ## end loop j
    ML.X <- paste0(ML.X, names.item[[i]][[j]],"\n")
  } # end loop i
  ML.X <- paste0(ML.X, paste0("\n", "level: 2", "\n"))
  for (i in 1: factor.no) {
    ML.X <- paste0(ML.X, names.lv[i], "b =~ ")
    for (j in 1:no.items[i]) { ML.X <- paste0(ML.X, paste0(names.item[[i]][[j]], " + ")) } ## end loop j
    ML.X <- paste0(ML.X, names.item[[i]][[j]],"\n")
  } # end loop i
  ML.X <- paste0(ML.X, "'","\n")
  eval(parse(text = ML.X))


## ===== Run Configural Model ===== ##

  Model.config <- lavaan::sem(Model.ML.X,
      data.source,
      estimator = 'MLR',
      cluster=Cluster,
      verbose = FALSE, optim.method = 'em')

  ## Request summary outputs
  cat(rep("\n", 2), "## ===== Configural Invariance Model ===== ##", "\n")
  print(lavaan::summary(Model.config, fit.measure = T, standardized = T, rsq = T))
  cat("\n")

## ===== End (Run Configural Model) ===== ##


  par.est <- lavaan::coef(Model.config)  # sample parameters

  temp <- lavaan::parameterEstimates(Model.config)
  no.markers <- matrix(1:no.factor, nrow = 1)  # location of marker items
  for (factor.no in 1:no.factor) {
    no.markers[factor.no] <- which(temp[,"lhs"] == paste0(names.lv[factor.no],"w") & temp[,"op"] == "=~" & temp[,"est"] == 1)
    if (factor.no > 1) {no.markers[factor.no] <- no.markers[factor.no] - sum(no.items[1:(factor.no-1)])}
  } ## end loop factor.no

  temp <- lavaan::parameterEstimates(Model.config, remove.nonfree=TRUE)
  simvcov <- lavInspect(Model.config, what="vcov")
  par.est <- lavaan::coef(Model.config)  # sample parameters

  ## Extraxt factor loadings ##
  ext <- c(which(temp[,"op"] == "=~"))
  par.est <- par.est[ext]
  simvcov <- simvcov[ext,ext]
  no.par.g <- length(par.est)/2  # number of estimated LX per level #


  cat(rep("\n", 3), "## ======= METRIC INVARIANCE ANALYSIS ======= ##", rep("\n", 2))  ## print heading

# Calculate alpha for factor[factor.no] #
  ERate <- matrix(1:no.factor, nrow = 1)  # error rate for each factor

  cat("# -- Overall Type I Error Rates for Invariance Tests for Each Factor = ", Type1, " -- #", rep("\n", 2))

  if (Type1Adj == "PFDR") {  ## print adjusted error rate
    cat("# -- Adjusted Type I Error Rates for Invariance Tests (Adjusted by Possible False Discovery Rate) -- #", rep("\n", 2))
    for (factor.no in 1: no.factor) {
      no.adjust <- (no.items[factor.no]-1)
      ERate[factor.no] <- Type1 / no.adjust
      cat(paste0("Factor ", factor.no, ": Number of items = ", no.items[factor.no], "; number of free factor loadings = ", no.adjust, "; adjusted error rate = ", sprintf("%.4f", ERate[factor.no]), "\n"))
    } # end (for factor.no)
  } else if (Type1Adj == "BON") {
    cat("# -- Adjusted Type I Error Rates for Invariance Tests (Bonforroni adjustment for number of pairwise comparisons) -- #", rep("\n", 2))
    for (factor.no in 1: no.factor) {
      no.adjust <- (no.items[factor.no]-1)
      ERate[factor.no] <- Type1 / no.adjust
      cat(paste0("Factor ", factor.no, ": Number of items = ", no.items[factor.no], "; number of pairwise comparisons = ", no.adjust, "; adjusted error rate = ", sprintf("%.4f", ERate[factor.no]), "\n"))
    } # end (for factor.no)
  } else if (Type1Adj == "NULL") {
    cat("# -- Error Rates for Invariance Tests (No adjustment is made) -- #", rep("\n", 2))
    cat(paste0("Type I Error Rate for Invariance Tests = ", Type1, "\n"))
    for (factor.no in 1: no.factor) {
      ERate[factor.no] <- Type1
    } # end (for factor.no)
  } # end (Type1Adj)
  cat(rep("\n",2))

  ## Monte Carlo Simulation ##
  mcmc <- MASS::mvrnorm(n=2000000, mu=par.est, Sigma=simvcov, tol = 1e-6)  # Run 2,000,000 simulations
  bootcoef <- mcmc
  bootno <- nrow(mcmc)  # No. of successful simulated samples
  cat(paste0("Number of Successful Simulated Samples = ", bootno, "\n"))

  ## == Start the factor.no loop for CompareLoadings == ##
  for (factor.no in 1: no.factor) {

    alpha <- ERate[factor.no]

    flY <<- matrix(" ",1, 4)
    flYY <<- matrix(" ",1, 4)
    EP <<- 1  # estimated parameter number

    no.k <- no.items[factor.no]*2  # number of factor loadings in FL.PMI
    FL.PMI <- matrix(1:no.k, nrow = 2)  # factor loading matrix
    colnames(FL.PMI) <- paste0(' Item ', 1:no.items[factor.no])
    rownames(FL.PMI) <- c(paste0("Level ", 1:2))
    for (j in 1:no.items[factor.no]) {
      if (j == no.markers[factor.no]) {
        FL.PMI[1, j] <- 1
      } else {
        FL.PMI[1, j] <- format(round(par.est[paste0(names.lv[factor.no],"w=~", names.item[[factor.no]][[j]])], digits=4), scientific = FALSE)
      }  ## end if j
    }  ## end loop j
    for (j in 1:no.items[factor.no]) {
      if (j == no.markers[factor.no]) {
        FL.PMI[2, j] <- 1
      } else {
        FL.PMI[2, j] <- format(round(par.est[paste0(names.lv[factor.no],"b=~", names.item[[factor.no]][[j]],".l2")], digits=4), scientific = FALSE)
      }  ## end if j
    }  ## end loop j

    class(FL.PMI) <- "numeric"


    for (Referent in 1:no.items[factor.no]) {  ## Loop referent item number for comparison

      ## ==  Print Factor Loadings of all levels  == ##
      no.k <- no.items[factor.no]*2  # number of factor loadings in FL
      FL <- matrix(1:no.k, nrow = 2)  # factor loading matrix
      colnames(FL) <- paste0(' Item ', 1:no.items[factor.no])
      rownames(FL) <- lavaan::lavInspect(Model.config, "group.label")

      if (Referent == 1) {  ## if referent is the first item
        for (i in 1:2) {
          for (j in 1:no.items[factor.no]) {
            if (j == Referent) {
              FL[i, j] <- 1
            } else {
              FL[i, j] <- FL.PMI[i, j]
            }  ## end if j
          }  ## end for j
        }  ## end for i
      } else {  ## referent is not the first item
        for (i in 1:2) {
          for (j in 1:no.items[factor.no]) {
            if (j == Referent) {
              FL[i, j] <- 1
            } else {
              FL[i, j] <- FL.PMI[i, j]/FL.PMI[i, Referent]
            }  ## end if j
          }  ## end for j
        }  ## end for i
      }  ## end if referent

#$      cat(rep("\n", 3), "Factor Loadings ", rep("\n", 2))
#$      print(round(FL[], digits=4))

    ## == End (Print Factor Loadings)  == ##

      for (nArg in 1: no.items[factor.no]) {  ## Argument item number for comparison
        no.dif <- 1 # No. of pairwise comparisons for each item
        boot.dif.lx <- matrix(0, nrow(bootcoef), no.dif)  # Create bootstrap difference matrix
        samp.dif.lx <- matrix(0, 2, 2) # Create sample difference matrix
        if (nArg == Referent) {
          next
        } else {
          Arg <- nArg
        } ## end if nArg


      ## == Calculate bootstrap difference and sample estimate difference == ##
        comp = 0
        comp = comp + 1
        if (Referent == 1) {
          boot.dif.lx[,comp] <- bootcoef[, paste0(names.lv[factor.no],"w=~", names.item[[factor.no]][[Arg]])] -
                                bootcoef[, paste0(names.lv[factor.no],"b=~", names.item[[factor.no]][[Arg]],".l2")]
          samp.dif.lx[1,2] <- par.est[paste0(names.lv[factor.no],"w=~", names.item[[factor.no]][[Arg]])] -
                              par.est[paste0(names.lv[factor.no],"b=~", names.item[[factor.no]][[Arg]],".l2")]
        } else {
          if (Arg == no.markers[factor.no]) {
            boot.dif.lx[,comp] <- 1/bootcoef[, paste0(names.lv[factor.no],"w=~", names.item[[factor.no]][[Referent]])] -
                                  1/bootcoef[, paste0(names.lv[factor.no],"b=~", names.item[[factor.no]][[Referent]],".l2")]
            samp.dif.lx[1,2] <- 1/par.est[paste0(names.lv[factor.no],"w=~", names.item[[factor.no]][[Referent]])] -
                                1/par.est[paste0(names.lv[factor.no],"b=~", names.item[[factor.no]][[Referent]],".l2")]
          } else { # Arg != no.markers[factor.no]
            boot.dif.lx[,comp] <- bootcoef[,paste0(names.lv[factor.no],"w=~", names.item[[factor.no]][[Arg]])]/
                                  bootcoef[,paste0(names.lv[factor.no],"w=~", names.item[[factor.no]][[Referent]])] -
                                  bootcoef[,paste0(names.lv[factor.no],"b=~", names.item[[factor.no]][[Arg]],".l2")]/
                                  bootcoef[,paste0(names.lv[factor.no],"b=~", names.item[[factor.no]][[Referent]],".l2")]
            samp.dif.lx[1,2] <- par.est[paste0(names.lv[factor.no],"w=~", names.item[[factor.no]][[Arg]])]/
                                par.est[paste0(names.lv[factor.no],"w=~", names.item[[factor.no]][[Referent]])] -
                                par.est[paste0(names.lv[factor.no],"b=~", names.item[[factor.no]][[Arg]],".l2")]/
                                par.est[paste0(names.lv[factor.no],"b=~", names.item[[factor.no]][[Referent]],".l2")]
          }  ## end if Arg
        }  ## end if Referent

#$ print(samp.dif.lx)


      ## == Calculate Percentile Probability == ##

      F1.comp.pp <- matrix(0, 2, 2)  ## Percentile Probability
      pno_nipair <- 0  ## Set number of non-invariant pairs to zero

      comp = 0
      comp = comp + 1
      if (quantile(boot.dif.lx[, comp], probs = 0.5, na.rm = TRUE) > 0) {
        F1.comp.pp[1,2] = 2*(sum(boot.dif.lx[, comp] < 0, na.rm=TRUE)/bootno)
      } else {
        F1.comp.pp[1,2] = 2*(sum(boot.dif.lx[, comp] > 0, na.rm=TRUE)/bootno)
      }  ## end if


      ## == Percentile Confidence Intervals of Pairwise Comparisons == ##

      PCI <- matrix(1:(no.dif*10), nrow = no.dif)
      colnames(PCI) <- c("Level 1","Level 2","0.5%","2.5%","5%","Estimate","95%","97.5%","99.5%","p-value")

      comp = 0
      comp = comp + 1
      PCI[comp, 1] <- 1
      PCI[comp, 2] <- 2
      PCI[comp, 3] <- quantile(boot.dif.lx[, comp],c(0.005), na.rm = TRUE)
      PCI[comp, 4] <- quantile(boot.dif.lx[, comp],c(0.025), na.rm = TRUE)
      PCI[comp, 5] <- quantile(boot.dif.lx[, comp],c(0.05), na.rm = TRUE)
      PCI[comp, 6] <- samp.dif.lx[1,2]
      PCI[comp, 7] <- quantile(boot.dif.lx[, comp],c(0.95), na.rm = TRUE)
      PCI[comp, 8] <- quantile(boot.dif.lx[, comp],c(0.975), na.rm = TRUE)
      PCI[comp, 9] <- quantile(boot.dif.lx[, comp],c(0.995), na.rm = TRUE)
      PCI[comp,10] <- F1.comp.pp[1,2]
      if (F1.comp.pp[1,2] < alpha) {
        if (pno_nipair == 0)	{
          pnipair <- c(1,2)
          pno_nipair <- pno_nipair + 1
        } else {
          pnipair <- c(pnipair,1,2)
          pno_nipair <- pno_nipair + 1
        }  ## end if pno_nipair
      }  ## end if F1.comp.pp


      ## == Run List and Delete == ##
#$      cat(rep("\n", 2), "Sets of Groups with Invariant Loadings", "\n")
        if (pno_nipair == no.dif) {

          flZ <<- matrix(" ", 1, 4)
          flX <<- matrix(1:4, 1, 4)
          flX[1, 1:no.group] <- 0
          for (temp.flZ in 1: 2) {
            flZ[1,temp.flZ] <- paste0("F", factor.no, "L", EP)
            EP <<- EP + 1
          }
          flX[1, 3] <- Referent
          flZ[1, 3] <- Referent
#?          if (Referent > Arg) {
            flX[1, 4] <- Arg
            flZ[1, 4] <- Arg
#?          } else {
#?            flX[1, 4] <- Arg + 1
#?            flZ[1, 4] <- Arg + 1
#?          }  ## end if Referent
#$        print(flX, quote=FALSE)
          flY <<- rbind(flY,flX)
#$        print(flZ, quote=FALSE)
          flYY <<- rbind(flYY,flZ)

        } else if (pno_nipair > 0) {

          flX <<- matrix(" ", 1, 4)
          flZ <<- matrix(" ", 1, 4)

          listanddelete.lx(factor.no, no.group, no_nipair=pno_nipair, nipair=pnipair, Referent, Arg) ## Run list and delete lx

        } else {  # (pno_nipair = 0)

          flX <- matrix(1:4, 1, 4)
          flZ <<- matrix(" ", 1, 4)
          flZ <- matrix(paste0("F", factor.no, "L", EP), 1, (no.group+2))
          EP <<- EP + 1
          flX[1, 3] <- Referent
          flZ[1, 3] <- Referent
          flX[1, 4] <- Arg
          flZ[1, 4] <- Arg
          flY <<- rbind(flY,flX)
#$        print(flZ, quote=FALSE)
          flYY <<- rbind(flYY,flZ)

        }  ## end if pno_nipair
      } ## end loop Arg
    } ## end loop Referent

  if (length(unique(flY[, 3])) == 2 & no.items[factor.no] == 2) { # Only 2 groups with more than 2 items
    flX <<- matrix(1:4, 1, 4)
    flZ <<- matrix(" ", 1, 4)
    flYa <<- flY[nrow(flY), 3]
    flYb <<- flY[nrow(flY), 4]
    flX[1, 3] <- flYb
    flX[1, 4] <- flYa
    flY <<- rbind(flY,flX)

    flZ[1, 3] <- flYb
    flZ[1, 4] <- flYa
    for (temp.flYY in 1: 2) {
      flZ[1, temp.flYY] <- paste0("F", factor.no, "I", EP)
      EP <<- EP + 1
    }
    flYY <<- rbind(flYY,flZ)
  }

  flY <<- flY[-c(1),]
#$  cat("\n")
#$  print(flY, quote = FALSE)

  flYY <<- flYY[-c(1),]
#$  cat("\n")
#$  print(flYY, quote = FALSE)

  if (is.matrix(flYY) == FALSE) { flYY <<- matrix(flYY, nrow=1) }

  NIcombine <- table(flY[, 3], flY[, 4])
#  NI.level <- unique(flY[, 3])
#  NIcombine <- NIcombine[NI.level, NI.level]

  tempR <- unique(flYY[, 3])

    for (aR in 1:length(tempR)) {
      R <- tempR[aR]
      for (r in 1:no.items[factor.no]) {
        if (r == 1) {
          assign(paste0("Rset", r), c(length(which(flYY[,3] < R))+1:(NIcombine[R, r])))
        } else if (r == R) {
          assign(paste0("Rset", r), 0)
        } else {
          assign(paste0("Rset", r), c(length(which(flYY[,3] < R))+sum(NIcombine[R, 1:(r-1)])+1:(NIcombine[R, r])))
        }
        if (R == 1) { assign("Rset1", 0) }
      }  ## end loop r

      ## expand.grid(Rset1, Rset2, Rset3, Rset4, Rset5 ...)

      temp.row <- eval(parse(text = paste0("expand.grid(", paste0(c("Rset"), 1: no.items[factor.no], collapse=","),")")))
      if (aR == 1) {
        Set.row <- temp.row
      } else {
        Set.row <- rbind(Set.row, temp.row)
      } ## end if R
    }  ## end loop aR

    Model.load <- array(" ", dim = c(no.items[factor.no], 2, nrow(Set.row)))  # all possible invariance models

    for (Xset in 1: nrow(Set.row)) {
      a <- matrix(" ", no.items[factor.no], 2)
      for (r in 1:no.items[factor.no]) {
        if (Set.row[Xset, r] == 0) {
          a[r, ] <- 1
        } else {
          a[r, ] <- flYY[Set.row[Xset, r], 1:2]
        } ## end if Set.row
#?        if (flYY[r, 3] > 1) { break }
      }  ## end loop r
      Model.load[, , Xset] <- a
    }  ## end loop Xset

#$ print(Model.load)

  ## == Model.summary == ##
  # Column 1 - number of estimated loadings
  # Column 2 - number of measurement invariance items
  # Column 3 - p-value of Model chi-square
  # Column 4 - CFI
  # Column 5 - RMSEA
  # Column 6 - SRMR

    Model.summary <- matrix(0, nrow(Set.row), 7)
    for (R in 1:nrow(Set.row)) {
      Model.summary[R, 1] <- length(unique(c(Model.load[,,R]))) - 1
      no.MI <- no.items[factor.no] - 1
      for (Ra in 1: no.items[factor.no]) {
        if (length(unique(c(Model.load[Ra,,R]))) > 1) {no.MI <- no.MI - 1}
      }  ## end loop Ra
      Model.summary[R, 2] <- no.MI
    }  ## end loop R



    ## == Recommended Models (Rec.Model) - Largest no. of MI items, then smallest number of estimated loadings == ##
    Rec.Model <- which(Model.summary[,2] == max(Model.summary[,2]))
    Rec.Model <- Rec.Model[which(Model.summary[Rec.Model,1] == min(Model.summary[Rec.Model,1]))]

    ## == Find recommended model if no.items > 2 == ##
    if (no.items[factor.no] > 2) {

      Rec.Model.load <- matrix(nrow = length(Rec.Model), ncol = length(lavaan::coef(Model.config)))
      MODFIT <- matrix(0, nrow = length(Rec.Model), ncol = 10)
      colnames(MODFIT) <-
        c("  chisq", "  df", "  pvalue", "  rmsea", "  cfi", "  tli", "  srmr_within", "  srmr_between", "  aic", "  bic")
      Rec.Model.FL <- array(0, dim = c(2, no.items[factor.no], length(Rec.Model)))
      colnames(Rec.Model.FL) <- paste0(" Item ", 1:no.items[factor.no])
      rownames(Rec.Model.FL) <- c("Level: 1", "Level: 2")

      for (Model.R in 1: length(Rec.Model)) {

        ## == Specify the partial metric invariance models (PMI.X) == ##
        PMI.X.1 <- paste0("PMI.X <- '", "\n", "level: 1", "\n",
                          names.lv[factor.no], "w =~ ", paste0(Model.load[1, 1, Rec.Model[Model.R]]), "*", names.item[[factor.no]][[1]], " + ")

        for (i in 2: (no.items[factor.no]-1)) {
          PMI.X.1 <- paste0(PMI.X.1, Model.load[i, 1, Rec.Model[Model.R]],"*", names.item[[factor.no]][[i]], " + ")
        } ## end loop i
        PMI.X.1 <- paste0(PMI.X.1, Model.load[i+1, 1, Rec.Model[Model.R]],"*", names.item[[factor.no]][[i+1]])

        PMI.X.1 <- paste0(PMI.X.1, paste0("\n", "level: 2", "\n",
                          names.lv[factor.no], "b =~ ", paste0(Model.load[1, 2, Rec.Model[Model.R]]), "*", names.item[[factor.no]][[1]], " + "))

        for (i in 2: (no.items[factor.no]-1)) {
          PMI.X.1 <- paste0(PMI.X.1, Model.load[i, 2, Rec.Model[Model.R]],"*", names.item[[factor.no]][[i]], " + ")
        } ## end loop i
        PMI.X.1 <- paste0(PMI.X.1, Model.load[i+1, 2, Rec.Model[Model.R]],"*", names.item[[factor.no]][[i+1]], " '")

        eval(parse(text = PMI.X.1))

        ## == Run PMI.X == ##
          eval(parse(text=   "PMI.X.fit <- lavaan::sem(PMI.X,
                             data.source,
                             auto.fix.first = FALSE,
                             estimator = 'MLR',
 #                            marker.int.zero = TRUE,
 #                            meanstructure = T,
                            information = 'observed',
                             cluster = Cluster,
                             verbose = FALSE, optim.method = 'em')"))

        ## == Request summary outputs == ##
        #$ eval(parse(text = "print(lavaan::summary(PMI.X.fit, fit.measure = T, standardized = T, rsq = T))"))

        ## == Save fit indices == ##
        XX <- lavaan::fitMeasures(PMI.X.fit, fit.measures=
              c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr_within", "srmr_between", "aic", "bic"))
        MODFIT[Model.R,] <- format(round((XX), digits = 4), nsmall = 4, scientific = FALSE)

        ## == Save Factor Loadings == ##
        Rec.Model.FL[1, , Model.R] <- eval(parse(text = paste0("lavaan::lavInspect(PMI.X.fit, what='est')$within$lambda")))
        Rec.Model.FL[2, , Model.R] <- eval(parse(text = paste0("lavaan::lavInspect(PMI.X.fit, what='est')$",Cluster,"$lambda")))

      } ## end loop Model.R


      ## == Print model fit indices == ##
      rownames(MODFIT) <- c(paste0("Model ", 1: length(Rec.Model)))
#$    cat(rep("\n", 3))
#$    cat("####    MODEL FIT INDICES   ####", rep("\n", 2))
#$    print(MODFIT, quote=FALSE, right=TRUE)

    ## == Print factor loadings == ##
#$    for (Model.R in 1: length(Rec.Model)) {
#$      cat(rep("\n", 2), "Factor Loadings of Model", Model.R, rep("\n", 2))
#$      print(round(Rec.Model.FL[,,Model.R], digits=4))
#$    }  ## end loop Model.R

      ## == First model with best model fit (R.Model) == ##
      class(MODFIT) <- "numeric"
      if (BMSC == "Chisquare") {
        Model.R <- which.min(MODFIT[,1])
      } else if (BMSC == "CFI") {
        Model.R <- which.max(MODFIT[,5])
      } else if (BMSC == "RMSEA") {
        Model.R <- which.min(MODFIT[,4])
      } else if (BMSC == "SRMR") {
        if (min(MODFIT[,7]) < .08 & min(MODFIT[,8]) > .08) {
          Model.R <- which.min(MODFIT[,8])
        } else if (min(MODFIT[,7]) > .08 & min(MODFIT[,8]) < .08) {
          Model.R <- which.min(MODFIT[,7])
        } else {
          Model.R <- which((MODFIT[,7] + MODFIT[,8]) == min((MODFIT[,7]+MODFIT[,8])))
        }
      }
      Model.R <- Model.R[1]
      S.Model <- Model.load[,,Rec.Model[Model.R]]  # Recommended Model
      un.load <- sort(subset(unique(c(S.Model)), unique(c(S.Model))[] != 1))
      no.un.load <- length(un.load)
      R.Model <- S.Model
      for (j in 1:no.un.load) { R.Model[S.Model == un.load[j]] <- paste0("F", factor.no, "L", j) }  ## end loop j

    ## == Save Recommended Model (Recommend.Model) == ##
      if (factor.no == 1) {Recommend.Model <- matrix("PMI.Model.R <- '", 1)}
      Recommend.Model <- rbind(Recommend.Model, paste0("level: 1", "\n"))
      Recommend.Model <- rbind(Recommend.Model, paste0("   ", names.lv[factor.no], "w =~ ", paste0(Model.load[1, 1, Rec.Model[Model.R]]), "*", names.item[[factor.no]][[1]], " + "))

      for (i in 2: (no.items[factor.no]-1)) {
        Recommend.Model <- rbind(Recommend.Model, paste0(paste0(R.Model[i,1],"*", names.item[[factor.no]][[i]]," + ")))
      }  ## end loop i
      Recommend.Model <- rbind(Recommend.Model, paste0(paste0(R.Model[i+1,1],"*", names.item[[factor.no]][[i+1]], "\n")))

      Recommend.Model <- rbind(Recommend.Model, paste0("level: 2", "\n"))
      Recommend.Model <- rbind(Recommend.Model, paste0("   ", names.lv[factor.no], "b =~ ", paste0(Model.load[1, 2, Rec.Model[Model.R]]), "*", names.item[[factor.no]][[1]], " + "))

      for (i in 2: (no.items[factor.no]-1)) {
        Recommend.Model <- rbind(Recommend.Model, paste0(paste0(R.Model[i,2],"*", names.item[[factor.no]][[i]]," + ")))
      }  ## end loop i
      Recommend.Model <- rbind(Recommend.Model, paste0(paste0(R.Model[i+1,2],"*", names.item[[factor.no]][[i+1]], "\n")))

    }  ## end if no.items[factor.no] > 2

    if (no.items[factor.no] < 3) {   ## Begin Models with 2 items ##
      R.Model <- Model.load[,,1]  # Recommended Model
      un.load <- subset(unique(c(R.Model)), unique(c(R.Model))[] != 1)
      no.un.load <- length(un.load)
      for (j in 1:no.un.load) { R.Model[R.Model == un.load[j]] <- paste0("F", factor.no, "L", j) }  ## end loop j

      ## == Save Recommended Model (Recommend.Model) == ##
      if (factor.no == 1) {Recommend.Model <- matrix("PMI.Model.R <- '", 1)}
      Recommend.Model <- rbind(Recommend.Model, paste0("level: 1", "\n"))
      Recommend.Model <- rbind(Recommend.Model, paste0("   ", names.lv[factor.no], "w =~ ", paste0(Model.load[1, 1, Rec.Model[Model.R]]), "*", names.item[[factor.no]][[1]], " + "))
      Recommend.Model <- rbind(Recommend.Model, paste0(paste0(R.Model[2,1],"*", names.item[[factor.no]][[2]], "\n")))
      Recommend.Model <- rbind(Recommend.Model, paste0("level: 2", "\n"))
      Recommend.Model <- rbind(Recommend.Model, paste0("   ", names.lv[factor.no], "b =~ ", paste0(Model.load[1, 2, Rec.Model[Model.R]]), "*", names.item[[factor.no]][[1]], " + "))
      Recommend.Model <- rbind(Recommend.Model, paste0(paste0(R.Model[2,2],"*", names.item[[factor.no]][[2]], "\n")))

    }  ## End (Models with 2 items)

  } ## End (loop factor.no)

  Recommend.Model <- rbind(Recommend.Model, "  '")  ## last line of Recommended Model

  cat("\n", "## =====  Recommended Model  ===== ##", rep("\n", 2))
  for (i in 1: nrow(Recommend.Model)) { cat(Recommend.Model[i], "\n") } ## end loop i
  cat(rep("\n", 3))

  ## == Start print recommended model to PMI.txt (Partial Metric Invariance Model) == ##
  sink('PMI.txt')

  cat("## =====  Recommended Model  ===== ##")
  cat("\n")
  for (i in 1: nrow(Recommend.Model)) { cat(Recommend.Model[i], "\n") } ## end loop i

  cat(rep("\n",2), "## Run PMI.Model.R")
  cat(rep("\n",2), paste0("  PMI.Model.fit <- lavaan::sem(PMI.Model.R,"))
  cat("\n", paste0("    ", arg2_char, ","))
  cat("\n", "    auto.fix.first = FALSE,")
  cat("\n", "    information = 'observed',")
  cat("\n", paste0("    cluster = ", arg4_char,","))
  cat("\n", "    estimator = 'MLR',")
  cat("\n", "    verbose = FALSE, optim.method = 'em')")

  cat(rep("\n", 3), "## Request summary outputs")
  cat(rep("\n", 3), paste0("  print('## ===== Final Model ===== ##')", "\n"))
  cat(rep("\n", 2), paste0("  print(lavaan::summary(PMI.Model.fit, fit.measure = T, standardized = T, rsq = T))"), "\n")

  sink()  ## Stop writing to file
  source('PMI.txt') ## Run the script

  cat(rep("\n", 2), "## Compare fit indices across models ##", "\n")
  DMODFIT <- matrix(0, nrow = 3, ncol = 10)
  colnames(DMODFIT) <-
        c("  chisq", "  df", "  pvalue", "  rmsea", "  cfi", "  tli", "  srmr_within", "  srmr_between", "  aic", "  bic")
  rownames(DMODFIT) <- c("Model.config", "PMI.Model.fit", "PMI.Model.fit - Model.config")
  XX1 <- lavaan::fitMeasures(Model.config, fit.measures=
      c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr_within", "srmr_between", "aic", "bic"))
  DMODFIT[1,] <- XX1
  XX2 <- lavaan::fitMeasures(PMI.Model.fit, fit.measures=
      c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr_within", "srmr_between", "aic", "bic"))
  DMODFIT[2,] <- XX2

  DMODFIT[3,1] <- DMODFIT[2,1] -DMODFIT[1,1]
  DMODFIT[3,2] <- DMODFIT[2,2] -DMODFIT[1,2]
  DMODFIT[3,3] <- pchisq(q=DMODFIT[3,1], df=DMODFIT[3,2], lower.tail=FALSE)
  DMODFIT[3,4] <- DMODFIT[2,4] -DMODFIT[1,4]
  DMODFIT[3,5] <- DMODFIT[2,5] -DMODFIT[1,5]
  DMODFIT[3,6] <- DMODFIT[2,6] -DMODFIT[1,6]
  DMODFIT[3,7] <- DMODFIT[2,7] -DMODFIT[1,7]
  DMODFIT[3,8] <- DMODFIT[2,8] -DMODFIT[1,8]
  DMODFIT[3,9] <- DMODFIT[2,9] -DMODFIT[1,9]
  DMODFIT[3,10] <- DMODFIT[2,10] -DMODFIT[1,10]
  print(round(DMODFIT, digits=4))


  ## == Print Factor Loadings == ##
  Final.Model.FL <- matrix(0, sum(no.items), 2)
  rownames(Final.Model.FL) <- paste0(" Item ", names.ind, "    ")
  colnames(Final.Model.FL) <- c("level: 1", "level: 2")

  Final.par.est <- lavaan::parameterEstimates(PMI.Model.fit)

  ## Extraxt factor loadings ##
  for (G in 1: 2) {
    ext <- c(which(Final.par.est[,"op"] == "=~" & Final.par.est[,"level"] == G))
    Final.Model.FL[, G] <- Final.par.est[ext, "est"]
  }
  cat("\n")
  cat("## === Factor Loadings in Final Model === ##")
  cat("\n")
  print(round(Final.Model.FL, digits=4))
  ## == End (Print Intercepts) == ##

  cat(rep("\n",2),"The recommended model PMI.Model.R is saved in the file 'PMI.txt'", "\n")

} ## End (Function MLCompareLoadings)

# ==================== Finish Function "MLCompareLoadings" ==================== #





# ==================== Create Function "LGCompareLoadings" ==================== #
#' Metric Invariance for Longitudinal Panel Data
#'
#' Conduct metric invariance test and identify a partial metric invariance model for longitudinal panel data.
#'
#' Requires defining the measurement model only once. All indicators should have the suffix _T1, _T2 and _T3 (e.g., x1_T1, x1_T2, x1_T3) in the data file to indicate when the indicators were measured.
#'
#' Residuals of indicators are covaried across time automatically.
#'
#' @param model User-specified measurement model.
#' @param data.source A data frame containing the observed variables used in the model.
#' @param no.waves Number of waves for comparisons.
#' @param Cluster Cluster variable for nested data. The Monte Carlo simulation method should be used for nested data.
#' @param Bootstrap Number of bootstrap samples, must be between 500 and 10,000. If not specified, 1 million Monte Carlo simulated samples (default) will be used.
#' @param Type1 Overall Type I error rate (default is 0.05) for identifying non-invariant items in the List and Delete method.
#' @param Type1Adj Adjustment of Type I error rate for multiple tests for each construct. Default is Possible False Discovery Rate (PFDR = Type I error rate/No. of freely estimated factor loadings/intercepts across all waves), can also use (BON)ferroni adjustment (Type I error rate/No. of pairwise comparisons), or "NULL".
#' @param BMSC Best Model Selection Criterion used to select among models with the same number of estimated parameters. "ChiSquare" = minimum Chi-Square value, "CFI" = maximum CFI value, "RMSEA" = minimum RMSEA value, and "SRMR" = minimum SRMR value (default).
#' @return Partial metric invariance model in the PMI.txt file.
#' @export
#' @examples
#'
#' ## == Example B - Panel Data in Longitudinal Studies == ##
#'
#' # Data file is "Example.B"
#'
#' ## Specify the measurement model - Model.B ##
#' Model.B <- '
#'   SWLS =~ x1 + x2 + x3 + x4 + x5
#' '
#'
#' ## ===== Compare Factor Loadings ===== ##
#' LGCompareLoadings(Model.B, Example.B, no.waves=3, Type1=0.05, Type1Adj="PFDR", BMSC="SRMR")
#'
#'
LGCompareLoadings <- function(model, data.source, Cluster="NULL", no.waves=3, Bootstrap=0, Type1 = 0.05, Type1Adj = "PFDR", BMSC="SRMR") {

#  model <<- model
#  Bootstrap = 0 # Number of bootstrap samples
#  model = Model.B
#  data.source = Data # Example.B
#  Type1 = 0.05
#  Type1Adj = "PFDR"
#  Cluster="NULL"
#  no.waves <- 3
#  BMSC="SRMR"

  model <- paste0("\n", "  ", model, "\n")
  Type1Adj <- toupper(Type1Adj)
  match.arg(Type1Adj, c("PFDR","BON", "NULL"))
  BMSC <- toupper(BMSC)
  match.arg(BMSC, c("CHISQUARE","CFI", "RMSEA", "SRMR"))

  arg2_char <<- deparse(substitute(data.source))
  arg4_char <<- deparse(substitute(Cluster))

  no.group <- no.waves
  parsed <- lavParseModelString(model, as.data.frame = TRUE)
  parsed <- parsed[, c("lhs", "op", "rhs")]
  parsed <- parsed[parsed$op == "=~",]
  names.lv <- unique(parsed[,"lhs"])
  names.ind <- unique(parsed[, "rhs"])
  no.factor <- length(names.lv)
  no.items.g <- length(names.ind) # number of items per group
  no.items <- matrix(1:no.factor, nrow = 1)  # number of items per factor
  for (factor.no in 1:no.factor) { no.items[factor.no] <- sum(parsed[,"lhs"] == names.lv[factor.no] & parsed[,"op"] == "=~") } ## end loop factor.no

  names.item <- list()
  for (i in 1:no.factor) {
    names.item[[i]] <- parsed[which(parsed[,"lhs"] == names.lv[i] & parsed[,"op"] == "=~"), "rhs"]
  }

  ## Check for bootstrap sample number (Bootstrap) ##
  if (Bootstrap !=0) {
    b.no.integer <- Bootstrap == round(Bootstrap)
    if (b.no.integer == "FALSE") stop("Bootstrap sample number must be an integer")
    if (Bootstrap > 10000) stop("Bootstrap sample number greater than 10,000 is not recommended")
    if (Bootstrap < 500) stop("Bootstrap sample number smaller than 500 is not recommended")
    TYPE = "Bootstrap"
    b.no <- Bootstrap
  } else {
    TYPE = "MonteCarlo"
  } # end (Bootstrap != 0)

  for (i in 1: no.factor) { model <- eval(parse(text = (paste0("sub('", names.lv[i],"', '", names.lv[i], "_T1', model, fixed=TRUE)")))) }
  for (i in 1: no.items.g) { model <- eval(parse(text = (paste0("sub('", names.ind[i],"', '", names.ind[i], "_T1', model, fixed=TRUE)")))) }


  ## -- Create Model.Long with all waves and covary residuals -- ##
  Model.Long <- model
  for (i in 2: no.waves) { Model.Long <- paste0(Model.Long, gsub("T1", paste0("T", i), model)) }
  Model.Long <- paste0(Model.Long, "\n")

  ## == Create Residual Covariance == ##
  for (i in 1:no.items.g) {
    for (j in 1:(no.waves - 1)) {
      for (k in (j+1):no.waves) {
        Model.Long <- paste0(Model.Long, "  ", names.ind[i], "_T", j, " ~~ ", names.ind[i], "_T", k, "\n")
      }
    }
  }
  ## -- Finish creating Model.Long -- ##


  ## ===== Run Configural Model ===== ##

  if (Cluster == "NULL") {
    Model.Long.config <- lavaan::sem(Model.Long,
                                data.source,
                                missing = 'fiml',
                                marker.int.zero = TRUE,
                                meanstructure = T,
                                estimator = 'MLR')
  } else {
    Model.Long.config <- lavaan::sem(Model.Long, cluster=Cluster,
                                data.source,
                                missing = 'fiml',
                                marker.int.zero = TRUE,
                                meanstructure = T,
                                estimator = 'MLR')
  } # end if cluster

  ## Request summary outputs
  cat(rep("\n", 2), "## ===== Configural Invariance Model ===== ##", "\n")
  print(lavaan::summary(Model.Long.config, fit.measure = T, standardized = T, rsq = T))
  #$  parameterEstimates(Model.Long.config)
  cat("\n")

  Model.Long.config <<- Model.Long.config

  ## ===== End (Run Configural Model) ===== ##


  par.est <- lavaan::coef(Model.Long.config)  # sample parameters

  temp <- lavaan::parameterEstimates(Model.Long.config)
  no.markers <- matrix(1:no.factor, nrow = 1)  # location of marker items
  for (factor.no in 1:no.factor) {
    no.markers[factor.no] <- which(temp[,"lhs"] == paste0(names.lv[factor.no],"_T1") & temp[,"op"] == "=~" & temp[,"est"] == 1)
    if (factor.no > 1) {no.markers[factor.no] <- no.markers[factor.no] - sum(no.items[1:(factor.no-1)])}
  } ## end loop factor.no

  temp <- lavaan::parameterEstimates(Model.Long.config, remove.nonfree=TRUE)
  simvcov <- lavInspect(Model.Long.config, what="vcov")
  par.est <- lavaan::coef(Model.Long.config)  # sample parameters

  ## Extraxt factor loadings ##
  ext <- c(which(temp[,"op"] == "=~"))
  par.est <- par.est[ext]
  simvcov <- simvcov[ext,ext]
  no.par.g <- length(par.est)/2  # number of estimated LX per group #


  cat(rep("\n", 3), "## ======= METRIC INVARIANCE ANALYSIS ======= ##", rep("\n", 2))  ## print heading

  # Calculate alpha for factor[factor.no] #
  ERate <- matrix(1:no.factor, nrow = 1)  # error rate for each factor

  cat("# -- Overall Type I Error Rates for Invariance Tests for Each Factor = ", Type1, " -- #", rep("\n", 2))

  if (Type1Adj == "PFDR") {  ## print adjusted error rate
    cat("# -- Adjusted Type I Error Rates for Invariance Tests (Adjusted by Possible False Discovery Rate) -- #", rep("\n", 2))
    for (factor.no in 1: no.factor) {
      no.adjust <- (no.items[factor.no]-1)*(no.group-1)
      ERate[factor.no] <- Type1 / no.adjust
      cat(paste0("Factor ", factor.no, ": Number of items = ", no.items[factor.no], "; number of free factor loadings = ", no.adjust, "; adjusted error rate = ", sprintf("%.4f", ERate[factor.no]), "\n"))
    } # end (for factor.no)
  } else if (Type1Adj == "BON") {
    cat("# -- Adjusted Type I Error Rates for Invariance Tests (Bonforroni adjustment for number of pairwise comparisons) -- #", rep("\n", 2))
    for (factor.no in 1: no.factor) {
      no.adjust <- (no.items[factor.no]-1)*(no.group*(no.group-1)/2)
      ERate[factor.no] <- Type1 / no.adjust
      cat(paste0("Factor ", factor.no, ": Number of items = ", no.items[factor.no], "; number of pairwise comparisons = ", no.adjust, "; adjusted error rate = ", sprintf("%.4f", ERate[factor.no]), "\n"))
    } # end (for factor.no)
  } else if (Type1Adj == "NULL") {
    cat("# -- Error Rates for Invariance Tests (No adjustment is made) -- #", rep("\n", 2))
    cat(paste0("Type I Error Rate for Invariance Tests = ", Type1, "\n"))
    for (factor.no in 1: no.factor) {
      ERate[factor.no] <- Type1
    } # end (for factor.no)
  } # end (Type1Adj)
  cat(rep("\n",2))

  if (TYPE == "MonteCarlo") {

    mcmc <- MASS::mvrnorm(n=2000000, mu=par.est, Sigma=simvcov, tol = 1e-6)  # Run 2,000,000 simulations
    bootcoef <- mcmc
    bootno <- nrow(mcmc)  # No. of successful simulated samples
    cat(paste0("Number of Successful Simulated Samples = ", bootno, "\n"))

  } else { # TYPE=Bootstrap

    ## == Simplified bootstrapping model == ##
    Model.Config.Boot <- lavaan::cfa(Model.Long, data = data.source,
                            meanstructure = TRUE,
                            auto.fix.first = FALSE,
                            marker.int.zero = TRUE,
                            ordered = FALSE,
                            missing = "fiml", se = "none", test = "none", check.start = FALSE, check.post = FALSE, check.gradient = FALSE)

    ## == Bootstrapping == ##
    bootcoef <- lavaan::bootstrapLavaan(Model.Config.Boot, R = b.no, ncpus = max(1L, parallel::detectCores() - 1L), parallel="snow")

    # Remove error and nonadmissible bootstrap samples #
    B.na <- attr(bootcoef,"nonadmissible")
    B.er <- attr(bootcoef,"error.idx")
    B.de <- c(B.na,B.er)
    if (length(B.de) != 0) {
      B.re <- bootcoef[-c(B.de),]
      bootcoef <- B.re
    }
    bootno <- nrow(bootcoef)  # number of successful bootstrap samples
    cat(paste0("Number of Successful Bootstrapped Samples = ", bootno, "\n"))
    bootcoef <- bootcoef[,ext]

  } ## end MonteCarlo or Bootstrap


  ## == Start the factor.no loop for CompareLoadings == ##
  for (factor.no in 1: no.factor) {

    alpha <- ERate[factor.no]

    flY <<- matrix(" ",1, (no.group+2))
    flYY <<- matrix(" ",1, (no.group+2))
    EP <<- 1  # estimated parameter number

    no.k <- no.items[factor.no]*no.group  # number of factor loadings in FL.PMI
    FL.PMI <- matrix(1:no.k, nrow = no.group)  # factor loading matrix
    colnames(FL.PMI) <- paste0(' Item ', 1:no.items[factor.no])
    rownames(FL.PMI) <- c(paste0("Time ", 1:no.group))

    for (i in 1:no.group) {
      for (j in 1:no.items[factor.no]) {
        if (j == no.markers[factor.no]) {
          FL.PMI[i, j] <- 1
        } else {
          FL.PMI[i, j] <- format(round(par.est[paste0(names.lv[factor.no],"_T",i,"=~", names.item[[factor.no]][[j]],"_T",i)], digits=4), scientific = FALSE)
        }  ## end if j
      }  ## end loop j
    }  ## end loop i

    class(FL.PMI) <- "numeric"

    for (Referent in 1:no.items[factor.no]) {  ## Loop referent item number for comparison

      ## ==  Print Factor Loadings of all groups  == ##
      no.k <- no.items[factor.no]*no.group  # number of factor loadings in FL
      FL <- matrix(1:no.k, nrow = no.group)  # factor loading matrix
      colnames(FL) <- paste0(' Item ', 1:no.items[factor.no])
      rownames(FL) <- c(paste0("Time ", 1:no.group))

      if (Referent == 1) {  ## if referent is the first item
        for (i in 1:no.group) {
          for (j in 1:no.items[factor.no]) {
            if (j == Referent) {
              FL[i, j] <- 1
            } else {
              FL[i, j] <- FL.PMI[i, j]
            }  ## end if j
          }  ## end for j
        }  ## end for i
      } else {  ## referent is not the first item
        for (i in 1:no.group) {
          for (j in 1:no.items[factor.no]) {
            if (j == Referent) {
              FL[i, j] <- 1
            } else {
              FL[i, j] <- FL.PMI[i, j]/FL.PMI[i,Referent]
            }  ## end if j
          }  ## end for j
        }  ## end for i
      }  ## end if referent

#$      cat(rep("\n", 3), "Factor Loadings ", rep("\n", 2))
#$      print(round(FL[], digits=4))

    ## == End (Print Factor Loadings)  == ##

      for (nArg in 1: no.items[factor.no]) {  ## Argument item number for comparison
        no.dif <- factorial(no.group) / factorial(no.group - 2) / 2  # No. of pairwise comparisons for each item
        boot.dif.lx <- matrix(0, nrow(bootcoef), no.dif)  # Create bootstrap difference matrix
        samp.dif.lx <- matrix(0, no.group, no.group) # Create sample difference matrix
        if (nArg == Referent) {
          next
        } else {
          Arg <- nArg
        } ## end if nArg


      ## == Calculate bootstrap difference and sample estimate difference == ##
        comp = 0
        for (r in 1:(no.group-1)) {  ## r is referent group
          for (a in (r+1):no.group) {  ## a is argument group
            comp = comp + 1
            if (Referent == 1) {
              boot.dif.lx[,comp] <- bootcoef[, paste0(names.lv[factor.no],"_T", r, "=~", names.item[[factor.no]][[Arg]],"_T", r)] -
                                    bootcoef[, paste0(names.lv[factor.no],"_T", a, "=~", names.item[[factor.no]][[Arg]],"_T", a)]
              samp.dif.lx[r,a] <- par.est[paste0(names.lv[factor.no],"_T", r, "=~", names.item[[factor.no]][[Arg]],"_T", r)] -
                                  par.est[paste0(names.lv[factor.no],"_T", a, "=~", names.item[[factor.no]][[Arg]],"_T", a)]
            } else {
              if (Arg == no.markers[factor.no]) {
                boot.dif.lx[,comp] <- 1/bootcoef[, paste0(names.lv[factor.no],"_T", r, "=~", names.item[[factor.no]][[Referent]],"_T", r)] -
                                      1/bootcoef[, paste0(names.lv[factor.no],"_T", a, "=~", names.item[[factor.no]][[Referent]],"_T", a)]
                samp.dif.lx[r,a] <- 1/par.est[paste0(names.lv[factor.no],"_T", "=~", names.item[[factor.no]][[Referent]],"_T", r)] -
                                    1/par.est[paste0(names.lv[factor.no],"_T", "=~", names.item[[factor.no]][[Referent]],"_T", a)]
              } else { # # Arg != no.markers[factor.no]
                boot.dif.lx[,comp] <- bootcoef[,paste0(names.lv[factor.no],"_T", r, "=~", names.item[[factor.no]][[Arg]],"_T", r)]/
                                      bootcoef[,paste0(names.lv[factor.no],"_T", r, "=~", names.item[[factor.no]][[Referent]],"_T", r)] -
                                      bootcoef[,paste0(names.lv[factor.no],"_T", a, "=~", names.item[[factor.no]][[Arg]],"_T", a)]/
                                      bootcoef[,paste0(names.lv[factor.no],"_T", a, "=~", names.item[[factor.no]][[Referent]],"_T", a)]
                samp.dif.lx[r,a] <- par.est[paste0(names.lv[factor.no],"_T", r, "=~", names.item[[factor.no]][[Arg]],"_T", r)]/
                                    par.est[paste0(names.lv[factor.no],"_T", r, "=~", names.item[[factor.no]][[Referent]],"_T", r)] -
                                    par.est[paste0(names.lv[factor.no],"_T", a, "=~", names.item[[factor.no]][[Arg]],"_T", a)]/
                                    par.est[paste0(names.lv[factor.no],"_T", a, "=~", names.item[[factor.no]][[Referent]],"_T", a)]
              }  ## end if Arg
            }  ## end if Referent
          }  ## end loop a
        }  ## end loop r

#$ print(samp.dif.lx)


      ## == Calculate Percentile Probability == ##

      F1.comp.pp <- matrix(0, no.group, no.group)  ## Percentile Probability
      pno_nipair <- 0  ## Set number of non-invariant pairs to zero

      comp = 0
      for (r in 1:(no.group-1)) {  ## r is the referent group
        for (a in (r+1):no.group) {  ## a is the argument
          comp = comp + 1
          if (quantile(boot.dif.lx[, comp], probs = 0.5, na.rm = TRUE) > 0) {
            F1.comp.pp[r,a] = 2*(sum(boot.dif.lx[, comp] < 0, na.rm=TRUE)/bootno)
          } else {
            F1.comp.pp[r,a] = 2*(sum(boot.dif.lx[, comp] > 0, na.rm=TRUE)/bootno)
          }  ## end if
        }  ## end loop a
      }  ## end loop r


      ## == Percentile Confidence Intervals of Pairwise Comparisons == ##

      PCI <- matrix(1:(no.dif*10), nrow = no.dif)
      colnames(PCI) <- c("Group 1","Group 2","0.5%","2.5%","5%","Estimate","95%","97.5%","99.5%","p-value")

      comp = 0
      for (r in 1:(no.group-1)) {
        for (a in (r+1):no.group) {
          comp = comp + 1
          PCI[comp, 1] <- r
          PCI[comp, 2] <- a
          PCI[comp, 3] <- quantile(boot.dif.lx[, comp],c(0.005), na.rm = TRUE)
          PCI[comp, 4] <- quantile(boot.dif.lx[, comp],c(0.025), na.rm = TRUE)
          PCI[comp, 5] <- quantile(boot.dif.lx[, comp],c(0.05), na.rm = TRUE)
          PCI[comp, 6] <- samp.dif.lx[r,a]
          PCI[comp, 7] <- quantile(boot.dif.lx[, comp],c(0.95), na.rm = TRUE)
          PCI[comp, 8] <- quantile(boot.dif.lx[, comp],c(0.975), na.rm = TRUE)
          PCI[comp, 9] <- quantile(boot.dif.lx[, comp],c(0.995), na.rm = TRUE)
          PCI[comp,10] <- F1.comp.pp[r,a]
          if (F1.comp.pp[r,a] < alpha) {
            if (pno_nipair == 0)	{
              pnipair <- c(r,a)
              pno_nipair <- pno_nipair + 1
            } else {
              pnipair <- c(pnipair,r,a)
              pno_nipair <- pno_nipair + 1
            }  ## end if pno_nipair
          }  ## end if F1.comp.pp
        }  ## end loop a
      }  ## end loop r

#$      cat("\n")
#$      if (Referent > Arg) {
#$        cat(paste0("Referent Item: ", Referent, "  Argument Item: ", Arg, "\n"))
#$      } else {
#$        cat(paste0("Referent Item: ", Referent, "  Argument Item: ", (Arg+1), "\n"))
#$      }  ## end if Referent

#$      cat("\n")

#$    print(paste0("Factor = ", factor.no, "  Referent = ", Referent, " Argument = ", Arg))
#$    print(round(PCI[], digits=4))

      ## == Run List and Delete == ##
#$      cat(rep("\n", 2), "Sets of Groups with Invariant Loadings", "\n")
        if (pno_nipair == no.dif) {

          flZ <<- matrix(" ", 1, (no.group+2))
          flX <<- matrix(1:(no.group+2), 1, (no.group+2))
          flX[1, 1:no.group] <- 0
          for (temp.flZ in 1: no.group) {
            flZ[1,temp.flZ] <- paste0("F", factor.no, "L", EP)
            EP <<- EP + 1
          }
          flX[1, (no.group+1)] <- Referent
          flZ[1, (no.group+1)] <- Referent
#?          if (Referent > Arg) {
            flX[1, (no.group+2)] <- Arg
            flZ[1, (no.group+2)] <- Arg
#?          } else {
#?            flX[1, (no.group+2)] <- Arg + 1
#?            flZ[1, (no.group+2)] <- Arg + 1
#?          }  ## end if Referent
#$        print(flX, quote=FALSE)
          flY <<- rbind(flY,flX)
#$        print(flZ, quote=FALSE)
          flYY <<- rbind(flYY,flZ)

        } else if (pno_nipair > 0) {

          flX <<- matrix(" ", 1, (no.group+2))
          flZ <<- matrix(" ", 1, (no.group+2))

          listanddelete.lx(factor.no, no.group, no_nipair=pno_nipair, nipair=pnipair, Referent, Arg) ## Run list and delete lx

        } else {  # (pno_nipair = 0)

          flX <- matrix(1:(no.group+2), 1, (no.group+2))
          flZ <<- matrix(" ", 1, (no.group+2))
          flZ <- matrix(paste0("F", factor.no, "L", EP), 1, (no.group+2))
          EP <<- EP + 1
          flX[1, (no.group+1)] <- Referent
          flZ[1, (no.group+1)] <- Referent
          flX[1, (no.group+2)] <- Arg
          flZ[1, (no.group+2)] <- Arg
          flY <<- rbind(flY,flX)
#$        print(flZ, quote=FALSE)
          flYY <<- rbind(flYY,flZ)

        }  ## end if pno_nipair
      } ## end loop Arg
    } ## end loop Referent



  if (length(unique(flY[, no.group+1])) == 2 & no.items[factor.no] == 2) { # Only 2 groups with more than 2 items
    flX <<- matrix(1:(no.group+2), 1, (no.group+2))
    flZ <<- matrix(" ", 1, (no.group+2))
    flYa <<- flY[nrow(flY), no.group+1]
    flYb <<- flY[nrow(flY), no.group+2]
    flX[1, no.group + 1] <- flYb
    flX[1, no.group + 2] <- flYa
    flY <<- rbind(flY,flX)

    flZ[1, no.group + 1] <- flYb
    flZ[1, no.group + 2] <- flYa
    for (temp.flYY in 1: no.group) {
      flZ[1, temp.flYY] <- paste0("F", factor.no, "I", EP)
      EP <<- EP + 1
    }
    flYY <<- rbind(flYY,flZ)
  }

  flY <<- flY[-c(1),]
#$  cat("\n")
#$  print(flY, quote = FALSE)

  flYY <<- flYY[-c(1),]
#$  cat("\n")
#$  print(flYY, quote = FALSE)

  if (is.matrix(flYY) == FALSE) { flYY <<- matrix(flYY, nrow=1) }

  NIcombine <- table(flY[, no.group+1], flY[, no.group+2])
#  NI.level <- unique(flY[, no.group+1])
#  NIcombine <- NIcombine[NI.level, NI.level]

  tempR <- unique(flYY[, no.group+1])

    for (aR in 1:length(tempR)) {
      R <- tempR[aR]
      for (r in 1:no.items[factor.no]) {
        if (r == 1) {
          assign(paste0("Rset", r), c(length(which(flYY[,no.group+1] < R))+1:(NIcombine[R, r])))
        } else if (r == R) {
          assign(paste0("Rset", r), 0)
        } else {
          assign(paste0("Rset", r), c(length(which(flYY[,no.group+1] < R))+sum(NIcombine[R, 1:(r-1)])+1:(NIcombine[R, r])))
        }
        if (R == 1) { assign("Rset1", 0) }
      }  ## end loop r

      ## expand.grid(Rset1, Rset2, Rset3, Rset4, Rset5 ...)

      temp.row <- eval(parse(text = paste0("expand.grid(", paste0(c("Rset"), 1: no.items[factor.no], collapse=","),")")))
      if (aR == 1) {
        Set.row <- temp.row
      } else {
        Set.row <- rbind(Set.row, temp.row)
      } ## end if R
    }  ## end loop aR

    Model.load <- array(" ", dim = c(no.items[factor.no], no.group, nrow(Set.row)))  # all possible invariance models

    for (Xset in 1: nrow(Set.row)) {
      a <- matrix(" ", no.items[factor.no], no.group)
      for (r in 1:no.items[factor.no]) {
        if (Set.row[Xset, r] == 0) {
          a[r, ] <- 1
        } else {
          a[r, ] <- flYY[Set.row[Xset, r], 1:no.group]
        } ## end if Set.row
#?        if (flYY[r, no.group+1] > 1) { break }
      }  ## end loop r
      Model.load[, , Xset] <- a
    }  ## end loop Xset

#$ print(Model.load)

  ## == Model.summary == ##
  # Column 1 - number of estimated loadings
  # Column 2 - number of measurement invariance items
  # Column 3 - p-value of Model chi-square
  # Column 4 - CFI
  # Column 5 - RMSEA
  # Column 6 - SRMR

    Model.summary <- matrix(0, nrow(Set.row), 7)
    for (R in 1:nrow(Set.row)) {
      Model.summary[R, 1] <- length(unique(c(Model.load[,,R]))) - 1
      no.MI <- no.items[factor.no] - 1
      for (Ra in 1: no.items[factor.no]) {
        if (length(unique(c(Model.load[Ra,,R]))) > 1) {no.MI <- no.MI - 1}
      }  ## end loop Ra
      Model.summary[R, 2] <- no.MI
    }  ## end loop R



    ## == Recommended Models (Rec.Model) - Largest no. of MI items, then smallest number of estimated loadings == ##
    Rec.Model <- which(Model.summary[,2] == max(Model.summary[,2]))
    Rec.Model <- Rec.Model[which(Model.summary[Rec.Model,1] == min(Model.summary[Rec.Model,1]))]

    ## == Find recommended model if no.items > 2 == ##
    if (no.items[factor.no] > 2) {

      Rec.Model.load <- matrix(nrow = length(Rec.Model), ncol = length(lavaan::coef(Model.Long.config)))
      MODFIT <- matrix(0, nrow = length(Rec.Model), ncol = 9)
      colnames(MODFIT) <-
        c("  chisq.scaled", "  df.scaled", "  pvalue.scaled", "  rmsea.scaled", "  cfi.scaled", "  tli.scaled", "  srmr", "  aic", "  bic")
      Rec.Model.FL <- array(0, dim = c(no.group, no.items[factor.no], length(Rec.Model)))
      colnames(Rec.Model.FL) <- paste0(" Item ", 1:no.items[factor.no])
      rownames(Rec.Model.FL) <- c(paste0("Time ", 1:no.group))

      for (Model.R in 1: length(Rec.Model)) {

        ## == Specify the partial metric invariance models (PMI.X) == ##

        final.load <- unique(as.vector(Model.load[,, Rec.Model[Model.R]]))
        for (j in 1:no.waves) {
          for (i in 1:no.items[factor.no]) {
            if (Model.load[i,j,Rec.Model[Model.R]] != 1) { Model.load[i,j, Rec.Model[Model.R]] <- paste0("F", factor.no, "L", which(Model.load[i,j,Rec.Model[Model.R]] == final.load)) }
          } # end for i
        } # end for j

        PMI.X.1 <- paste0("PMI.X <- '", "\n")
        for (j in 1:(no.waves)) {
          PMI.X.1 <- paste0(PMI.X.1, "  ",
                          names.lv[factor.no], "_T",j," =~ ", paste0(Model.load[1, j, Rec.Model[Model.R]]), "*", names.item[[factor.no]][[1]], "_T",j," + ")
          for (i in 2: (no.items[factor.no]-1)) {
            PMI.X.1 <- paste0(PMI.X.1, Model.load[i, j, Rec.Model[Model.R]],"*", names.item[[factor.no]][[i]], "_T",j," + ")
          } ## end for i
          PMI.X.1 <- paste0(PMI.X.1, Model.load[i+1, j, Rec.Model[Model.R]],"*", names.item[[factor.no]][[i+1]], "_T",j, "\n")
        } # end for j
        ## -- Create Residual Covariance -- ##

        PMI.X.1 <- paste0(PMI.X.1, "\n")
        for (i in 1:no.items[factor.no]) {
          for (j in 1:(no.waves - 1)) {
            for (k in (j+1):no.waves) {
              PMI.X.1 <- paste0(PMI.X.1, "  ", names.item[[factor.no]][[i]], "_T", j, " ~~ ", names.item[[factor.no]][[i]], "_T", k, "\n")
            } # end for k
          } # end for j
        } # end for i

        PMI.X.1 <- paste0(PMI.X.1, "'", "\n")
        eval(parse(text = PMI.X.1))


        ## == Run PMI.X == ##
        if (Cluster == "NULL") {
          eval(parse(text=   "PMI.X.fit <- lavaan::sem(PMI.X,
                             data.source,
                             missing = 'fiml',
                             auto.fix.first = FALSE,
                             marker.int.zero = TRUE,
                             meanstructure = T,
                             information = 'observed',
                             estimator = 'MLR')"))
        } else {
          eval(parse(text=   "PMI.X.fit <- lavaan::sem(PMI.X,
                             data.source,
                             missing = 'fiml',
                             auto.fix.first = FALSE,
                             marker.int.zero = TRUE,
                             meanstructure = T,
                             information = 'observed',
                             cluster = Cluster,
                             estimator = 'MLR')"))
        }


        ## == Request summary outputs == ##
 #$     eval(parse(text = "print(lavaan::summary(PMI.X.fit, fit.measure = T, standardized = T, rsq = T))"))

        ## == Save fit indices == ##
        XX <- lavaan::fitMeasures(PMI.X.fit, fit.measures=
              c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr", "aic", "bic"))
        MODFIT[Model.R,] <- format(round((XX), digits = 4), nsmall = 4, scientific = FALSE)


        aa <- lavaan::lavInspect(PMI.X.fit, what='est')$lambda
        ## == Save Factor Loadings == ##
        for (G in 1: no.group) {
          Rec.Model.FL[G, , Model.R] <- aa[(no.items[factor.no]*(G-1)+1):(no.items[factor.no]*G),G]
        } ## end loop G

      } ## end loop Model.R


      ## == Print model fit indices == ##
      rownames(MODFIT) <- c(paste0("Model ", 1: length(Rec.Model)))
#$    cat(rep("\n", 3))
#$    cat("####    MODEL FIT INDICES   ####", rep("\n", 2))
#$    print(MODFIT, quote=FALSE, right=TRUE)

    ## == Print factor loadings == ##
#$    for (Model.R in 1: length(Rec.Model)) {
#$      cat(rep("\n", 2), "Factor Loadings of Model", Model.R, rep("\n", 2))
#$      print(round(Rec.Model.FL[,,Model.R], digits=4))
#$    }  ## end loop Model.R

      ## == First model with best model fit (R.Model) == ##
      class(MODFIT) <- "numeric"
      if (BMSC == "Chisquare") {
        Model.R <- which.min(MODFIT[,1])
      } else if (BMSC == "CFI") {
        Model.R <- which.max(MODFIT[,5])
      } else if (BMSC == "RMSEA") {
        Model.R <- which.min(MODFIT[,4])
      } else if (BMSC == "SRMR") {
        Model.R <- which.min(MODFIT[,7])
      }
      Model.R <- Model.R[1]
      S.Model <- Model.load[,,Rec.Model[Model.R]]  # Recommended Model
      un.load <- sort(subset(unique(c(S.Model)), unique(c(S.Model))[] != 1))
      no.un.load <- length(un.load)
      R.Model <- S.Model
      for (j in 1:no.un.load) { R.Model[S.Model == un.load[j]] <- paste0("F", factor.no, "L", j) }  ## end loop j

    ## == Save Recommended Model (Recommend.Model) == ##
      if (factor.no == 1) {
        Recommend.Model <- matrix("PMI.Model.R <- '", 1)
      }
      Recommend.Model <- rbind(Recommend.Model, "\n")

      for (j in 1:no.waves) {
        PMI <- paste0("   ", names.lv[factor.no], "_T",j," =~ ", paste0(Model.load[1, j, Rec.Model[Model.R]]), "*", names.item[[factor.no]][[1]], "_T",j," + ")
        for (i in 2: (no.items[factor.no]-1)) {
          PMI <- paste0(PMI, R.Model[i,j],"*", names.item[[factor.no]][[i]], "_T",j," + ")
        } ## end for i
        PMI <- paste0(PMI, R.Model[i+1,j],"*", names.item[[factor.no]][[i+1]], "_T",j, "\n")
      Recommend.Model <- rbind(Recommend.Model, PMI)
      } # end for j

      ## -- Create Residual Covariance -- ##
      Recommend.Model <- rbind("\n", Recommend.Model, "\n")
      for (i in 1:no.items[factor.no]) {
        for (j in 1:(no.waves - 1)) {
          for (k in (j+1):no.waves) {
            Recommend.Model <- rbind(Recommend.Model, paste0("   ", names.item[[factor.no]][[i]], "_T", j, " ~~ ", names.item[[factor.no]][[i]], "_T", k, "\n"))
          } # end for k
        } # end for j
      } # end for i

    }  ## end if no.items[factor.no] > 2


    if (no.items[factor.no] < 3) {   ## Begin Models with 2 items ##
      R.Model <- Model.load[,,1]  # Recommended Model
      un.load <- subset(unique(c(R.Model)), unique(c(R.Model))[] != 1)
      no.un.load <- length(un.load)
      for (j in 1:no.un.load) { R.Model[R.Model == un.load[j]] <- paste0("F", factor.no, "L", j) }  ## end loop j

      ## == Save Recommended Model (Recommend.Model) == ##
      if (factor.no == 1) {
        Recommend.Model <- matrix("PMI.Model.R <- '", 1)
      }
      Recommend.Model <- rbind(Recommend.Model, "\n")

      for (j in 1:no.waves) {
        PMI <- paste0("   ", names.lv[factor.no], "_T",j," =~ ", paste0(Model.load[1, j, Rec.Model[Model.R]]), "*", names.item[[factor.no]][[1]], "_T",j," + ")
        PMI <- paste0(PMI, R.Model[2,j],"*", names.item[[factor.no]][[2]], "_T",j, "\n")
      Recommend.Model <- rbind(Recommend.Model, PMI)
      } # end for j

      ## -- Create Residual Covariance -- ##
      Recommend.Model <- rbind("\n", Recommend.Model, "\n")
      for (i in 1:2) {
        for (j in 1:(no.waves - 1)) {
          for (k in (j+1):no.waves) {
            Recommend.Model <- rbind(Recommend.Model, paste0("   ", names.item[[factor.no]][[i]], "_T", j, " ~~ ", names.item[[factor.no]][[i]], "_T", k, "\n"))
          } # end for k
        } # end for j
      } # end for i
      Recommend.Model <- rbind(Recommend.Model, "\n")

    }  ## End (Models with 2 items)

  } ## End (loop factor.no)

  Recommend.Model <- rbind(Recommend.Model, "  '", "\n")  ## last line of Recommended Model

  cat("\n", "## =====  Recommended Model  ===== ##", rep("\n", 2))
  cat(Recommend.Model)
  cat(rep("\n", 2))
  cat("## ===== lavaan Outputs of PMI.Model.R ===== ##")
  cat(rep("\n", 2))


  ## == Start print recommended model to PMI.txt (Partial Metric Invariance Model) == ##
  sink('PMI.txt')

  cat("## =====  Recommended Model  ===== ##")
  cat("\n")
  cat(Recommend.Model)
  cat(rep("\n",2), "## Run PMI.Model.R")
  cat(rep("\n",2), paste0("  PMI.Model.fit <- lavaan::sem(PMI.Model.R,"))
  cat("\n", paste0("    ", arg2_char, ","))
  cat("\n", "    missing = 'fiml',")
  cat("\n", "    auto.fix.first = FALSE,")
  cat("\n", "    information = 'observed',")
  if (Cluster != "NULL") { cat("\n", paste0("    cluster = ", arg4_char,",")) }
  cat("\n", "    estimator = 'MLR')")

  cat(rep("\n", 3), "## Request summary outputs")
  cat(rep("\n", 3), paste0("  print('## ===== Final Model -- PMI.Model.R ===== ##')", "\n"))
  cat(rep("\n", 2), paste0("  print(lavaan::summary(PMI.Model.fit, fit.measure = T, standardized = T, rsq = T))"), "\n")

  sink()  ## Stop writing to file

  PMI.Model.R <<- eval(parse(text=Recommend.Model))

  ## == Run PMI.Model.R == ##
  if (Cluster == "NULL") {
    eval(parse(text=   "PMI.Model.fit <- lavaan::sem(PMI.Model.R,
                       data.source,
                       missing = 'fiml',
                       auto.fix.first = FALSE,
                       marker.int.zero = TRUE,
                       meanstructure = T,
                       information = 'observed',
                       estimator = 'MLR')"))
  } else {
    eval(parse(text=   "PMI.Model.fit <- lavaan::sem(PMI.Model.R,
                       data.source,
                       missing = 'fiml',
                       auto.fix.first = FALSE,
                       marker.int.zero = TRUE,
                       meanstructure = T,
                       information = 'observed',
                       cluster = Cluster,
                       estimator = 'MLR')"))
  }

  print(lavaan::summary(PMI.Model.fit, fit.measure = T, standardized = T, rsq = T))


  cat(rep("\n", 2), "## Compare fit indices across models ##", "\n")
  DMODFIT <- matrix(0, nrow = 3, ncol = 9)
  colnames(DMODFIT) <-
        c("  chisq.scaled", "  df.scaled", "  pvalue.scaled", "  rmsea.scaled", "  cfi.scaled", "  tli.scaled", "  srmr", "  aic", "  bic")
  rownames(DMODFIT) <- c("Model.config", "PMI.Model.fit", "PMI.Model.fit - Model.config")
  XX1 <- lavaan::fitMeasures(Model.Long.config, fit.measures=
      c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr", "aic", "bic"))
  DMODFIT[1,] <- XX1
  XX2 <- lavaan::fitMeasures(PMI.Model.fit, fit.measures=
      c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr", "aic", "bic"))
  DMODFIT[2,] <- XX2

  DMODFIT[3,1] <- DMODFIT[2,1] -DMODFIT[1,1]
  DMODFIT[3,2] <- DMODFIT[2,2] -DMODFIT[1,2]
  DMODFIT[3,3] <- pchisq(q=DMODFIT[3,1], df=DMODFIT[3,2], lower.tail=FALSE)
  DMODFIT[3,4] <- DMODFIT[2,4] -DMODFIT[1,4]
  DMODFIT[3,5] <- DMODFIT[2,5] -DMODFIT[1,5]
  DMODFIT[3,6] <- DMODFIT[2,6] -DMODFIT[1,6]
  DMODFIT[3,7] <- DMODFIT[2,7] -DMODFIT[1,7]
  DMODFIT[3,8] <- DMODFIT[2,8] -DMODFIT[1,8]
  DMODFIT[3,9] <- DMODFIT[2,9] -DMODFIT[1,9]
  print(round(DMODFIT, digits=4))


  ## == Print Factor Loadings == ##
  Final.Model.FL <- matrix(0, sum(no.items), no.waves)
  rownames(Final.Model.FL) <- paste0(" Item ", names.ind, "    ")
  colnames(Final.Model.FL) <- c(paste0("Time ", 1:no.group))

  Final.par.est <- lavaan::parameterEstimates(PMI.Model.fit)
  ext <- c(which(Final.par.est[,"op"] == "=~"))
  aa <- Final.par.est[ext, "est"]

  ## Extract factor loadings ##
  n <- 1
  for (j in 1: no.waves) {
    Final.Model.FL[1:no.items[1], j] <- aa[(no.items[1]*(j-1)+1):(no.items[1]*j)]
  }
  n <- n + no.items[1]*no.waves

  if(no.factor > 1){
    for (i in 2:no.factor) {
      for (j in 1: no.waves) {
        Final.Model.FL[(sum(no.items[1:(i-1)])+1):sum(no.items[1:i]), j] <- aa[(n+no.items[i]*(j-1)):(n+(no.items[i]*j-1))]
      }
      n <- n + no.items[i]*no.waves
    }
  }

  cat("\n")
  cat("## === Factor Loadings in Final Model === ##")
  cat("\n")
  print(round(Final.Model.FL, digits=4))
  ## == End (Print Factor Loadings) == ##

  cat(rep("\n",2),"The recommended model PMI.Model.R is saved in the file 'PMI.txt'", "\n")

} ## End (Function LGCompareLoadings)

# ==================== Finish Function "LGCompareLoadings" ==================== #






# ==================== Create Function "LGCompareMeans" ==================== #
#' Scalar Invariance Test and Compare Latent Means for Longitudinal Panel Data
#'
#' Conduct scalar invariance test, identify a partial scalar invariance model, and compare latent means for longitudinal panel data.
#'
#' Requires defining the measurement model only once. All indicators should have the suffix _T1, _T2 and _T3 (e.g., x1_T1, x1_T2, x1_T3) in the data file to indicate when the indicators were measured.
#'
#' Residuals of indicators are covaried across time automatically.
#'
#' Reference for correct.cfi: Widaman, K. F., & Thompson, J. S. (2003). On specifying the null model for incremental fit indices in structural equation modeling. Psychological Methods, 8(1), 16-37.
#'
#' @param model.PMI Partial metric invariance model (PMI.Model.R from LGCompareLoadings() or user-specified).
#' @param data.source A data frame containing the observed variables used in the model.
#' @param no.waves Number of waves for comparisons.
#' @param Cluster Cluster variable for nested data. The Monte Carlo simulation method should be used for nested data.
#' @param Bootstrap Number of bootstrap samples, must be between 500 and 10,000. If not specified, 1 million Monte Carlo simulated samples (default) will be used.
#' @param Type1 Overall Type I error rate (default is 0.05) for identifying non-invariant items in the List and Delete method.
#' @param Type1Adj Adjustment of Type I error rate for multiple tests for each construct. Default is Possible False Discovery Rate (PFDR = Type I error rate/No. of freely estimated factor loadings/intercepts across all waves), can also use (BON)ferroni adjustment (Type I error rate/No. of pairwise comparisons), or "NULL".
#' @param BMSC Best Model Selection Criterion used to select among models with the same number of estimated parameters. "ChiSquare" = minimum Chi-Square value, "CFI" = maximum CFI value, "RMSEA" = minimum RMSEA value, and "SRMR" = minimum SRMR value (default).
#' @param correct.cfi If TRUE (default), baseline model will be corrected for calculation of incremental fit indices based on Widaman & Thompson (2003).
#'
#' @return Partial scalar invariance model in PSI.txt file and results of latent means comparisons.
#' @export
#' @examples
#'
#'
#' ## == Example B - Panel Data in Longitudinal Studies == ##
#'
#' # Data file is "Example.B"
#'
#' ## Not run:
#' ## Specify the measurement model - Model.B ##
#' Model.B <- '
#'   SWLS =~ x1 + x2 + x3 + x4 + x5
#' '
#'
#' ## ===== Compare Factor Loadings ===== ##
#' LGCompareLoadings(Model.B, Example.B, no.waves=3, Type1=0.05, Type1Adj="PFDR", BMSC="SRMR")
#' ## End (Not run)
#'
#' ## ===== Compare Means ===== ##
#' LGCompareMeans(PMI.Model.R, Example.B, no.waves=3, Type1=0.05, Type1Adj="PFDR", BMSC="SRMR")
#'
#'
LGCompareMeans <- function(model.PMI, data.source, Cluster="NULL", no.waves=3, Bootstrap=0, Type1 = 0.05, Type1Adj = "PFDR", BMSC="SRMR", correct.cfi=TRUE) {

  Type1Adj <- toupper(Type1Adj)
  match.arg(Type1Adj, c("PFDR","BON", "NULL"))
  BMSC <- toupper(BMSC)
  match.arg(BMSC, c("CHISQUARE","CFI", "RMSEA", "SRMR"))
  if (is.logical(correct.cfi) == FALSE) stop("correct.cfi (corrected baseline model for calculating incremental fit indices) can only be TRUE or FALSE")

  arg1_char <- deparse(substitute(model.PMI))
  arg2_char <- deparse(substitute(data.source))
  arg4_char <- deparse(substitute(Cluster))

  no.group <- no.waves
  parsed <- lavParseModelString(model.PMI, as.data.frame = TRUE)
  parsed <- parsed[, c("lhs", "op", "rhs")]
  parsed <- parsed[parsed$op == "=~",]
  parsed[,"lhs"] <- gsub("_T.*", "", parsed[,"lhs"])
  parsed[,"rhs"] <- gsub("_T.*", "", parsed[,"rhs"])
  names.lv <- unique(parsed[,"lhs"])
  names.ind <- unique(parsed[, "rhs"])
  no.factor <- length(names.lv)
  no.items.g <- length(names.ind) # number of items per group
  no.items <- matrix(1:no.factor, nrow = 1)  # number of items per factor
  for (factor.no in 1:no.factor) { no.items[factor.no] <- sum(parsed[,"lhs"] == names.lv[factor.no] & parsed[,"op"] == "=~")/no.waves } ## end loop factor.no

  names.item <- list()
  for (i in 1:no.factor) {
    names.item[[i]] <- unique(parsed[which(parsed[,"lhs"] == names.lv[i] & parsed[,"op"] == "=~"), "rhs"])
  }


  ## Check for bootstrap sample number (Bootstrap) ##
  if (Bootstrap !=0) {
    b.no.integer <- Bootstrap == round(Bootstrap)
    if (b.no.integer == "FALSE") stop("Bootstrap sample number must be an integer")
    if (Bootstrap > 10000) stop("Bootstrap sample number greater than 10,000 is not recommended")
    if (Bootstrap < 500) stop("Bootstrap sample number smaller than 500 is not recommended")
    TYPE = "Bootstrap"
    b.no <- Bootstrap
  } else {
    TYPE = "MonteCarlo"
  } # end (Bootstrap != 0)

  count.tx <- 0


  ## ========== Run model.PMI ========== ##

  if (Cluster == "NULL") {
    PMI.Model.fit <- lavaan::sem(model.PMI,
      data.source,
      missing = 'fiml',
      auto.fix.first = FALSE,
      marker.int.zero = TRUE,
      meanstructure = T,
      information = 'observed',
      estimator = 'MLR')
  } else {
      PMI.Model.fit <- lavaan::sem(model.PMI,
      data.source,
      missing = 'fiml',
      auto.fix.first = FALSE,
      marker.int.zero = TRUE,
      meanstructure = T,
      information = 'observed',
      cluster = Cluster,
      estimator = 'MLR')
   }  # end Cluster

  ## Request summary outputs
  cat(rep("\n", 2), "## ===== Partial Metric Invariance Model ===== ##", "\n")
  print(lavaan::summary(PMI.Model.fit, fit.measure = T, standardized = T, rsq = T))
  cat("\n")

  ## ===== End (Run model.PMI) ===== ##

  par.est <- lavaan::coef(PMI.Model.fit)  # sample parameters
  group.names <<- c(paste0("Time ", 1:no.group))  # group names

  temp <- lavaan::parameterEstimates(PMI.Model.fit)
  no.markers <- matrix(1:no.factor, nrow = 1)  # location of marker items

  for (factor.no in 1:no.factor) {
    no.markers[factor.no] <- which(names.item[[factor.no]] == sub("_T1$", "", temp[which(temp[,"lhs"] == paste0(names.lv[factor.no],"_T1") & temp[,"op"] == "=~" & temp[,"est"] == 1), "rhs"]))
  } ## end for factor.no

  temp <- lavaan::parameterEstimates(PMI.Model.fit, remove.nonfree=TRUE)
  simvcov <- lavInspect(PMI.Model.fit, what="vcov")
  par.est <- lavaan::coef(PMI.Model.fit)  # sample parameters

  ## -- Change labels to names -- ##
  names(par.est) <- paste0(temp[,"lhs"], temp[,"op"], temp[,"rhs"])
  colnames(simvcov) <- paste0(temp[,"lhs"], temp[,"op"], temp[,"rhs"])
  rownames(simvcov) <- paste0(temp[,"lhs"], temp[,"op"], temp[,"rhs"])


  ## Extract factor loadings and indicator intercepts ##
  mext <- which(temp[,"op"] == "=~" | temp[, "op"] == "~1")
  par.est <- par.est[mext]
  simvcov <- simvcov[mext,mext]
  no.par.g <- length(par.est)/no.waves  # number of estimated parameters per group #
  no.lx.g <- sum(temp[, "op"] == "=~")/no.waves # number of estimated LX per group #


  cat(rep("\n", 3), "## ======= SCALAR INVARIANCE ANALYSIS ======= ##", rep("\n", 2))  ## print heading

  # Calculate alpha for factor[factor.no] #
  ERate <- matrix(1:no.factor, nrow = 1)  # error rate for each factor

  cat("# -- Overall Type I Error Rates for Invariance Tests for Each Factor = ", Type1, " -- #", rep("\n", 2))

  if (Type1Adj == "PFDR") {  ## print adjusted error rate
    cat("# -- Adjusted Type I Error Rates for Invariance Tests (Adjusted by Possible False Discovery Rate) -- #", rep("\n", 2))
    for (factor.no in 1: no.factor) {
      no.adjust <- (no.items[factor.no]-1)*(no.group-1)
      ERate[factor.no] <- Type1 / no.adjust
      cat(paste0("Factor ", factor.no, ": Number of items = ", no.items[factor.no], "; number of free intercepts = ", no.adjust, "; adjusted error rate = ", sprintf("%.4f", ERate[factor.no]), "\n"))
    } # end (for factor.no)
  } else if (Type1Adj == "BON") {
    cat("# -- Adjusted Type I Error Rates for Invariance Tests (Bonforroni adjustment for number of pairwise comparisons) -- #", rep("\n", 2))
    for (factor.no in 1: no.factor) {
      no.adjust <- (no.items[factor.no]-1)*(no.group*(no.group-1)/2)
      ERate[factor.no] <- Type1 / no.adjust
      cat(paste0("Factor ", factor.no, ": Number of items = ", no.items[factor.no], "; number of pairwise comparisons = ", no.adjust, "; adjusted error rate = ", sprintf("%.4f", ERate[factor.no]), "\n"))
    } # end (for factor.no)
  } else if (Type1Adj == "NULL") {
    cat("# -- Error Rates for Invariance Tests (No adjustment is made) -- #", rep("\n", 2))
    cat(paste0("Type I Error Rate for Invariance Tests = ", Type1, "\n"))
    for (factor.no in 1: no.factor) {
      ERate[factor.no] <- Type1
    } # end (for factor.no)
  } # end (Type1Adj)
  cat(rep("\n",2))

  if (TYPE == "MonteCarlo") {

    mcmc <- MASS::mvrnorm(n=2000000, mu=par.est, Sigma=simvcov, tol = 1e-6)  # Run 2,000,000 simulations
    bootcoef <- mcmc
    bootno <- nrow(mcmc)  # No. of successful simulated samples
    cat(paste0("Number of Successful Simulated Samples = ", bootno, "\n"))

  } else { # TYPE=Bootstrap

    ## == Simplified bootstrapping model == ##
    PMI.boot <- lavaan::cfa(model.PMI, data = data.source,
                meanstructure = TRUE,
                auto.fix.first = FALSE,
                marker.int.zero = TRUE,
                ordered = FALSE,
                missing = "fiml", se = "none", test = "none", check.start = FALSE, check.post = FALSE, check.gradient = FALSE)

    ## == Bootstrapping == ##
    bootcoef <- lavaan::bootstrapLavaan(PMI.boot, R = b.no, ncpus = max(1L, parallel::detectCores() - 1L), parallel="snow")

    # Remove error and nonadmissible bootstrap samples #
    B.na <- attr(bootcoef,"nonadmissible")
    B.er <- attr(bootcoef,"error.idx")
    B.de <- c(B.na,B.er)
    if (length(B.de) != 0) {
      B.re <- bootcoef[-c(B.de),]
      bootcoef <- B.re
    }
    bootno <- nrow(bootcoef)  # number of successful bootstrap samples
    cat(paste0("Number of Successful Bootstrapped Samples = ", bootno, "\n"))
    bootcoef <- bootcoef[,ext]

  } ## end MonteCarlo or Bootstrap


  ## == Start the factor.no loop for CompareIntercepts == ##
  for (factor.no in 1: no.factor) {

    alpha <- ERate[factor.no]

    txY <<- matrix(" ",1, (no.group+2))
    txYY <<- matrix(" ",1, (no.group+2))
    EP <<- 1  # estimated parameter number

    no.k <- no.items[factor.no]*no.group  # number of factor loadings in FL.PMI
    FL.PMI <- matrix(1:no.k, nrow = no.group)  # factor loading matrix
    colnames(FL.PMI) <- paste0(' Item ', 1:no.items[factor.no])
    rownames(FL.PMI) <- c(paste0("Time ", 1:no.group))

    TX.PMI <- matrix(1:no.k, nrow = no.group)  # Intercept matrix
    colnames(TX.PMI) <- paste0(' Item ', 1:no.items[factor.no])
    rownames(TX.PMI) <- c(paste0("Time ", 1:no.group))

    for (i in 1:no.group) {
      for (j in 1:no.items[factor.no]) {
        if (j == no.markers[factor.no]) {
          FL.PMI[i, j] <- 1
          TX.PMI[i, j] <- 0
        } else {
          FL.PMI[i, j] <- format(round(par.est[paste0(names.lv[factor.no],"_T",i,"=~", names.item[[factor.no]][[j]],"_T",i)], digits=4), scientific = FALSE)
          TX.PMI[i, j] <- format(round(par.est[paste0(names.item[[factor.no]][[j]],"_T", i, "~1")], digits=4), scientific = FALSE)
        }  ## end if j
      }  ## end loop j
    }  ## end loop i

    class(FL.PMI) <- "numeric"
    class(TX.PMI) <- "numeric"

    for (Referent in 1:no.items[factor.no]) {  ## Loop referent item number for comparison

      if(length(unique(FL.PMI[,Referent])) != 1) {
        txZ <<- matrix(" ", 1, (no.group+2))
        txX <<- matrix(1:(no.group+2), 1, (no.group+2))
        txX[1, 1:no.group] <- 0
        uni.fl <- unique(FL.PMI[, Referent])
        Luni.fl <- length(uni.fl)
        for (uni in 1: Luni.fl) {
          flfl <- uni.fl[uni]
          for (temp.txZ in 1: no.group) {
            if (txZ[1, temp.txZ] == flfl) {
              txZ[1,temp.txZ] <- paste0("F", factor.no, "I", EP)
            }
          }
          EP <<- EP + 1
        }
        txX[1, (no.group+1)] <- Referent
        txZ[1, (no.group+1)] <- Referent

        Arg <- 1

        if (Referent > Arg) {
          txX[1, (no.group+2)] <- Arg
          txZ[1, (no.group+2)] <- Arg
        } else {
          txX[1, (no.group+2)] <- Arg + 1
          txZ[1, (no.group+2)] <- Arg + 1
        }  ## end if Referent
#$      print(txX, quote=FALSE)
#        txY <<- rbind(txY,txX)
#$      print(txZ, quote=FALSE)
#        txYY <<- rbind(txYY,txZ)

        next
       } # skip if unequal factor loadings


      ## ==  Print Factor Loadings FL of all groups  == ##
      FL <- matrix(1:no.k, nrow = no.group)  # factor loading matrix
      colnames(FL) <- paste0(' Item ', 1:no.items[factor.no])
      rownames(FL) <- c(paste0("Time ", 1:no.group))

      if (Referent == no.markers[factor.no]) {  ## if referent is the marker item
        for (i in 1:no.group) {
          for (j in 1:no.items[factor.no]) {
            if (j == Referent) {
              FL[i, j] <- 1
            } else {
              FL[i, j] <- FL.PMI[i, j]
            }  ## end if j
          }  ## end for j
        }  ## end for i
      } else {  ## referent is not the marker item
        for (i in 1:no.group) {
          for (j in 1:no.items[factor.no]) {
            if (j == Referent) {
              FL[i, j] <- 1
            } else {
              FL[i, j] <- FL.PMI[i, j]/FL.PMI[i, Referent]
            }  ## end if j
          }  ## end for j
        }  ## end for i
      }  ## end if referent

#$      cat(rep("\n", 3), "Factor Loadings ", rep("\n", 2))
#$      print(round(FL[], digits=4))

      ## == End (Print Factor Loadings)  == ##

      ## ==  Print Intercepts of all groups  == ##
      no.k <- no.items[factor.no]*no.group  # number of intercepts in TX
      TX <- matrix(1:no.k, nrow = no.group)  # intercept matrix
      colnames(TX) <- paste0(' Item ', 1:no.items[factor.no])
      rownames(TX) <- c(paste0("Time ", 1:no.group))

      if (Referent == no.markers[factor.no]) {  ## if referent is the marker item
        for (i in 1:no.group) {
          for (j in 1:no.items[factor.no]) {
            if (j == Referent) {
              TX[i, j] <- 0
            } else {
              TX[i, j] <- TX.PMI[i, j]
            }  ## end if j
          }  ## end for j
        }  ## end for i
      } else {  ## referent is not the marker item
        for (i in 1:no.group) {
          for (j in 1:no.items[factor.no]) {
            if (j == Referent) {
              TX[i, j] <- 0
            } else {
              TX[i, j] <- TX.PMI[i, j] - FL.PMI[i, j]/FL.PMI[i, Referent]*TX.PMI[i, Referent]
            }  ## end if j
          }  ## end for j
        }  ## end for i
      }  ## end if referent

#$      cat(rep("\n", 3), "Intercepts ", rep("\n", 2))
#$      print(round(TX[], digits=4))

      ## == End (Print Intercepts)  == ##


      for (nArg in 1: no.items[factor.no]) {  ## Argument item number for comparison
        no.dif <- factorial(no.group) / factorial(no.group - 2) / 2  # No. of pairwise comparisons for each item
        boot.dif.tx <- matrix(0, nrow(bootcoef), no.dif)  # Create bootstrap difference matrix
        samp.dif.tx <- matrix(0, no.group, no.group) # Create sample difference matrix
        if (nArg == Referent) {
          next
        } else {
          Arg <- nArg
        } ## end if nArg

        ## == Calculate bootstrap difference and sample estimate difference == ##
        comp = 0
        for (r in 1:(no.group-1)) {  ## r is referent group
          for (a in (r+1):no.group) {  ## a is argument group
            comp = comp + 1
            if (Referent == no.markers[factor.no]) {
              boot.dif.tx[,comp] <- bootcoef[, paste0(names.item[[factor.no]][[Arg]],"_T", r, "~1")] -
                                    bootcoef[, paste0(names.item[[factor.no]][[Arg]],"_T", a, "~1")]
              samp.dif.tx[r,a] <- par.est[paste0(names.item[[factor.no]][[Arg]],"_T", r, "~1")] -
                                  par.est[paste0(names.item[[factor.no]][[Arg]],"_T", a, "~1")]
#$ print(paste0("Referent = ", Referent, " & Arg = ", Arg))
#$ print(par.est[paste0(names.item[[factor.no]][[Arg]],"_T", r, "~1")])
#$ print(par.est[paste0(names.item[[factor.no]][[Arg]],"_T", a, "~1")])
            } else if (Referent != no.markers[factor.no]) {
              if (Arg == no.markers[factor.no]) {
                boot.dif.tx[, comp] <- (-1*bootcoef[, paste0(names.item[[factor.no]][[Referent]],"_T", r, "~1")]/
                                           bootcoef[, paste0(names.lv[factor.no],"_T", r, "=~", names.item[[factor.no]][[Referent]],"_T", r)]) -
                                       (-1*bootcoef[, paste0(names.item[[factor.no]][[Referent]],"_T", a, "~1")]/
                                           bootcoef[, paste0(names.lv[factor.no],"_T", a, "=~", names.item[[factor.no]][[Referent]],"_T", a)])
                samp.dif.tx[r,a] <- (-1*par.est[paste0(names.item[[factor.no]][[Referent]],"_T", r, "~1")]/
                                        par.est[paste0(names.lv[factor.no],"_T", r, "=~", names.item[[factor.no]][[Referent]],"_T", r)]) -
                                    (-1*par.est[paste0(names.item[[factor.no]][[Referent]],"_T", a, "~1")]/
                                        par.est[paste0(names.lv[factor.no],"_T", a, "=~", names.item[[factor.no]][[Referent]],"_T", a)])
#$ print(paste0("Referent = ", Referent, " & Arg = ", Arg))
#$ print(par.est[paste0(names.item[[factor.no]][[Referent]],"_T", r, "~1")])
#$ print(par.est[paste0(names.lv[factor.no],"_T", r, "=~", names.item[[factor.no]][[Referent]],"_T", r)])
#$ print(par.est[paste0(names.item[[factor.no]][[Referent]],"_T", a, "~1")])
#$ print(par.est[paste0(names.lv[factor.no],"_T", a, "=~", names.item[[factor.no]][[Referent]],"_T", a)])
              } else { ## (Arg != no.markers[factor.no])
                boot.dif.tx[,comp] <- (bootcoef[, paste0(names.item[[factor.no]][[Arg]],"_T", r, "~1")]-
                                       bootcoef[, paste0(names.lv[factor.no],"_T", r, "=~", names.item[[factor.no]][[Arg]],"_T", r)]*
                                       bootcoef[, paste0(names.item[[factor.no]][[Referent]],"_T", r, "~1")]/
                                       bootcoef[, paste0(names.lv[factor.no],"_T", r, "=~", names.item[[factor.no]][[Referent]],"_T", r)]) -
                                      (bootcoef[, paste0(names.item[[factor.no]][[Arg]],"_T", a, "~1")]-
                                       bootcoef[, paste0(names.lv[factor.no],"_T", a, "=~", names.item[[factor.no]][[Arg]],"_T", a)]*
                                       bootcoef[, paste0(names.item[[factor.no]][[Referent]],"_T", a, "~1")]/
                                       bootcoef[, paste0(names.lv[factor.no],"_T", a, "=~", names.item[[factor.no]][[Referent]],"_T", a)])
                samp.dif.tx[r,a] <- (par.est[paste0(names.item[[factor.no]][[Arg]],"_T", r, "~1")]-
                                     par.est[paste0(names.lv[factor.no],"_T", r, "=~", names.item[[factor.no]][[Arg]],"_T", r)]*
                                     par.est[paste0(names.item[[factor.no]][[Referent]],"_T", r, "~1")]/
                                     par.est[paste0(names.lv[factor.no],"_T", r, "=~", names.item[[factor.no]][[Referent]],"_T", r)]) -
                                    (par.est[paste0(names.item[[factor.no]][[Arg]],"_T", a, "~1")]-
                                     par.est[paste0(names.lv[factor.no],"_T", a, "=~", names.item[[factor.no]][[Arg]],"_T", a)]*
                                     par.est[paste0(names.item[[factor.no]][[Referent]],"_T", a, "~1")]/
                                     par.est[paste0(names.lv[factor.no],"_T", a, "=~", names.item[[factor.no]][[Referent]],"_T", a)])
#$ print(paste0("Referent = ", Referent, " & Arg = ", Arg))
#$ print(par.est[paste0(names.item[[factor.no]][[Arg]],"_T", a, "~1")])
#$ print(par.est[paste0(names.lv[factor.no],"_T", a, "=~", names.item[[factor.no]][[Arg]],"_T", a)])
#$ print(par.est[paste0(names.item[[factor.no]][[Referent]],"_T", a, "~1")])
#$ print(par.est[paste0(names.lv[factor.no],"_T", a, "=~", names.item[[factor.no]][[Referent]],"_T", a)])
              }  ## end if nArg
            }  ## end if Referent
          }  ## end loop a
        }  ## end loop r

#$ print(samp.dif.tx)


        ## == Calculate Percentile Probability == ##

        F1.comp.pp <- matrix(0, no.group, no.group)  ## Percentile Probability
        pno_nipair <- 0  ## Set number of non-invariant pairs to zero

        comp = 0
        for (r in 1:(no.group-1)) {  ## r is the referent group
          for (a in (r+1):no.group) {  ## a is the argument
            comp = comp + 1
            if (quantile(boot.dif.tx[, comp], probs = 0.5, na.rm = TRUE) > 0) {
              F1.comp.pp[r,a] = 2*(sum(boot.dif.tx[, comp] < 0, na.rm=TRUE)/bootno)
            } else {
              F1.comp.pp[r,a] = 2*(sum(boot.dif.tx[, comp] > 0, na.rm=TRUE)/bootno)
            }  ## end if
          }  ## end loop a
        }  ## end loop r


        ## == Percentile Confidence Intervals of Pairwise Comparisons == ##

        PCI <- matrix(1:(no.dif*10), nrow = no.dif)
        colnames(PCI) <- c("Group 1","Group 2","0.5%","2.5%","5%","Estimate","95%","97.5%","99.5%","p-value")

        comp = 0
        for (r in 1:(no.group-1)) {
          for (a in (r+1):no.group) {
            comp = comp + 1
            PCI[comp, 1] <- r
            PCI[comp, 2] <- a
            PCI[comp, 3] <- quantile(boot.dif.tx[, comp],c(0.005), na.rm = TRUE)
            PCI[comp, 4] <- quantile(boot.dif.tx[, comp],c(0.025), na.rm = TRUE)
            PCI[comp, 5] <- quantile(boot.dif.tx[, comp],c(0.05), na.rm = TRUE)
            PCI[comp, 6] <- samp.dif.tx[r,a]
            PCI[comp, 7] <- quantile(boot.dif.tx[, comp],c(0.95), na.rm = TRUE)
            PCI[comp, 8] <- quantile(boot.dif.tx[, comp],c(0.975), na.rm = TRUE)
            PCI[comp, 9] <- quantile(boot.dif.tx[, comp],c(0.995), na.rm = TRUE)
            PCI[comp,10] <- F1.comp.pp[r,a]
            if (FL[r,Arg] != FL[a,Arg]) {
              if (pno_nipair == 0)	{
                pnipair <- c(r,a)
                pno_nipair <- pno_nipair + 1
              } else {
                pnipair <- c(pnipair,r,a)
                pno_nipair <- pno_nipair + 1
              }  ## end if pno_nipair
            } else {
              if (F1.comp.pp[r,a] < alpha) {
                if (pno_nipair == 0)	{
                  pnipair <- c(r,a)
                  pno_nipair <- pno_nipair + 1
                } else {
                  pnipair <- c(pnipair,r,a)
                  pno_nipair <- pno_nipair + 1
                }  ## end if pno_nipair
              }  ## end if F1.comp.pp
            }  ## end if FL[, Arg]
          }  ## end loop a
        }  ## end loop r


#$      cat("\n")
#$      if (Referent > Arg) {
#$        cat(paste0("Referent Item: ", Referent, "  Argument Item: ", Arg, "\n"))
#$      } else {
#$        cat(paste0("Referent Item: ", Referent, "  Argument Item: ", (Arg+1), "\n"))
#$      }  ## end if Referent
#$      cat("\n")
#$      print(round(PCI[], digits=4))

        count.tx <- count.tx + 1

        ## == Run List and Delete == ##
#$      cat(rep("\n", 2), "Sets of Groups with Invariant Intercepts", "\n")
        if (pno_nipair == no.dif) {

          txZ <<- matrix(" ", 1, (no.group+2))
          txX <<- matrix(1:(no.group+2), 1, (no.group+2))
          txX[1, 1:no.group] <- 0
          for (temp.txZ in 1: no.group) {
            txZ[1,temp.txZ] <- paste0("F", factor.no, "I", EP)
            EP <<- EP + 1
          }
          txX[1, (no.group+1)] <- Referent
          txZ[1, (no.group+1)] <- Referent
###          if (Referent > Arg) {
            txX[1, (no.group+2)] <- Arg
            txZ[1, (no.group+2)] <- Arg
###          } else {
###            txX[1, (no.group+2)] <- Arg + 1
###            txZ[1, (no.group+2)] <- Arg + 1
###          }  ## end if Referent
#$        print(txX, quote=FALSE)
          txY <<- rbind(txY,txX)
#$        print(txZ, quote=FALSE)
          txYY <<- rbind(txYY,txZ)

        } else if (pno_nipair > 0) {
#          txX <<- matrix(" ", 1, (no.group+2))
#          txZ <<- matrix(" ", 1, (no.group+2))

          listanddelete.tx(factor.no, no.group, no_nipair=pno_nipair, nipair=pnipair, Referent, Arg, txY, txYY, EP) ## Run list and delete tx

#$        print(txX, quote=FALSE)
#$        print(txZ, quote=FALSE)

        } else {  # (pno_nipair = 0)
          txX <<- matrix(1:(no.group+2), 1, (no.group+2))
          txZ <<- matrix(" ", 1, (no.group+2))
          txZ <- matrix(paste0("F", factor.no, "I", EP), 1, (no.group+2))
          EP <<- EP + 1
          txX[1, (no.group+1)] <- Referent
          txZ[1, (no.group+1)] <- Referent
###          if (Referent > Arg) {
            txX[1, (no.group+2)] <- Arg
            txZ[1, (no.group+2)] <- Arg
###          } else {
###            txX[1, (no.group+2)] <- Arg + 1
###            txZ[1, (no.group+2)] <- Arg + 1
###          }  ## end if Referent
#$        print(txX, quote=FALSE)
          txY <<- rbind(txY,txX)
#$        print(txZ, quote=FALSE)
          txYY <<- rbind(txYY,txZ)

        }  ## end if pno_nipair
      } ## end loop Arg
    } ## end loop Referent

    if (length(unique(txY[, no.group+1])) == 2 & no.items[factor.no] == 2) { # Only 2 groups with more than 2 items
      txX <<- matrix(1:(no.group+2), 1, (no.group+2))
      txZ <<- matrix(" ", 1, (no.group+2))
      txYa <<- txY[nrow(txY), no.group+1]
      txYb <<- txY[nrow(txY), no.group+2]
      txX[1, no.group + 1] <- txYb
      txX[1, no.group + 2] <- txYa
      txY <<- rbind(txY,txX)

      txZ[1, no.group + 1] <- txYb
      txZ[1, no.group + 2] <- txYa
      for (temp.txYY in 1: no.group) {
        txZ[1, temp.txYY] <- paste0("F", factor.no, "I", EP)
        EP <<- EP + 1
      }
      txYY <<- rbind(txYY,txZ)
    }

  txY <<- txY[-c(1),]
#$  cat("\n")
#$  print(txY, quote = FALSE)

  txYY <<- txYY[-c(1),]
#$  cat("\n")
#$  print(txYY, quote = FALSE)

  if (is.matrix(txY) == FALSE) { txY <<- matrix(txY, nrow=1) }
  if (is.matrix(txYY) == FALSE) { txYY <<- matrix(txYY, nrow=1) }

    NIcombine <- table(txY[, no.group+1], txY[, no.group+2])
#    NI.level <- unique(txY[, no.group+1])
#    NIcombine <- NIcombine[NI.level, NI.level]


    tempR <- unique(txYY[, no.group+1])

    for (aR in 1:length(tempR)) {
      R <- tempR[aR]
      for (r in 1:(no.items[factor.no])) {
        if (r == 1) {
          assign(paste0("Rset", r), c(length(which(txYY[,no.group+1] < R))+1:(NIcombine[R, r])))
        } else if (r == R) {
          assign(paste0("Rset", r), 0)
        } else {
          assign(paste0("Rset", r), c(length(which(txYY[,no.group+1] < R))+sum(NIcombine[R, 1:(r-1)])+1:(NIcombine[R, r])))
        }
        if (R == 1) { assign("Rset1", 0) }
      }  ## end loop r

      ## expand.grid(Rset1, Rset2, Rset3, Rset4, Rset5 ...)

      temp.row <- eval(parse(text = paste0("expand.grid(", paste0(c("Rset"), 1:no.items[factor.no], collapse=","),")")))
      if (aR == 1) {
        Set.row <- temp.row
      } else {
        Set.row <- rbind(Set.row, temp.row)
      } ## end if R
    }  ## end loop aR

    Model.load <- array(" ", dim = c(no.items[factor.no], no.group, nrow(Set.row)))  # all possible invariance models

    for (Xset in 1: nrow(Set.row)) {
      a <- matrix(" ", no.items[factor.no], no.group)
      for (r in 1:no.items[factor.no]) {
        if (Set.row[Xset, r] == 0) {
          a[r, ] <- 0
        } else {
          a[r, ] <- txYY[Set.row[Xset, r], 1:no.group]
        } ## end if Set.row
###@@@        if (txYY[r, no.group+1] > 1) { break }
      }  ## end loop r
      Model.load[, , Xset] <- a

    }  ## end loop Xset

#$ print(Model.load)

##########################

  ## == Model.summary == ##
  # Column 1 - number of estimated intercepts
  # Column 2 - number of measurement invariance items
  # Column 3 - p-value of Model chi-square
  # Column 4 - CFI
  # Column 5 - RMSEA
  # Column 6 - SRMR

    Model.summary <- matrix(0, nrow(Set.row), 7)
    for (R in 1:nrow(Set.row)) {
      Model.summary[R, 1] <- length(unique(c(Model.load[,,R]))) - 1
      no.MI <- no.items[factor.no] - 1
      for (Ra in 1: no.items[factor.no]) {
        if (length(unique(c(Model.load[Ra,,R]))) > 1) {no.MI <- no.MI - 1}
      }  ## end loop Ra
      Model.summary[R, 2] <- no.MI
    }  ## end loop R


    ## == Recommended Models (Rec.Model) - Largest no. of MI items, then smallest number of estimated loadings == ##
    Rec.Model <- which(Model.summary[,2] == max(Model.summary[,2]))
    Rec.Model <- Rec.Model[which(Model.summary[Rec.Model,1] == min(Model.summary[Rec.Model,1]))]


    ## == Find recommended model if no.items > 2 == ##
    if (no.items[factor.no] > 2) {

      Rec.Model.load <- matrix(nrow = length(Rec.Model), ncol = length(lavaan::coef(PMI.Model.fit)))
      MODFIT <- matrix(0, nrow = length(Rec.Model), ncol = 9)
      colnames(MODFIT) <-
        c("  chisq.scaled", "  df.scaled", "  pvalue.scaled", "  rmsea.scaled", "  cfi.scaled", "  tli.scaled", "  srmr_bentler", "  aic", "  bic")
      Rec.Model.FL <- array(0, dim = c(no.group, no.items[factor.no], length(Rec.Model)))
      colnames(Rec.Model.FL) <- paste0(" Item ", 1:no.items[factor.no])
      rownames(Rec.Model.FL) <- c(paste0("Time ", 1:no.group))

      Rec.Model.TX <- array(0, dim = c(no.group, no.items[factor.no], length(Rec.Model)))
      colnames(Rec.Model.TX) <- paste0(" Item ", 1:no.items[factor.no])
      rownames(Rec.Model.TX) <- c(paste0("Time ", 1:no.group))

      PSI.XX <- array(" ", length(Rec.Model))
      for (Model.R in 1: length(Rec.Model)) {

        ## == Specify the factor loadings for partial scalar invariance models (PSI.X) == ##

        PSI.marker <- which(Model.load[,1,Rec.Model[Model.R]] == 0)
        if (length(PSI.marker) > 1) {PSI.marker <<- PSI.marker[1]}
        PSI.FL <- matrix(1:no.k, nrow = no.group)  # factor loading matrix
        PSI.load <- matrix(" ", no.items[factor.no], no.group)  # factor loading matrix

        for (FL.item.g in 1:no.group) {
          for (FL.item in 1:no.items[factor.no]) {
            PSI.FL[FL.item.g, FL.item] <- FL.PMI[FL.item.g,FL.item]/FL.PMI[FL.item.g,PSI.marker]
          }  ## end loop FL.item
        }  ## end loop FL.item.g

        FL.count <- 1
        for (FL.item in 1:no.items[factor.no]) {
          if (PSI.FL[FL.item.g, FL.item] == 1) {
            PSI.load[FL.item,] <- 1
          } else {
            temp.fl <- unique(PSI.FL[,FL.item])
            for (temp.x in 1:length(temp.fl)) {
              for (FL.item.g in 1:no.group) {
                if (PSI.FL[FL.item.g,FL.item] == temp.fl[temp.x]) {
                  PSI.load[FL.item,FL.item.g] <- paste0("F",factor.no,"L",FL.count)
                }  ## end if
              }  ## end loop FL.item.g
              FL.count <- FL.count + 1
            } ## end loop temp.x
          } ## end if PSI.FL
        }  ## end loop FL.item


        ## == Specify the partial metric invariance models (PMI.X) == ##
        PMI.X.1 <- paste0("PSI.X <- '", "\n")
        for (j in 1:(no.waves)) {
          PMI.X.1 <- paste0(PMI.X.1, "  ", names.lv[factor.no], "_T", j," =~ ", PSI.load[1, j], "*", names.item[[factor.no]][[1]], "_T",j," + ")
          for (i in 2: (no.items[factor.no]-1)) {
            PMI.X.1 <- paste0(PMI.X.1, PSI.load[i, j], "*", names.item[[factor.no]][[i]], "_T",j," + ")
          } ## end loop i
          PMI.X.1 <- paste0(PMI.X.1, PSI.load[(i+1), j], "*", names.item[[factor.no]][[i+1]], "_T",j, "\n")

        } # end for j

        ## -- Create Residual Covariance -- ##
        PMI.X.1 <- paste0(PMI.X.1, "\n")
        for (i in 1:no.items[factor.no]) {
          for (j in 1:(no.waves - 1)) {
            for (k in (j+1):no.waves) {
              PMI.X.1 <- paste0(PMI.X.1, "  ", names.item[[factor.no]][[i]], "_T", j, " ~~ ", names.item[[factor.no]][[i]], "_T", k, "\n")
            } # end for k
          } # end for j
        } # end for i


        ## == Specify the partial scalar invariance models (PSI.X) == ##

        final.intercept <- unique(as.vector(Model.load[,, Rec.Model[Model.R]]))
        for (j in 1:no.waves) {
          for (i in 1:no.items[factor.no]) {
            if (Model.load[i,j,Rec.Model[Model.R]] != 0) { Model.load[i,j, Rec.Model[Model.R]] <- paste0("F", factor.no, "I", which(Model.load[i,j,Rec.Model[Model.R]] == final.intercept)) }
          } # end for i
        } # end for j

        PSI.XX[Model.R] <- paste0(PMI.X.1, "\n")
        for (j in 1:(no.waves)) {
          for (i in 1:no.items[factor.no]) {
            PSI.XX[Model.R] <- paste0(PSI.XX[Model.R], "  ", names.item[[factor.no]][[i]], "_T", j, " ~ ", Model.load[i, j, Rec.Model[Model.R]], "*1", "\n")
          } ## end for i
        } # end for j

        PSI.XX[Model.R] <- paste0(PSI.XX[Model.R], "'", "\n")
        eval(parse(text = PSI.XX[Model.R]))

#$ cat(PSI.XX[Model.R])

        ## == Run PSI.X == ##
        if (Cluster == "NULL") {
          eval(parse(text=   "PSI.X.fit <- lavaan::sem(PSI.X,
                           data.source,
                           missing = 'fiml',
                           auto.fix.first = FALSE,
                           marker.int.zero = TRUE,
                           meanstructure = T,
                           information = 'observed',
                           estimator = 'MLR')"))
        } else {
          eval(parse(text=   "PSI.X.fit <- lavaan::sem(PSI.X,
                           data.source,
                           missing = 'fiml',
                           auto.fix.first = FALSE,
                           marker.int.zero = TRUE,
                           meanstructure = T,
                           information = 'observed',
                           cluster = Cluster,
                           estimator = 'MLR')"))
        }

        ## == Request summary outputs == ##
#$      eval(parse(text = "print(lavaan::summary(PSI.X.fit, fit.measure = T, standardized = T, rsq = T))"))


        ## == Save fit indices == ##
        if (correct.cfi == TRUE) {
          fitmeasures <- LGcorrected.cfi(PSI.object=PSI.X.fit, data=data.source, cluster=Cluster)
          XX <- fitmeasures[c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr_bentler", "aic", "bic")]
        } else {
          XX <- lavaan::fitMeasures(PSI.X.fit,
                                  c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr_bentler", "aic", "bic"))
        } # end (if correct.cfi)

        MODFIT[Model.R,] <- format(round((XX), digits = 4), nsmall = 4, scientific = FALSE)

        ## == Save Intercepts == ##
        temp.tx <- lavaan::lavInspect(PSI.X.fit, what='est')$nu
        for (G in 1: no.group) {
          Rec.Model.TX[G, , Model.R] <- temp.tx[(no.items[factor.no]*(G-1)+1):(no.items[factor.no]*G)]
        }  ## end loop G
      } ## end loop Model.R


      ## == Print model fit indices == ##
      rownames(MODFIT) <- c(paste0("Model ", 1: length(Rec.Model)))
#$    cat(rep("\n", 3))
#$    cat("####    MODEL FIT INDICES   ####", rep("\n", 2))
#$    print(MODFIT, quote=FALSE, right=TRUE)

      ## == Print factor loadings == ##
#$    for (Model.R in 1: length(Rec.Model)) {
#$      cat(rep("\n", 2), "Factor Loadings of Model", Model.R, rep("\n", 2))
#$      print(round(Rec.Model.FL[,,Model.R], digits=4))
#$    }  ## end loop Model.R


      ## == Save Recommended Model (Recommend.Model) -- First model with best model fit (R.Model) == ##
      class(MODFIT) <- "numeric"
      if (BMSC == "Chisquare") {
        Model.R <- which.min(MODFIT[,1])
      } else if (BMSC == "CFI") {
        Model.R <- which.max(MODFIT[,5])
      } else if (BMSC == "RMSEA") {
        Model.R <- which.min(MODFIT[,4])
      } else if (BMSC == "SRMR") {
        Model.R <- which.min(MODFIT[,7])
      }
      R.Model <- Model.R[1] # Select the first model if 2 or more models have the same fit
      if (factor.no == 1) {Recommend.Model <- matrix("PSI.Model.R <- '", 1)}  ## Reset Recommend.Model

      Recommend.Model <- rbind(Recommend.Model, PSI.XX[R.Model])

    }  ## end if no.items[factor.no] > 2


    if (no.items[factor.no] < 3) {   ## Begin Models with 2 items ##
        Rec.Model <- Rec.Model[1] # Select the first model if 2 or more models have the same fit
        PSI.marker <- which(Model.load[,1,Rec.Model] == 0)
        if (length(PSI.marker) > 1) {PSI.marker <<- PSI.marker[1]}
        PSI.FL <- matrix(1:no.k, nrow = no.group)  # factor loading matrix
        PSI.load <- matrix(" ", no.items[factor.no], no.group)  # factor loading matrix

        for (FL.item.g in 1:no.group) {
          for (FL.item in 1:no.items[factor.no]) {
            PSI.FL[FL.item.g, FL.item] <- FL.PMI[FL.item.g,FL.item]/FL.PMI[FL.item.g,PSI.marker]
          }  ## end loop FL.item
        }  ## end loop FL.item.g

        FL.count <- 1
        for (FL.item in 1:no.items[factor.no]) {
          if (PSI.FL[FL.item.g, FL.item] == 1) {
            PSI.load[FL.item,] <- 1
          } else {
            temp.fl <- unique(PSI.FL[,FL.item])
            for (temp.x in 1:length(temp.fl)) {
              for (FL.item.g in 1:no.group) {
                if (PSI.FL[FL.item.g,FL.item] == temp.fl[temp.x]) {
                  PSI.load[FL.item,FL.item.g] <- paste0("F",factor.no,"L",FL.count)
                }  ## end if
              }  ## end loop FL.item.g
              FL.count <- FL.count + 1
            } ## end loop temp.x
          } ## end if PSI.FL
        }  ## end loop FL.item

      R.Model <- Model.load[,,Rec.Model]  # Recommended Model

      Recommend.Model <- rbind(Recommend.Model, "\n")


      ## == Save Recommended Model (Recommend.Model) == ##
        if (factor.no == 1) {Recommend.Model <- matrix("PSI.Model.R <- '", 1)}

        for (j in 1:(no.waves)) {
          Recommend.Model <- rbind(Recommend.Model, paste0("  ", names.lv[factor.no], "_T", j," =~ ", PSI.load[1, j], "*", names.item[[factor.no]][[1]], "_T",j," + "))
          Recommend.Model <- rbind(Recommend.Model, paste0(PSI.load[2, j], "*", names.item[[factor.no]][[2]], "_T", j, "\n"))
        } # end for j

        ## -- Create Residual Covariance -- ##
        Recommend.Model <- rbind(Recommend.Model, "\n")
        for (i in 1:no.items[factor.no]) {
          for (j in 1:(no.waves - 1)) {
            for (k in (j+1):no.waves) {
              Recommend.Model <- rbind(Recommend.Model, paste0("  ", names.item[[factor.no]][[i]], "_T", j, " ~~ ", names.item[[factor.no]][[i]], "_T", k, "\n"))
            } # end for k
          } # end for j
        } # end for i

        ## == Specify the partial scalar invariance models (PSI.X) == ##

        final.intercept <- unique(as.vector(Model.load[,, Rec.Model[Model.R]]))
        for (j in 1:no.waves) {
          for (i in 1:no.items[factor.no]) {
            if (Model.load[i,j,Rec.Model[Model.R]] != 0) { Model.load[i,j, Rec.Model[Model.R]] <- paste0("F", factor.no, "I", which(Model.load[i,j,Rec.Model[Model.R]] == final.intercept)) }
          } # end for i
        } # end for j

        Recommend.Model <- rbind(Recommend.Model, "\n")
        for (j in 1:(no.waves)) {
          for (i in 1:no.items[factor.no]) {
            Recommend.Model <- rbind(Recommend.Model, paste0("  ", names.item[[factor.no]][[i]], "_T", j, " ~ ", Model.load[i, j, Rec.Model[Model.R]], "*1", "\n"))
          } # end for i
        } # end for j

#        Recommend.Model <- rbind(Recommend.Model, "'", "\n")


    }  ## End (Models with 2 items)

  } ## End (loop factor.no)



  Recommend.Model[1] <- sub("''","'", Recommend.Model[1])
  for (i in 1:no.factor) {
    Recommend.Model[i+1] <-  sub("PSI.X <- '", " ", Recommend.Model[i+1])
    Recommend.Model[i+1] <- gsub("\n'", " ", Recommend.Model[i+1])
  } # end loop i
  Recommend.Model <- rbind(Recommend.Model, "  '", "\n")  ## last line of Recommended Model


  cat(rep("\n", 2), "## =====  Recommended Model  ===== ##", rep("\n", 2))
  for (i in 1: nrow(Recommend.Model)) {
    if (i == 1) {
      cat(Recommend.Model[i], "\n")
    } else {
      cat(Recommend.Model[i])
    }
  } ## end loop i

  cat(rep("\n", 2))
  cat("## ===== lavaan Outputs of PSI.Model.R ===== ##")
  cat(rep("\n", 2))

  ## == Start print recommended model to PSI.txt (Partial Metric Invariance Model) == ##
  sink('PSI.txt')

  cat("## =====  Recommended Model  ===== ##")
  cat("\n")
    for (i in 1: nrow(Recommend.Model)) {
      if (i == 1) {
        cat(Recommend.Model[i], "\n")
      } else {
        cat(Recommend.Model[i])
      }
    } ## end loop i


  cat(rep("\n",2), "## Run PSI.Model.R")
  cat(rep("\n",2), paste0("  PSI.Model.fit <- lavaan::sem(PSI.Model.R,"))
  cat("\n", paste0("    ", arg2_char, ","))
  cat("\n", "    missing = 'fiml',")
  cat("\n", "    auto.fix.first = FALSE,")
  cat("\n", "    marker.int.zero = TRUE,")
  cat("\n", "    meanstructure = T,")
  cat("\n", "    information = 'observed',")
  if (Cluster != "NULL") { cat("\n", paste0("    cluster = ", arg4_char,",")) }
  cat("\n", "    estimator = 'MLR')")

  cat(rep("\n", 3), "## Request summary outputs")
  cat(rep("\n", 2), paste0("  print(lavaan::summary(PSI.Model.fit, fit.measure = T, standardized = T, rsq = T))"), "\n")

  sink()  ## Stop writing to file

#  source('PSI.txt') ## Run the script
  PSI.Model.R <<- eval(parse(text=Recommend.Model))

  ## == Run PSI.Model.R == ##
  if (Cluster == "NULL") {
    eval(parse(text=   "PSI.Model.fit <- lavaan::sem(PSI.Model.R,
                       data.source,
                       missing = 'fiml',
                       auto.fix.first = FALSE,
                       marker.int.zero = TRUE,
                       meanstructure = T,
                       information = 'observed',
                       estimator = 'MLR')"))
  } else {
    eval(parse(text=   "PSI.Model.fit <- lavaan::sem(PSI.Model.R,
                       data.source,
                       missing = 'fiml',
                       auto.fix.first = FALSE,
                       marker.int.zero = TRUE,
                       meanstructure = T,
                       information = 'observed',
                       cluster = Cluster,
                       estimator = 'MLR')"))
  }

  print(lavaan::summary(PSI.Model.fit, fit.measure = T, standardized = T, rsq = T))

  cat(rep("\n", 2), "## Compare fit indices across models ##", "\n")
  DMODFIT <- matrix(0, nrow = 3, ncol = 9)
  colnames(DMODFIT) <-
        c("  chisq.scaled", "  df.scaled", "  pvalue.scaled", "  rmsea.scaled", "  cfi.scaled", "  tli.scaled", "  srmr", "  aic", "  bic")
  rownames(DMODFIT) <- c("PMI.Model.fit", "PSI.Model.fit", "PSI.Model.fit - PMI.Model.fit")
  XX1 <- lavaan::fitMeasures(PMI.Model.fit, fit.measures=
      c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr", "aic", "bic"))
  DMODFIT[1,] <- XX1

  if (correct.cfi == TRUE) {
    fitmeasures <- LGcorrected.cfi(PSI.object=PSI.Model.fit, data=data.source, cluster=Cluster)
    XX2 <- fitmeasures[c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr_bentler", "aic", "bic")]
  } else {
    XX2 <- lavaan::fitMeasures(PSI.Model.fit, fit.measures=
                             c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr", "aic", "bic"))
  } # end (if correct.cfi)

  DMODFIT[2,] <- XX2

  DMODFIT[3,1] <- DMODFIT[2,1] -DMODFIT[1,1]
  DMODFIT[3,2] <- DMODFIT[2,2] -DMODFIT[1,2]
  DMODFIT[3,3] <- pchisq(q=DMODFIT[3,1], df=DMODFIT[3,2], lower.tail=FALSE)
  DMODFIT[3,4] <- DMODFIT[2,4] -DMODFIT[1,4]
  DMODFIT[3,5] <- DMODFIT[2,5] -DMODFIT[1,5]
  DMODFIT[3,6] <- DMODFIT[2,6] -DMODFIT[1,6]
  DMODFIT[3,7] <- DMODFIT[2,7] -DMODFIT[1,7]
  DMODFIT[3,8] <- DMODFIT[2,8] -DMODFIT[1,8]
  DMODFIT[3,9] <- DMODFIT[2,9] -DMODFIT[1,9]
  print(round(DMODFIT, digits=4))


  ## == Print Factor Loadings == ##
  Final.Model.FL <- matrix(0, sum(no.items), no.waves)
  rownames(Final.Model.FL) <- paste0(" Item ", names.ind[1:sum(no.items)])
  colnames(Final.Model.FL) <- c(paste0("Time ", 1:no.group))

  Final.par.est <- lavaan::parameterEstimates(PSI.Model.fit)
  ext <- c(which(Final.par.est[,"op"] == "=~"))
  aa <- Final.par.est[ext, "est"]

  ## Extract factor loadings ##
  n <- 1
  for (j in 1: no.waves) {
    Final.Model.FL[1:no.items[1], j] <- aa[(no.items[1]*(j-1)+1):(no.items[1]*j)]
  }
  n <- n + no.items[1]*no.waves

  if(no.factor > 1){
    for (i in 2:no.factor) {
      for (j in 1: no.waves) {
        Final.Model.FL[(sum(no.items[1:(i-1)])+1):sum(no.items[1:i]), j] <- aa[(n+no.items[i]*(j-1)):(n+(no.items[i]*j-1))]
      }
      n <- n + no.items[i]*no.waves
    }
  }

  cat("\n")
  cat("## === Factor Loadings in Final Model === ##")
  cat("\n")
  print(round(Final.Model.FL, digits=4))
  ## == End (Print Factor Loadings) == ##



  ## == Print Intercepts == ##
  Final.Model.TX <- matrix(0, sum(no.items), no.group)
  rownames(Final.Model.TX) <- paste0(" Item ", names.ind[1:sum(no.items)])
  colnames(Final.Model.TX) <- c(paste0("Time ", 1:no.group))

  temp.tx <- lavaan::lavInspect(PSI.Model.fit, what='est')$nu
  for (G in 1: no.group) {
    Final.Model.TX[1:no.items[1], G] <- temp.tx[(no.items[1]*(G-1)+1):(no.items[1]*G)]
  }
  if (no.factor>1) {
    for (i in 2: no.factor) {
      for (G in 1: no.group) {
        Final.Model.TX[(sum(no.items[1:(i-1)])+1):(sum(no.items[1:(i-1)])+no.items[i]),G] <-
         temp.tx[(sum(no.items[1:(i-1)])*no.waves+no.items[i]*(G-1)+1):(sum(no.items[1:(i-1)])*no.waves+no.items[i]*G)]
      }
    }
  }
  cat("\n")
  cat("## === Intercepts in Final Model === ##")
  cat("\n")
  print(round(Final.Model.TX, digits=4))

  ## == End (Print Intercepts) == ##


  ## == Print Latent Means == ##
  Final.Model.LM <- matrix(0, no.factor, no.group)
  rownames(Final.Model.LM) <- paste0(" ", names.lv)
  colnames(Final.Model.LM) <- c(paste0("Time ", 1:no.group))

  temp.LM <- lavaan::lavInspect(PSI.Model.fit, what='est')$alpha

  for (G in 1: factor.no) { Final.Model.LM[G,] <- temp.LM[(no.group*(G-1)+1):(no.group*G)]}
  cat("\n")
  cat("## === Latent Means in Final Model === ##")
  cat("\n")
  print(round(Final.Model.LM, digits=4))
  ## == End (Print Latent Means) == ##


  ## ========== Compare Latent Means ========== ##

  ## Extraxt latent means ##
  par.est <- lavaan::coef(PSI.Model.fit)  # sample parameters
  temp <- lavaan::parameterEstimates(PSI.Model.fit, remove.nonfree=TRUE)
  simvcov <- lavInspect(PSI.Model.fit, what="vcov")

#  ext <- c(which(temp$op %in% "~1" & temp$lhs %in% names.lv))

  parsed <- temp[temp$op == "=~",]
  names.lm <- unique(parsed[,"lhs"])

  ext <- c(which(temp$op %in% "~1" & temp$lhs %in% names.lm))


  par.est <- par.est[ext]
  simvcov <- simvcov[ext,ext]
  no.lm.g <- length(names.lv)  # number of estimated latent means per group #

  cat(rep("\n", 3), "## ======= COMPARISON OF LATENT MEANS ======= ##", rep("\n", 2))  ## print heading

  if (TYPE == "MonteCarlo") {

    mcmc <- MASS::mvrnorm(n=2000000, mu=par.est, Sigma=simvcov, tol = 1e-6)  # Run 2,000,000 simulations
    bootcoef <- mcmc
    bootno <- nrow(mcmc)  # No. of successful simulated samples
    cat(paste0("Number of Successful Simulated Samples = ", bootno, "\n"))

  } else { # TYPE=Bootstrap

    ## == Simplified bootstrapping model == ##
    PSI.boot <- lavaan::cfa(PSI.Model.R, data = data.source,
                  meanstructure = TRUE,
                  auto.fix.first = FALSE,
                  marker.int.zero = TRUE,
                  ordered = FALSE,
                  missing = "fiml", se = "none", test = "none", check.start = FALSE, check.post = FALSE, check.gradient = FALSE)

    ## == Bootstrapping == ##
    bootcoef <- lavaan::bootstrapLavaan(PSI.boot, R = b.no, ncpus = max(1L, parallel::detectCores() - 1L), parallel="snow")

    # Remove error and nonadmissible bootstrap samples #
    B.na <- attr(bootcoef,"nonadmissible")
    B.er <- attr(bootcoef,"error.idx")
    B.de <- c(B.na,B.er)
    if (length(B.de) != 0) {
      B.re <- bootcoef[-c(B.de),]
      bootcoef <- B.re
    }
    bootno <- nrow(bootcoef)  # number of successful bootstrap samples
    cat(paste0("Number of Successful Bootstrapped Samples = ", bootno, "\n"))
    bootcoef <- bootcoef[,ext]  ## Extract the bootstrapped intercepts

  } ## end MonteCarlo or Bootstrap


  ## == Start the factor.no loop for CompareMeans == ##
  for (factor.no in 1: no.factor) {

    boot.dif.lm <- matrix(0, nrow(bootcoef), no.dif)  # Create bootstrap difference matrix
    samp.dif.lm <- matrix(0, no.group, no.group) # Create sample difference matrix

    ## == Calculate bootstrap difference and sample estimate difference == ##
    comp = 0
    for (r in 1:(no.group-1)) {  ## r is referent group
      for (a in (r+1):no.group) {  ## a is argument group
        comp = comp + 1
        boot.dif.lm[,comp] <- bootcoef[, paste0(names.lv[factor.no], "_T", r, "~1")] - bootcoef[,paste0(names.lv[factor.no], "_T", a, "~1")]
        samp.dif.lm[r,a] <- par.est[paste0(names.lv[factor.no], "_T", a, "~1")] - par.est[paste0(names.lv[factor.no], "_T", r, "~1")]
        samp.dif.lm[a,r] <- par.est[paste0(names.lv[factor.no], "_T", r, "~1")] - par.est[paste0(names.lv[factor.no], "_T", a, "~1")]
      }  ## end loop a
    }  ## end loop r

    colnames(samp.dif.lm) <- c(paste0("Time ", 1:no.group))
    rownames(samp.dif.lm) <- c(paste0("Time ", 1:no.group))


    ## == Calculate Number of Invariant Intercepts == ##
    LM.inv.tau <- matrix(0, no.group, no.group)
    for (r in 1:(no.group-1)) {  ## r is the referent group
      for (a in (r+1):no.group) {  ## a is the argument
        if (factor.no == 1) {
          LM.inv.tau[r,a] <- sum(round(Final.Model.TX[1:no.items[factor.no],r],6) == round(Final.Model.TX[1:no.items[factor.no],a],6)) 
        } else {
          LM.inv.tau[r,a] <- sum(round(Final.Model.TX[(sum(no.items[1:(factor.no-1)])+1):sum(no.items[1:factor.no]),r], 6) == 
                                 round(Final.Model.TX[(sum(no.items[1:(factor.no-1)])+1):sum(no.items[1:factor.no]),a], 6)) 
        }
        LM.inv.tau[a,r] <- LM.inv.tau[r,a]
      } # end for a
    } # end for r 


    ## == Calculate Percentile Probability == ##
    LM.comp.pp <- matrix(0, no.group, no.group)  ## Percentile Probability
    colnames(LM.comp.pp) <- c(paste0("Time ", 1:no.group))
    rownames(LM.comp.pp) <- c(paste0("Time ", 1:no.group))

    comp = 0
    for (r in 1:(no.group-1)) {  ## r is the referent group
      for (a in (r+1):no.group) {  ## a is the argument
        comp = comp + 1
        if (quantile(boot.dif.lm[, comp], probs = 0.5, na.rm = TRUE) > 0) {
          LM.comp.pp[r,a] = 2*(sum(boot.dif.lm[, comp] < 0, na.rm=TRUE)/bootno)
          LM.comp.pp[a,r] = 2*(sum(boot.dif.lm[, comp] < 0, na.rm=TRUE)/bootno)
        } else {
          LM.comp.pp[r,a] = 2*(sum(boot.dif.lm[, comp] > 0, na.rm=TRUE)/bootno)
          LM.comp.pp[a,r] = 2*(sum(boot.dif.lm[, comp] > 0, na.rm=TRUE)/bootno)
        }  ## end if
      }  ## end loop a
    }  ## end loop r

    ## == Percentile Confidence Intervals of Pairwise Comparisons == ##

    PCI <- matrix(1:(no.dif*10), nrow = no.dif)
    colnames(PCI) <- c("Group 1","Group 2","0.5%","2.5%","5%","Estimate","95%","97.5%","99.5%","p-value")

    comp = 0
    for (r in 1:(no.group-1)) {
      for (a in (r+1):no.group) {
        comp = comp + 1
        PCI[comp, 1] <- paste0("Time ", r)
        PCI[comp, 2] <- paste0("Time ", a)
        PCI[comp, 3] <- round(quantile(boot.dif.lm[, comp],c(0.005), na.rm = TRUE), digits=4)
        PCI[comp, 4] <- round(quantile(boot.dif.lm[, comp],c(0.025), na.rm = TRUE), digits=4)
        PCI[comp, 5] <- round(quantile(boot.dif.lm[, comp],c(0.05), na.rm = TRUE), digits=4)
        PCI[comp, 6] <- round(samp.dif.lm[r,a], digits=4)
        PCI[comp, 7] <- round(quantile(boot.dif.lm[, comp],c(0.95), na.rm = TRUE), digits=4)
        PCI[comp, 8] <- round(quantile(boot.dif.lm[, comp],c(0.975), na.rm = TRUE), digits=4)
        PCI[comp, 9] <- round(quantile(boot.dif.lm[, comp],c(0.995), na.rm = TRUE), digits=4)
        PCI[comp,10] <- round(LM.comp.pp[r,a], digits=4)
      }  ## end loop a
    }  ## end loop r

    cat("\n")
#$  print(PCI[],quote=F, nsmall=4, scientific=FALSE)

    # == Combine LM.comp.pp with LM.inv.tau == #
    LM.comp.pp <- matrix(sprintf("%.4f", LM.comp.pp), nrow = nrow(LM.comp.pp))
    LM.inv.tau <- matrix(sprintf("%.0f", LM.inv.tau), nrow = nrow(LM.inv.tau))
    LM.comp.pp <- matrix(paste0(LM.comp.pp, " (", LM.inv.tau, ")  "), nrow = nrow(LM.comp.pp), ncol = ncol(LM.comp.pp))
    colnames(LM.comp.pp) <- c(paste0("Time ", 1:no.group))
    rownames(LM.comp.pp) <- c(paste0("Time ", 1:no.group))
    diag(LM.comp.pp) <- ""

    samp.dif.lm <- matrix(sprintf("%.4f", samp.dif.lm), nrow = nrow(samp.dif.lm))
    diag(samp.dif.lm) <- ""
    colnames(samp.dif.lm) <- c(paste0("Time ", 1:no.group))
    rownames(samp.dif.lm) <- c(paste0("Time ", 1:no.group))


    cat("\n")
    cat("## ===== Latent Variable: ", names.lv[factor.no], " ===== ##")
    cat("\n")
    cat("== Factor Loadings ==")
    cat("\n")
    if (factor.no == 1) {
      print(format(round(Final.Model.FL[1:no.items[factor.no],],4)), quote=F)
    } else {
      print(format(round(Final.Model.FL[(sum(no.items[1:(factor.no-1)])+1):sum(no.items[1:factor.no]),],4)), quote=F)
    }
    cat("\n")
    cat("== Intercepts ==")
    cat("\n")
    if (factor.no == 1) {
      print(format(round(Final.Model.TX[1:no.items[factor.no],],4)), quote=F)
    } else {
      print(format(round(Final.Model.TX[(sum(no.items[1:(factor.no-1)])+1):sum(no.items[1:factor.no]),],4)), quote=F)
    }
    cat("\n")
    cat("== Latent Means ==")
    cat("\n")
    print(Final.Model.LM[factor.no,], quote=F, digits=5, nsmall=4, scientific=FALSE)
    cat("\n")
    cat("== Pairwise Difference in Latent Means ==")
    cat("\n")
    print(format(samp.dif.lm, justify="right"), quote=F)
    cat("\n")
    cat("== p-values for pairwise comparisons ==")
    cat("\n")
    print(format(LM.comp.pp, justify="right"), quote=F)
    cat("\n")
    cat("Note:", "\n")
    cat("Numbers in parentheses are numbers of items with invariant intercepts.", "\n")
    cat("There are ", no.items[factor.no], " items for this factor. Latent means should only be compared with at least ", ceiling((no.items[factor.no]+1)/2), " items with invariant intercepts.", "\n")
    cat("There are ", choose(no.group, 2), "pairwise comparisons.", "\n")
    cat("A critical p-value of", sprintf("%.4f", 0.05/choose(no.group, 2)), "is recommended for an overall Type I error rate of 0.05 with Bonferroni adjustment.", "\n")
    cat("A critical p-value of", sprintf("%.4f", 0.01/choose(no.group, 2)), "is recommended for an overall Type I error rate of 0.01 with Bonferroni adjustment.", "\n")

  } # End (factor.no loop)

  cat(rep("\n",2),"The recommended model PSI.Model.R is saved in the file 'PSI.txt'", "\n")

} ## End (Function LGCompareMeans)

# ==================== Finish Function "LGCompareMeans" ==================== #







# ==================== Create Function "PAIRCompareLoadings" ==================== #
#' Metric Invariance Test Between Two Sources
#'
#' Conduct metric invariance test and identify a partial metric invariance model for data from two non-independent sources.
#'
#' Requires defining the measurement model only once. All indicators should have the suffix _G1 and _G2 (e.g., x1_G1, x1_G2, x2_G1, x2_G2) in the data file to indicate the source from which the indicators were measured.
#'
#' Residuals of indicators are covaried across sources automatically.
#'
#' @param model User-specified measurement model.
#' @param data.source A data frame containing the observed variables used in the model.
#' @param Cluster Cluster variable for nested data. The Monte Carlo simulation method should be used for nested data.
#' @param Bootstrap Number of bootstrap samples, must be between 500 and 10,000. If not specified, 1 million Monte Carlo simulated samples (default) will be used.
#' @param Type1 Overall Type I error rate (default is 0.05) for identifying non-invariant items in the List and Delete method.
#' @param Type1Adj Adjustment of Type I error rate for multiple tests for each construct. Default is Possible False Discovery Rate (PFDR = Type I error rate/No. of freely estimated factor loadings/intercepts across all waves), can also use (BON)ferroni adjustment (Type I error rate/No. of pairwise comparisons), or "NULL".
#' @param BMSC Best Model Selection Criterion used to select among models with the same number of estimated parameters. "ChiSquare" = minimum Chi-Square value, "CFI" = maximum CFI value, "RMSEA" = minimum RMSEA value, and "SRMR" = minimum SRMR value (default).
#' @return Partial metric invariance model in the PMI.txt file.
#' @export
#' @examples
#'
#' ## == Example C - Non-independent Data from Two Sources == ##
#'
#' # Data file is "Example.C"
#'
#' ## Specify the measurement model - Model.C ##
#' Model.C <- '
#'      External =~ x1 + x2 + x3 + x4 + x5 + x6
#'      Internal =~ x7 + x8 + x9 + x10
#' '
#'
#' ## ===== Compare Factor Loadings ===== ##
#' PAIRCompareLoadings(Model.C, Example.C, Type1=0.01, Type1Adj="PFDR", BMSC="SRMR")

#'
PAIRCompareLoadings <- function(model, data.source, Cluster="NULL", Bootstrap=0, Type1 = 0.01, Type1Adj = "PFDR", BMSC="SRMR") {

# model <- Model.C
# data.source <- Example.C
# Cluster="NULL"
# Bootstrap=0
# Type1 = 0.01
# Type1Adj = "PFDR"
# BMSC="SRMR"

  model <<- paste0("\n", "  ", model, "\n")

  Type1Adj <- toupper(Type1Adj)
  match.arg(Type1Adj, c("PFDR","BON", "NULL"))
  BMSC <- toupper(BMSC)
  match.arg(BMSC, c("CHISQUARE","CFI", "RMSEA", "SRMR"))

  arg2_char <<- deparse(substitute(data.source))
  arg4_char <<- deparse(substitute(Cluster))

  no.group <- 2
  parsed <- lavParseModelString(model, as.data.frame = TRUE)
  parsed <- parsed[, c("lhs", "op", "rhs")]
  parsed <- parsed[parsed$op == "=~",]
  names.lv <- unique(parsed[,"lhs"])
  names.ind <- unique(parsed[, "rhs"])
  no.factor <- length(names.lv)
  no.items.g <- length(names.ind) # number of items per group
  no.items <- matrix(1:no.factor, nrow = 1)  # number of items per factor
  for (factor.no in 1:no.factor) { no.items[factor.no] <- sum(parsed[,"lhs"] == names.lv[factor.no] & parsed[,"op"] == "=~") } ## end loop factor.no

  names.item <- list()
  for (i in 1:no.factor) {
    names.item[[i]] <- parsed[which(parsed[,"lhs"] == names.lv[i] & parsed[,"op"] == "=~"), "rhs"]
  }

  ## Check for bootstrap sample number (Bootstrap) ##
  if (Bootstrap !=0) {
    b.no.integer <- Bootstrap == round(Bootstrap)
    if (b.no.integer == "FALSE") stop("Bootstrap sample number must be an integer")
    if (Bootstrap > 10000) stop("Bootstrap sample number greater than 10,000 is not recommended")
    if (Bootstrap < 500) stop("Bootstrap sample number smaller than 500 is not recommended")
    TYPE = "Bootstrap"
    b.no <- Bootstrap
  } else {
    TYPE = "MonteCarlo"
  } # end (Bootstrap != 0)

  for (i in 1: no.factor) { model <- eval(parse(text = (paste0("sub('", names.lv[i],"', '", names.lv[i], "_G1', model, fixed=TRUE)")))) }
  for (i in 1: no.items.g) { model <- eval(parse(text = (paste0("sub('", names.ind[i],"', '", names.ind[i], "_G1', model, fixed=TRUE)")))) }

  ## -- Create Model.Long with all waves -- ##
  Model.Long <- paste0("\n", model, "\n")
  Model.Long <- paste0(Model.Long, gsub("_G1", paste0("_G2"), model), "\n")

  ## == Create Residual Covariance == ##
  for (i in 1:no.items.g) {
    Model.Long <- paste0(Model.Long, "  ", names.ind[i], "_G1 ~~ ", names.ind[i], "_G2", "\n")
  }

  ## Finish creating Model.Long ##


  ## ===== Run Configural Model ===== ##

  if (Cluster == "NULL") {
    Model.Long.config <- lavaan::sem(Model.Long,
                                data.source,
                                missing = 'fiml',
                                marker.int.zero = TRUE,
                                meanstructure = T,
                                estimator = 'MLR')
  } else {
    Model.Long.config <- lavaan::sem(Model.Long, cluster=Cluster,
                                data.source,
                                missing = 'fiml',
                                marker.int.zero = TRUE,
                                meanstructure = T,
                                estimator = 'MLR')
  } # end if cluster

  ## Request summary outputs
  cat(rep("\n", 2), "## ===== Configural Invariance Model ===== ##", "\n")
  print(lavaan::summary(Model.Long.config, fit.measure = T, standardized = T, rsq = T))
  #$  parameterEstimates(Model.Long.config)
  cat("\n")

  Model.Long.config <<- Model.Long.config

  ## ===== End (Run Configural Model) ===== ##


  par.est <- lavaan::coef(Model.Long.config)  # sample parameters

  temp <- lavaan::parameterEstimates(Model.Long.config)
  no.markers <- matrix(1:no.factor, nrow = 1)  # location of marker items
  for (factor.no in 1:no.factor) {
    no.markers[factor.no] <- which(temp[,"lhs"] == paste0(names.lv[factor.no],"_G1") & temp[,"op"] == "=~" & temp[,"est"] == 1)
    if (factor.no > 1) {no.markers[factor.no] <- no.markers[factor.no] - sum(no.items[1:(factor.no-1)])}
  } ## end loop factor.no

  temp <- lavaan::parameterEstimates(Model.Long.config, remove.nonfree=TRUE)
  simvcov <- lavInspect(Model.Long.config, what="vcov")
  par.est <- lavaan::coef(Model.Long.config)  # sample parameters

  ## Extraxt factor loadings ##
  ext <- c(which(temp[,"op"] == "=~"))
  par.est <- par.est[ext]
  simvcov <- simvcov[ext,ext]
  no.par.g <- length(par.est)/2  # number of estimated LX per group #



  cat(rep("\n", 3), "## ======= METRIC INVARIANCE ANALYSIS ======= ##", rep("\n", 2))  ## print heading

  # Calculate alpha for factor[factor.no] #
  ERate <- matrix(1:no.factor, nrow = 1)  # error rate for each factor

  cat("# -- Overall Type I Error Rates for Invariance Tests for Each Factor = ", Type1, " -- #", rep("\n", 2))

  if (Type1Adj == "PFDR") {  ## print adjusted error rate
    cat("# -- Adjusted Type I Error Rates for Invariance Tests (Adjusted by Possible False Discovery Rate) -- #", rep("\n", 2))
    for (factor.no in 1: no.factor) {
      no.adjust <- (no.items[factor.no]-1)*(no.group-1)
      ERate[factor.no] <- Type1 / no.adjust
      cat(paste0("Factor ", factor.no, ": Number of items = ", no.items[factor.no], "; number of free factor loadings = ", no.adjust, "; adjusted error rate = ", sprintf("%.4f", ERate[factor.no]), "\n"))
    } # end (for factor.no)
  } else if (Type1Adj == "BON") {
    cat("# -- Adjusted Type I Error Rates for Invariance Tests (Bonforroni adjustment for number of pairwise comparisons) -- #", rep("\n", 2))
    for (factor.no in 1: no.factor) {
      no.adjust <- (no.items[factor.no]-1)*(no.group*(no.group-1)/2)
      ERate[factor.no] <- Type1 / no.adjust
      cat(paste0("Factor ", factor.no, ": Number of items = ", no.items[factor.no], "; number of pairwise comparisons = ", no.adjust, "; adjusted error rate = ", sprintf("%.4f", ERate[factor.no]), "\n"))
    } # end (for factor.no)
  } else if (Type1Adj == "NULL") {
    cat("# -- Error Rates for Invariance Tests (No adjustment is made) -- #", rep("\n", 2))
    cat(paste0("Type I Error Rate for Invariance Tests = ", Type1, "\n"))
    for (factor.no in 1: no.factor) {
      ERate[factor.no] <- Type1
    } # end (for factor.no)
  } # end (Type1Adj)
  cat(rep("\n",2))

  if (TYPE == "MonteCarlo") {

    mcmc <- MASS::mvrnorm(n=2000000, mu=par.est, Sigma=simvcov, tol = 1e-6)  # Run 2,000,000 simulations
    bootcoef <- mcmc
    bootno <- nrow(mcmc)  # No. of successful simulated samples
    cat(paste0("Number of Successful Simulated Samples = ", bootno, "\n"))

  } else { # TYPE=Bootstrap

    ## == Simplified bootstrapping model == ##
    Model.Config.Boot <- lavaan::cfa(Model.Long, data = data.source,
                            meanstructure = TRUE,
                            auto.fix.first = FALSE,
                            marker.int.zero = TRUE,
                            ordered = FALSE,
                            missing = "fiml", se = "none", test = "none", check.start = FALSE, check.post = FALSE, check.gradient = FALSE)

    ## == Bootstrapping == ##
    bootcoef <- lavaan::bootstrapLavaan(Model.Config.Boot, R = b.no, ncpus = max(1L, parallel::detectCores() - 1L), parallel="snow")

    # Remove error and nonadmissible bootstrap samples #
    B.na <- attr(bootcoef,"nonadmissible")
    B.er <- attr(bootcoef,"error.idx")
    B.de <- c(B.na,B.er)
    if (length(B.de) != 0) {
      B.re <- bootcoef[-c(B.de),]
      bootcoef <- B.re
    }
    bootno <- nrow(bootcoef)  # number of successful bootstrap samples
    cat(paste0("Number of Successful Bootstrapped Samples = ", bootno, "\n"))
    bootcoef <- bootcoef[,ext]

  } ## end MonteCarlo or Bootstrap


  ## == Start the factor.no loop for CompareLoadings == ##
  for (factor.no in 1: no.factor) {

    alpha <- ERate[factor.no]

    flY <<- matrix(" ",1, (no.group+2))
    flYY <<- matrix(" ",1, (no.group+2))
    EP <<- 1  # estimated parameter number

    no.k <- no.items[factor.no]*no.group  # number of factor loadings in FL.PMI
    FL.PMI <- matrix(1:no.k, nrow = no.group)  # factor loading matrix
    colnames(FL.PMI) <- paste0(' Item ', 1:no.items[factor.no])
    rownames(FL.PMI) <- c(paste0("Group ", 1:no.group))

    for (i in 1:no.group) {
      for (j in 1:no.items[factor.no]) {
        if (j == no.markers[factor.no]) {
          FL.PMI[i, j] <- 1
        } else {
          FL.PMI[i, j] <- format(round(par.est[paste0(names.lv[factor.no],"_G",i,"=~", names.item[[factor.no]][[j]],"_G",i)], digits=4), scientific = FALSE)
        }  ## end if j
      }  ## end loop j
    }  ## end loop i

    class(FL.PMI) <- "numeric"

    for (Referent in 1:no.items[factor.no]) {  ## Loop referent item number for comparison

      ## ==  Print Factor Loadings of all groups  == ##
      no.k <- no.items[factor.no]*no.group  # number of factor loadings in FL
      FL <- matrix(1:no.k, nrow = no.group)  # factor loading matrix
      colnames(FL) <- paste0(' Item ', 1:no.items[factor.no])
      rownames(FL) <- c(paste0("Group ", 1:no.group))

      if (Referent == 1) {  ## if referent is the first item
        for (FL.item.g in 1:no.group) {
          for (FL.item in 1:no.items[factor.no]) {
            if (FL.item == Referent) {
              FL[FL.item.g, FL.item] <- 1
            } else {
              FL[FL.item.g, FL.item] <- FL.PMI[FL.item.g,FL.item]
            }  ## end if FL.item
          }  ## end loop FL.item
        }  ## end loop FL.item.g
      } else {  ## referent is not the first item
        for (FL.item.g in 1:no.group) {
          for (FL.item in 1:no.items[factor.no]) {
            if (FL.item == Referent) {
              FL[FL.item.g, FL.item] <- 1
            } else {
              FL[FL.item.g, FL.item] <- FL.PMI[FL.item.g,FL.item]/FL.PMI[FL.item.g,Referent]
            }  ## end if FL.item
          }  ## end loop FL.item
        }  ## end loop FL.item.g
      }  ## end if referent

#$      cat(rep("\n", 3), "Factor Loadings ", rep("\n", 2))
#$      print(round(FL[], digits=4))

    ## == End (Print Factor Loadings)  == ##

      for (nArg in 1: no.items[factor.no]) {  ## Argument item number for comparison
        no.dif <- factorial(no.group) / factorial(no.group - 2) / 2  # No. of pairwise comparisons for each item
        boot.dif.lx <- matrix(0, nrow(bootcoef), no.dif)  # Create bootstrap difference matrix
        samp.dif.lx <- matrix(0, no.group, no.group) # Create sample difference matrix
        if (nArg == Referent) {
          next
        } else {
          Arg <- nArg
        } ## end if nArg


      ## == Calculate bootstrap difference and sample estimate difference == ##
        comp = 0
        for (r in 1:(no.group-1)) {  ## r is referent group
          for (a in (r+1):no.group) {  ## a is argument group
          comp = comp + 1
          if (Referent == 1) {
            boot.dif.lx[,comp] <- bootcoef[, paste0(names.lv[factor.no],"_G1=~", names.item[[factor.no]][[Arg]],"_G1")] -
                                  bootcoef[, paste0(names.lv[factor.no],"_G2=~", names.item[[factor.no]][[Arg]],"_G2")]
            samp.dif.lx[r,a] <- par.est[paste0(names.lv[factor.no],"_G1=~", names.item[[factor.no]][[Arg]],"_G1")] -
                                par.est[paste0(names.lv[factor.no],"_G2=~", names.item[[factor.no]][[Arg]],"_G2")]
          } else {
            if (Arg == no.markers[factor.no]) {
               boot.dif.lx[,comp] <- 1/bootcoef[, paste0(names.lv[factor.no],"_G1=~", names.item[[factor.no]][[Referent]],"_G1")] -
                                     1/bootcoef[, paste0(names.lv[factor.no],"_G2=~", names.item[[factor.no]][[Referent]],"_G2")]
               samp.dif.lx[r,a] <- 1/par.est[paste0(names.lv[factor.no],"_G1=~", names.item[[factor.no]][[Referent]],"_G1")] -
                                   1/par.est[paste0(names.lv[factor.no],"_G2=~", names.item[[factor.no]][[Referent]],"_G2")]
            } else { # Arg != no.markers[factor.no]
              boot.dif.lx[,comp] <- bootcoef[,paste0(names.lv[factor.no],"_G1=~", names.item[[factor.no]][[Arg]],"_G1")]/
                                    bootcoef[,paste0(names.lv[factor.no],"_G1=~", names.item[[factor.no]][[Referent]],"_G1")] -
                                    bootcoef[,paste0(names.lv[factor.no],"_G2=~", names.item[[factor.no]][[Arg]],"_G2")]/
                                    bootcoef[,paste0(names.lv[factor.no],"_G2=~", names.item[[factor.no]][[Referent]],"_G2")]
              samp.dif.lx[r,a] <- par.est[paste0(names.lv[factor.no],"_G1=~", names.item[[factor.no]][[Arg]],"_G1")]/
                                  par.est[paste0(names.lv[factor.no],"_G1=~", names.item[[factor.no]][[Referent]],"_G1")] -
                                  par.est[paste0(names.lv[factor.no],"_G2=~", names.item[[factor.no]][[Arg]],"_G2")]/
                                  par.est[paste0(names.lv[factor.no],"_G2=~", names.item[[factor.no]][[Referent]],"_G2")]
            }  ## end if Arg
          }  ## end if Referent
        }  ## end loop a
      }  ## end loop r

#$ print(samp.dif.lx)


      ## == Calculate Percentile Probability == ##

      F1.comp.pp <- matrix(0, no.group, no.group)  ## Percentile Probability
      pno_nipair <- 0  ## Set number of non-invariant pairs to zero

      if (quantile(boot.dif.lx[, 1], probs = 0.5, na.rm = TRUE) > 0) {
        F1.comp.pp[1,2] = 2*(sum(boot.dif.lx[, 1] < 0, na.rm=TRUE)/bootno)
      } else {
        F1.comp.pp[1,2] = 2*(sum(boot.dif.lx[, 1] > 0, na.rm=TRUE)/bootno)
      }  ## end if


      ## == Percentile Confidence Intervals of Pairwise Comparisons == ##

      PCI <- matrix(1:(no.dif*10), nrow = no.dif)
      colnames(PCI) <- c("Group 1","Group 2","0.5%","2.5%","5%","Estimate","95%","97.5%","99.5%","p-value")

      comp = 0
      for (r in 1:(no.group-1)) {
        for (a in (r+1):no.group) {
          comp = comp + 1
          PCI[comp, 1] <- r
          PCI[comp, 2] <- a
          PCI[comp, 3] <- quantile(boot.dif.lx[, comp],c(0.005), na.rm = TRUE)
          PCI[comp, 4] <- quantile(boot.dif.lx[, comp],c(0.025), na.rm = TRUE)
          PCI[comp, 5] <- quantile(boot.dif.lx[, comp],c(0.05), na.rm = TRUE)
          PCI[comp, 6] <- samp.dif.lx[r,a]
          PCI[comp, 7] <- quantile(boot.dif.lx[, comp],c(0.95), na.rm = TRUE)
          PCI[comp, 8] <- quantile(boot.dif.lx[, comp],c(0.975), na.rm = TRUE)
          PCI[comp, 9] <- quantile(boot.dif.lx[, comp],c(0.995), na.rm = TRUE)
          PCI[comp,10] <- F1.comp.pp[r,a]
          if (F1.comp.pp[r,a] < alpha) {
            if (pno_nipair == 0)	{
              pnipair <- c(r,a)
              pno_nipair <- pno_nipair + 1
            } else {
              pnipair <- c(pnipair,r,a)
              pno_nipair <- pno_nipair + 1
            }  ## end if pno_nipair
          }  ## end if F1.comp.pp
        }  ## end loop a
      }  ## end loop r

#$      cat("\n")
#$      if (Referent > Arg) {
#$        cat(paste0("Referent Item: ", Referent, "  Argument Item: ", Arg, "\n"))
#$      } else {
#$        cat(paste0("Referent Item: ", Referent, "  Argument Item: ", (Arg+1), "\n"))
#$      }  ## end if Referent

#$      cat("\n")

#$    print(paste0("Factor = ", factor.no, "  Referent = ", Referent, " Argument = ", Arg))
#$    print(round(PCI[], digits=4))

      ## == Run List and Delete == ##
#$      cat(rep("\n", 2), "Sets of Groups with Invariant Loadings", "\n")
        if (pno_nipair == no.dif) {

          flZ <<- matrix(" ", 1, (no.group+2))
          flX <<- matrix(1:(no.group+2), 1, (no.group+2))
          flX[1, 1:no.group] <- 0
          for (temp.flZ in 1: no.group) {
            flZ[1,temp.flZ] <- paste0("F", factor.no, "L", EP)
            EP <<- EP + 1
          }
          flX[1, (no.group+1)] <- Referent
          flZ[1, (no.group+1)] <- Referent
#?          if (Referent > Arg) {
            flX[1, (no.group+2)] <- Arg
            flZ[1, (no.group+2)] <- Arg
#?          } else {
#?            flX[1, (no.group+2)] <- Arg + 1
#?            flZ[1, (no.group+2)] <- Arg + 1
#?          }  ## end if Referent
#$        print(flX, quote=FALSE)
          flY <<- rbind(flY,flX)
#$        print(flZ, quote=FALSE)
          flYY <<- rbind(flYY,flZ)

        } else if (pno_nipair > 0) {

          flX <<- matrix(" ", 1, (no.group+2))
          flZ <<- matrix(" ", 1, (no.group+2))

          listanddelete.lx(factor.no, no.group, no_nipair=pno_nipair, nipair=pnipair, Referent, Arg) ## Run list and delete lx

        } else {  # (pno_nipair = 0)

          flX <- matrix(1:(no.group+2), 1, (no.group+2))
          flZ <<- matrix(" ", 1, (no.group+2))
          flZ <- matrix(paste0("F", factor.no, "L", EP), 1, (no.group+2))
          EP <<- EP + 1
          flX[1, (no.group+1)] <- Referent
          flZ[1, (no.group+1)] <- Referent
          flX[1, (no.group+2)] <- Arg
          flZ[1, (no.group+2)] <- Arg
          flY <<- rbind(flY,flX)
#$        print(flZ, quote=FALSE)
          flYY <<- rbind(flYY,flZ)

        }  ## end if pno_nipair
      } ## end loop Arg
    } ## end loop Referent



  if (length(unique(flY[, no.group+1])) == 2 & no.items[factor.no] == 2) { # Only 2 groups with more than 2 items
    flX <<- matrix(1:(no.group+2), 1, (no.group+2))
    flZ <<- matrix(" ", 1, (no.group+2))
    flYa <<- flY[nrow(flY), no.group+1]
    flYb <<- flY[nrow(flY), no.group+2]
    flX[1, no.group + 1] <- flYb
    flX[1, no.group + 2] <- flYa
    flY <<- rbind(flY,flX)

    flZ[1, no.group + 1] <- flYb
    flZ[1, no.group + 2] <- flYa
    for (temp.flYY in 1: no.group) {
      flZ[1, temp.flYY] <- paste0("F", factor.no, "I", EP)
      EP <<- EP + 1
    }
    flYY <<- rbind(flYY,flZ)
  }

  flY <<- flY[-c(1),]
#$  cat("\n")
#$  print(flY, quote = FALSE)

  flYY <<- flYY[-c(1),]
#$  cat("\n")
#$  print(flYY, quote = FALSE)

  if (is.matrix(flYY) == FALSE) { flYY <<- matrix(flYY, nrow=1) }

  NIcombine <- table(flY[, no.group+1], flY[, no.group+2])
#  NI.level <- unique(flY[, no.group+1])
#  NIcombine <- NIcombine[NI.level, NI.level]


  tempR <- unique(flYY[, no.group+1])

    for (aR in 1:length(tempR)) {
      R <- tempR[aR]
      for (r in 1:no.items[factor.no]) {
        if (r == 1) {
          assign(paste0("Rset", r), c(length(which(flYY[,no.group+1] < R))+1:(NIcombine[R, r])))
        } else if (r == R) {
          assign(paste0("Rset", r), 0)
        } else {
          assign(paste0("Rset", r), c(length(which(flYY[,no.group+1] < R))+sum(NIcombine[R, 1:(r-1)])+1:(NIcombine[R, r])))
        }
        if (R == 1) { assign("Rset1", 0) }
      }  ## end loop r

      ## expand.grid(Rset1, Rset2, Rset3, Rset4, Rset5 ...)

      temp.row <- eval(parse(text = paste0("expand.grid(", paste0(c("Rset"), 1: no.items[factor.no], collapse=","),")")))
      if (aR == 1) {
        Set.row <- temp.row
      } else {
        Set.row <- rbind(Set.row, temp.row)
      } ## end if R
    }  ## end loop aR

    Model.load <- array(" ", dim = c(no.items[factor.no], no.group, nrow(Set.row)))  # all possible invariance models

    for (Xset in 1: nrow(Set.row)) {
      a <- matrix(" ", no.items[factor.no], no.group)
      for (r in 1:no.items[factor.no]) {
        if (Set.row[Xset, r] == 0) {
          a[r, ] <- 1
        } else {
          a[r, ] <- flYY[Set.row[Xset, r], 1:no.group]
        } ## end if Set.row
#?        if (flYY[r, no.group+1] > 1) { break }
      }  ## end loop r
      Model.load[, , Xset] <- a
    }  ## end loop Xset

#$ print(Model.load)

  ## == Model.summary == ##
  # Column 1 - number of estimated loadings
  # Column 2 - number of measurement invariance items
  # Column 3 - p-value of Model chi-square
  # Column 4 - CFI
  # Column 5 - RMSEA
  # Column 6 - SRMR

    Model.summary <- matrix(0, nrow(Set.row), 7)
    for (R in 1:nrow(Set.row)) {
      Model.summary[R, 1] <- length(unique(c(Model.load[,,R]))) - 1
      no.MI <- no.items[factor.no] - 1
      for (Ra in 1: no.items[factor.no]) {
        if (length(unique(c(Model.load[Ra,,R]))) > 1) {no.MI <- no.MI - 1}
      }  ## end loop Ra
      Model.summary[R, 2] <- no.MI
    }  ## end loop R



    ## == Recommended Models (Rec.Model) - Largest no. of MI items, then smallest number of estimated loadings == ##
    Rec.Model <- which(Model.summary[,2] == max(Model.summary[,2]))
    Rec.Model <- Rec.Model[which(Model.summary[Rec.Model,1] == min(Model.summary[Rec.Model,1]))]


    ## == Find recommended model if no.items > 2 == ##
    if (no.items[factor.no] > 2) {

      Rec.Model.load <- matrix(nrow = length(Rec.Model), ncol = length(lavaan::coef(Model.Long.config)))
      MODFIT <- matrix(0, nrow = length(Rec.Model), ncol = 9)
      colnames(MODFIT) <-
        c("  chisq.scaled", "  df.scaled", "  pvalue.scaled", "  rmsea.scaled", "  cfi.scaled", "  tli.scaled", "  srmr", "  aic", "  bic")
      Rec.Model.FL <- array(0, dim = c(no.group, no.items[factor.no], length(Rec.Model)))
      colnames(Rec.Model.FL) <- paste0(" Item ", 1:no.items[factor.no])
      rownames(Rec.Model.FL) <- c(paste0("Group ", 1:no.group))

      for (Model.R in 1: length(Rec.Model)) {

        ## == Specify the partial metric invariance models (PMI.X) == ##
        PMI.X.1 <- paste0("PMI.X <- '", "\n")
        for (j in 1:2) {
          PMI.X.1 <- paste0(PMI.X.1, "  ",
                          names.lv[factor.no], "_G",j," =~ ", paste0(Model.load[1, j, Rec.Model[Model.R]]), "*", names.item[[factor.no]][[1]], "_G",j," + ")
          for (i in 2: (no.items[factor.no]-1)) {
            PMI.X.1 <- paste0(PMI.X.1, Model.load[i, j, Rec.Model[Model.R]],"*", names.item[[factor.no]][[i]], "_G",j," + ")
          } ## end loop i
          PMI.X.1 <- paste0(PMI.X.1, Model.load[i+1, j, Rec.Model[Model.R]],"*", names.item[[factor.no]][[i+1]], "_G",j, "\n")
        }
        ## -- Create Residual Covariance -- ##

  ## == Create Residual Covariance == ##
        PMI.X.1 <- paste0(PMI.X.1, "\n")
        for (i in 1:no.items[factor.no]) {
          PMI.X.1 <- paste0(PMI.X.1, "  ", names.item[[factor.no]][[i]], "_G1 ~~ ", names.item[[factor.no]][[i]], "_G2", "\n")
        }

        PMI.X.1 <- paste0(PMI.X.1, "'", "\n")
        eval(parse(text = PMI.X.1))


        ## == Run PMI.X == ##
        if (Cluster == "NULL") {
          eval(parse(text=   "PMI.X.fit <- lavaan::sem(PMI.X,
                             data.source,
                             missing = 'fiml',
                             auto.fix.first = FALSE,
                             marker.int.zero = TRUE,
                             meanstructure = T,
                             information = 'observed',
                             estimator = 'MLR')"))
        } else {
          eval(parse(text=   "PMI.X.fit <- lavaan::sem(PMI.X,
                             data.source,
                             missing = 'fiml',
                             auto.fix.first = FALSE,
                             marker.int.zero = TRUE,
                             meanstructure = T,
                             information = 'observed',
                             cluster = Cluster,
                             estimator = 'MLR')"))
        }


        ## == Request summary outputs == ##
 #$     eval(parse(text = "print(lavaan::summary(PMI.X.fit, fit.measure = T, standardized = T, rsq = T))"))

        ## == Save fit indices == ##
        XX <- lavaan::fitMeasures(PMI.X.fit, fit.measures=
              c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr", "aic", "bic"))
        MODFIT[Model.R,] <- format(round((XX), digits = 4), nsmall = 4, scientific = FALSE)


        aa <- lavaan::lavInspect(PMI.X.fit, what='est')$lambda
        ## == Save Factor Loadings == ##
        for (G in 1: no.group) {
          Rec.Model.FL[G, , Model.R] <- aa[(no.items[factor.no]*(G-1)+1):(no.items[factor.no]*G),G]
        } ## end loop G

      } ## end loop Model.R


      ## == Print model fit indices == ##
      rownames(MODFIT) <- c(paste0("Model ", 1: length(Rec.Model)))
#$    cat(rep("\n", 3))
#$    cat("####    MODEL FIT INDICES   ####", rep("\n", 2))
#$    print(MODFIT, quote=FALSE, right=TRUE)

    ## == Print factor loadings == ##
#$    for (Model.R in 1: length(Rec.Model)) {
#$      cat(rep("\n", 2), "Factor Loadings of Model", Model.R, rep("\n", 2))
#$      print(round(Rec.Model.FL[,,Model.R], digits=4))
#$    }  ## end loop Model.R

      ## == First model with best model fit (R.Model) == ##
      class(MODFIT) <- "numeric"
      if (BMSC == "Chisquare") {
        Model.R <- which.min(MODFIT[,1])
      } else if (BMSC == "CFI") {
        Model.R <- which.max(MODFIT[,5])
      } else if (BMSC == "RMSEA") {
        Model.R <- which.min(MODFIT[,4])
      } else if (BMSC == "SRMR") {
        Model.R <- which.min(MODFIT[,7])
      }
      Model.R <- Model.R[1]
      S.Model <- Model.load[,,Rec.Model[Model.R]]  # Recommended Model
      un.load <- sort(subset(unique(c(S.Model)), unique(c(S.Model))[] != 1))
      no.un.load <- length(un.load)
      R.Model <- S.Model
      for (j in 1:no.un.load) { R.Model[S.Model == un.load[j]] <- paste0("F", factor.no, "L", j) }  ## end loop j

    ## == Save Recommended Model (Recommend.Model) == ##
      if (factor.no == 1) {
        Recommend.Model <- matrix("PMI.Model.R <- '", 1)
      }
      Recommend.Model <- rbind(Recommend.Model, "\n")
      for (j in 1:2) {
        PMI <- paste0("   ", names.lv[factor.no], "_G",j," =~ ", paste0(Model.load[1, j, Rec.Model[Model.R]]), "*", names.item[[factor.no]][[1]], "_G",j," + ")
        for (i in 2: (no.items[factor.no]-1)) {
          PMI <- paste0(PMI, R.Model[i,j],"*", names.item[[factor.no]][[i]], "_G",j," + ")
        } ## end loop i
        PMI <- paste0(PMI, R.Model[i+1,j],"*", names.item[[factor.no]][[i+1]], "_G",j, "\n")
      Recommend.Model <- rbind(Recommend.Model, PMI)
      }

      ## -- Create Residual Covariance -- ##
      Recommend.Model <- rbind(Recommend.Model, "\n")
      for (i in 1:no.items[factor.no]) {
        Recommend.Model <- rbind(Recommend.Model, paste0("   ", names.item[[factor.no]][[i]], "_G1 ~~ ", names.item[[factor.no]][[i]], "_G2", "\n"))
      }
    }  ## end if no.items[factor.no] > 2


    if (no.items[factor.no] < 3) {   ## Begin Models with 2 items ##
      R.Model <- Model.load[,,1]  # Recommended Model
      un.load <- subset(unique(c(R.Model)), unique(c(R.Model))[] != 1)
      no.un.load <- length(un.load)
      for (j in 1:no.un.load) { R.Model[R.Model == un.load[j]] <- paste0("F", factor.no, "L", j) }  ## end loop j

      ## == Save Recommended Model (Recommend.Model) == ##

      if (factor.no == 1) {
        Recommend.Model <- matrix("PMI.Model.R <- '", 1)
      }
      Recommend.Model <- rbind(Recommend.Model, "\n")
      for (j in 1:2) {
        PMI <- paste0("   ", names.lv[factor.no], "_G",j," =~ ", paste0(Model.load[1, j, Rec.Model[Model.R]]), "*", names.item[[factor.no]][[1]], "_G",j," + ")
        PMI <- paste0(PMI, R.Model[2,j],"*", names.item[[factor.no]][[2]], "_G",j, "\n")
      Recommend.Model <- rbind(Recommend.Model, PMI)
      }

      ## -- Create Residual Covariance -- ##
      Recommend.Model <- rbind(Recommend.Model, "\n")
      for (i in 1:2) {
        Recommend.Model <- rbind(Recommend.Model, paste0("   ", names.item[[factor.no]][[i]], "_G1 ~~ ", names.item[[factor.no]][[i]], "_G2", "\n"))
      }
      Recommend.Model <- rbind(Recommend.Model, "\n")

    }  ## End (Models with 2 items)

  } ## End (loop factor.no)

  Recommend.Model <- rbind(Recommend.Model, "  '", "\n")  ## last line of Recommended Model


  cat("\n", "## =====  Recommended Model  ===== ##", rep("\n", 2))
  cat(Recommend.Model)
#  for (i in 1: nrow(Recommend.Model)) { cat(Recommend.Model[i], "\n") } ## end loop i
  cat(rep("\n", 2))
  cat("## ===== lavaan Outputs of PMI.Model.R ===== ##")
  cat(rep("\n", 2))


  ## == Start print recommended model to PMI.txt (Partial Metric Invariance Model) == ##
  sink('PMI.txt')

  cat("## =====  Recommended Model  ===== ##")
  cat("\n")
  cat(Recommend.Model)
#  for (i in 1: nrow(Recommend.Model)) { cat(Recommend.Model[i], "\n") } ## end loop i

  cat(rep("\n",2), "## Run PMI.Model.R")
  cat(rep("\n",2), paste0("  PMI.Model.fit <- lavaan::sem(PMI.Model.R,"))
  cat("\n", paste0("    ", arg2_char, ","))
  cat("\n", "    missing = 'fiml',")
  cat("\n", "    auto.fix.first = FALSE,")
  cat("\n", "    information = 'observed',")
  if (Cluster != "NULL") { cat("\n", paste0("    cluster = ", arg4_char,",")) }
  cat("\n", "    estimator = 'MLR')")

  cat(rep("\n", 3), "## Request summary outputs")
  cat(rep("\n", 3), paste0("  print('## ===== Final Model -- PMI.Model.R ===== ##')", "\n"))
  cat(rep("\n", 2), paste0("  print(lavaan::summary(PMI.Model.fit, fit.measure = T, standardized = T, rsq = T))"), "\n")

  sink()  ## Stop writing to file

#  source("PMI.txt")
  PMI.Model.R <<- eval(parse(text=Recommend.Model))

  ## == Run PMI.Model.R == ##
  if (Cluster == "NULL") {
    eval(parse(text=   "PMI.Model.fit <- lavaan::sem(PMI.Model.R,
                       data.source,
                       missing = 'fiml',
                       auto.fix.first = FALSE,
                       marker.int.zero = TRUE,
                       meanstructure = T,
                       information = 'observed',
                       estimator = 'MLR')"))
  } else {
    eval(parse(text=   "PMI.Model.fit <- lavaan::sem(PMI.Model.R,
                       data.source,
                       missing = 'fiml',
                       auto.fix.first = FALSE,
                       marker.int.zero = TRUE,
                       meanstructure = T,
                       information = 'observed',
                       cluster = Cluster,
                       estimator = 'MLR')"))
  }

  print(lavaan::summary(PMI.Model.fit, fit.measure = T, standardized = T, rsq = T))


  cat(rep("\n", 2), "## Compare fit indices across models ##", "\n")
  DMODFIT <- matrix(0, nrow = 3, ncol = 9)
  colnames(DMODFIT) <-
        c("  chisq.scaled", "  df.scaled", "  pvalue.scaled", "  rmsea.scaled", "  cfi.scaled", "  tli.scaled", "  srmr", "  aic", "  bic")
  rownames(DMODFIT) <- c("Model.config", "PMI.Model.fit", "PMI.Model.fit - Model.config")
  XX1 <- lavaan::fitMeasures(Model.Long.config, fit.measures=
      c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr", "aic", "bic"))
  DMODFIT[1,] <- XX1
  XX2 <- lavaan::fitMeasures(PMI.Model.fit, fit.measures=
      c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr", "aic", "bic"))
  DMODFIT[2,] <- XX2

  DMODFIT[3,1] <- DMODFIT[2,1] -DMODFIT[1,1]
  DMODFIT[3,2] <- DMODFIT[2,2] -DMODFIT[1,2]
  DMODFIT[3,3] <- pchisq(q=DMODFIT[3,1], df=DMODFIT[3,2], lower.tail=FALSE)
  DMODFIT[3,4] <- DMODFIT[2,4] -DMODFIT[1,4]
  DMODFIT[3,5] <- DMODFIT[2,5] -DMODFIT[1,5]
  DMODFIT[3,6] <- DMODFIT[2,6] -DMODFIT[1,6]
  DMODFIT[3,7] <- DMODFIT[2,7] -DMODFIT[1,7]
  DMODFIT[3,8] <- DMODFIT[2,8] -DMODFIT[1,8]
  DMODFIT[3,9] <- DMODFIT[2,9] -DMODFIT[1,9]
  print(round(DMODFIT, digits=4))


  ## == Print Factor Loadings == ##
  Final.Model.FL <- matrix(0, sum(no.items), 2)
  rownames(Final.Model.FL) <- paste0(" Item ", names.ind, "    ")
  colnames(Final.Model.FL) <- c(paste0("Group ", 1:no.group))

  Final.par.est <- lavaan::parameterEstimates(PMI.Model.fit)
  ext <- c(which(Final.par.est[,"op"] == "=~"))
  aa <- Final.par.est[ext, "est"]

  ## Extract factor loadings ##
  n <- 1
  for (j in 1: 2) {
    Final.Model.FL[1:no.items[1], j] <- aa[(no.items[1]*(j-1)+1):(no.items[1]*j)]
  }
  n <- n + no.items[1]*2

  if(no.factor > 1){
    for (i in 2:no.factor) {
      for (j in 1: 2) {
        Final.Model.FL[(sum(no.items[1:(i-1)])+1):sum(no.items[1:i]), j] <- aa[(n+no.items[i]*(j-1)):(n+(no.items[i]*j-1))]
      }
      n <- n + no.items[i]*2
    }
  }

  cat("\n")
  cat("## === Factor Loadings in Final Model === ##")
  cat("\n")
  print(round(Final.Model.FL, digits=4))
  ## == End (Print Factor Loadings) == ##

  cat(rep("\n",2),"The recommended model PMI.Model.R is saved in the file 'PMI.txt'", "\n")

} ## End (Function PAIRCompareLoadings)

# ==================== Finish Function "PAIRCompareLoadings" ==================== #





# ==================== Create Function "PAIRCompareMeans" ==================== #
#' Scalar Invariance Test and Compare Latent Means in Congruence/Fit/Agreement Studies
#'
#' Conduct scalar invariance test, identify a partial scalar invariance model, and compare latent means for data from two non-independent sources.
#'
#' Requires defining the measurement model only once. All indicators should have the suffix _G1 and _G2 (e.g., x1_G1, x1_G2, x2_G1, and x2_G2) in the data file to indicate from which source the indicators were measured.
#'
#' Residuals of indicators are covaried across sources automatically.
#'
#' Reference for correct.cfi: Widaman, K. F., & Thompson, J. S. (2003). On specifying the null model for incremental fit indices in structural equation modeling. Psychological Methods, 8(1), 16-37.
#'
#' @param model.PMI Partial metric invariance model (PMI.Model.R from PAIRCompareLoadings() or user-specified).
#' @param data.source A data frame containing the observed variables used in the model.
#' @param Cluster Cluster variable for nested data. The Monte Carlo simulation method should be used for nested data.
#' @param Bootstrap Number of bootstrap samples, must be between 500 and 10,000. If not specified, 1 million Monte Carlo simulated samples (default) will be used.
#' @param Type1 Overall Type I error rate (default is 0.05) for identifying non-invariant items in the List and Delete method.
#' @param Type1Adj Adjustment of Type I error rate for multiple tests for each construct. Default is Possible False Discovery Rate (PFDR = Type I error rate/No. of freely estimated factor loadings/intercepts across all waves), can also use (BON)ferroni adjustment (Type I error rate/No. of pairwise comparisons), or "NULL".
#' @param BMSC Best Model Selection Criterion used to select among models with the same number of estimated parameters. "ChiSquare" = minimum Chi-Square value, "CFI" = maximum CFI value, "RMSEA" = minimum RMSEA value, and "SRMR" = minimum SRMR value (default).
#' @param correct.cfi If TRUE (default), baseline model will be corrected for calculation of incremental fit indices based on Widaman & Thompson (2003).
#' @return Partial scalar invariance model in PSI.txt file and results of latent means comparisons.
#' @export
#' @examples
#'
#'
#' ## == Example C - Non-independent Data from Two Sources == ##
#'
#' # Data file is "Example.C"
#'
#' ## Not run:
#' ## Specify the measurement model - Model.C ##
#' Model.C <- '
#'      External =~ x1 + x2 + x3 + x4 + x5 + x6
#'      Internal =~ x7 + x8 + x9 + x10
#' '
#'
#' ## ===== Compare Factor Loadings ===== ##
#' PAIRCompareLoadings(Model.C, Example.C, Type1=0.01, Type1Adj="PFDR", BMSC="SRMR")
#' ## End (Not run)
#'
#' ## ===== Compare Means ===== ##
#' PAIRCompareMeans(PMI.Model.R, Example.C, Type1=0.01, Type1Adj="PFDR", BMSC="SRMR")
#'
PAIRCompareMeans <- function(model.PMI, data.source, Cluster="NULL", Bootstrap=0, Type1 = 0.01, Type1Adj = "PFDR", BMSC="SRMR", correct.cfi=TRUE) {

# model.PMI <- PMI.Model.R
# data.source <- Example.C
# Cluster="NULL"
# Bootstrap=0
# Type1 = 0.01
# Type1Adj = "PFDR"
# BMSC="SRMR"
# correct.cfi=TRUE

  Type1Adj <- toupper(Type1Adj)
  match.arg(Type1Adj, c("PFDR","BON", "NULL"))
  BMSC <- toupper(BMSC)
  match.arg(BMSC, c("CHISQUARE","CFI", "RMSEA", "SRMR"))
  if (is.logical(correct.cfi) == FALSE) stop("correct.cfi (corrected baseline model for calculating incremental fit indices) can only be TRUE or FALSE")

  arg1_char <- deparse(substitute(model.PMI))
  arg2_char <- deparse(substitute(data.source))
  arg4_char <- deparse(substitute(Cluster))

  no.group <- 2

  parsed <- lavParseModelString(model.PMI, as.data.frame = TRUE)
  parsed <- parsed[, c("lhs", "op", "rhs")]
  parsed <- parsed[parsed$op == "=~",]
  parsed[,"lhs"] <- gsub("_G.*", "", parsed[,"lhs"])
  parsed[,"rhs"] <- gsub("_G.*", "", parsed[,"rhs"])
  names.lv <- unique(parsed[,"lhs"])
  names.ind <- unique(parsed[, "rhs"])
  no.factor <- length(names.lv)
  no.items.g <- length(names.ind) # number of items per group
  no.items <- matrix(1:no.factor, nrow = 1)  # number of items per factor
  for (factor.no in 1:no.factor) { no.items[factor.no] <- sum(parsed[,"lhs"] == names.lv[factor.no] & parsed[,"op"] == "=~")/2 } ## end loop factor.no

  names.item <- list()
  for (i in 1:no.factor) {
    names.item[[i]] <- unique(parsed[which(parsed[,"lhs"] == names.lv[i] & parsed[,"op"] == "=~"), "rhs"])
  }

  ## Check for bootstrap sample number (Bootstrap) ##
  if (Bootstrap !=0) {
    b.no.integer <- Bootstrap == round(Bootstrap)
    if (b.no.integer == "FALSE") stop("Bootstrap sample number must be an integer")
    if (Bootstrap > 10000) stop("Bootstrap sample number greater than 10,000 is not recommended")
    if (Bootstrap < 500) stop("Bootstrap sample number smaller than 500 is not recommended")
    TYPE = "Bootstrap"
    b.no <- Bootstrap
  } else {
    TYPE = "MonteCarlo"
  } # end (Bootstrap != 0)

  count.tx <- 0


  ## ========== Run model.PMI ========== ##

  if (Cluster == "NULL") {
    PMI.Model.fit <- lavaan::sem(model.PMI,
      data.source,
      missing = 'fiml',
      auto.fix.first = FALSE,
      marker.int.zero = TRUE,
      meanstructure = T,
      information = 'observed',
      estimator = 'MLR')
  } else {
      PMI.Model.fit <- lavaan::sem(model.PMI,
      data.source,
      missing = 'fiml',
      auto.fix.first = FALSE,
      marker.int.zero = TRUE,
      meanstructure = T,
      information = 'observed',
      cluster = Cluster,
      estimator = 'MLR')
   }  # end Cluster

  ## Request summary outputs
  cat(rep("\n", 2), "## ===== Partial Metric Invariance Model ===== ##", "\n")
  print(lavaan::summary(PMI.Model.fit, fit.measure = T, standardized = T, rsq = T))
  cat("\n")

  ## ===== End (Run model.PMI) ===== ##

  par.est <- lavaan::coef(PMI.Model.fit)  # sample parameters
  group.names <<- c(paste0("Group ", 1:no.group))  # group names

  # Find out the marker item per factor #
  temp <- lavaan::parameterEstimates(PMI.Model.fit)
  no.markers <- matrix(1:no.factor, nrow = 1)  # location of marker items
  for (factor.no in 1:no.factor) {
    no.markers[factor.no] <- which(names.item[[factor.no]] == sub("_G1$", "", temp[which(temp[,"lhs"] == paste0(names.lv[factor.no],"_G1") & temp[,"op"] == "=~" & temp[,"est"] == 1), "rhs"]))
  } ## end for factor.no

  temp <- lavaan::parameterEstimates(PMI.Model.fit, remove.nonfree=TRUE)
  simvcov <- lavInspect(PMI.Model.fit, what="vcov")
  par.est <- lavaan::coef(PMI.Model.fit)  # sample parameters

  ## -- Is change names required ? -- ##
  names(par.est) <- paste0(temp[,"lhs"], temp[,"op"], temp[,"rhs"])
  colnames(simvcov) <- paste0(temp[,"lhs"], temp[,"op"], temp[,"rhs"])
  rownames(simvcov) <- paste0(temp[,"lhs"], temp[,"op"], temp[,"rhs"])

  mext <- which(temp[,"op"] == "=~" | temp[, "op"] == "~1")
  par.est <- par.est[mext]
  simvcov <- simvcov[mext,mext]
  no.par.g <- length(par.est)/no.group  # number of estimated parameters per group #
  no.lx.g <- sum(temp[, "op"] == "=~")/2 # number of estimated LX per group #

  cat(rep("\n", 3), "## ======= SCALAR INVARIANCE ANALYSIS ======= ##", rep("\n", 2))  ## print heading

  # Calculate alpha for factor[factor.no] #
  ERate <- matrix(1:no.factor, nrow = 1)  # error rate for each factor

  cat("# -- Overall Type I Error Rates for Invariance Tests for Each Factor = ", Type1, " -- #", rep("\n", 2))

  if (Type1Adj == "PFDR") {  ## print adjusted error rate
    cat("# -- Adjusted Type I Error Rates for Invariance Tests (Adjusted by Possible False Discovery Rate) -- #", rep("\n", 2))
    for (factor.no in 1: no.factor) {
      no.adjust <- (no.items[factor.no]-1)*(no.group-1)
      ERate[factor.no] <- Type1 / no.adjust
      cat(paste0("Factor ", factor.no, ": Number of items = ", no.items[factor.no], "; number of free intercepts = ", no.adjust, "; adjusted error rate = ", sprintf("%.4f", ERate[factor.no]), "\n"))
    } # end (for factor.no)
  } else if (Type1Adj == "BON") {
    cat("# -- Adjusted Type I Error Rates for Invariance Tests (Bonforroni adjustment for number of pairwise comparisons) -- #", rep("\n", 2))
    for (factor.no in 1: no.factor) {
      no.adjust <- (no.items[factor.no]-1)*(no.group*(no.group-1)/2)
      ERate[factor.no] <- Type1 / no.adjust
      cat(paste0("Factor ", factor.no, ": Number of items = ", no.items[factor.no], "; number of pairwise comparisons = ", no.adjust, "; adjusted error rate = ", sprintf("%.4f", ERate[factor.no]), "\n"))
    } # end (for factor.no)
  } else if (Type1Adj == "NULL") {
    cat("# -- Error Rates for Invariance Tests (No adjustment is made) -- #", rep("\n", 2))
    cat(paste0("Type I Error Rate for Invariance Tests = ", Type1, "\n"))
    for (factor.no in 1: no.factor) {
      ERate[factor.no] <- Type1
    } # end (for factor.no)
  } # end (Type1Adj)
  cat(rep("\n",2))

  if (TYPE == "MonteCarlo") {

    mcmc <- MASS::mvrnorm(n=2000000, mu=par.est, Sigma=simvcov, tol = 1e-6)  # Run 2,000,000 simulations
    bootcoef <- mcmc
    bootno <- nrow(mcmc)  # No. of successful simulated samples
    cat(paste0("Number of Successful Simulated Samples = ", bootno, "\n"))

  } else { # TYPE=Bootstrap

    ## == Simplified bootstrapping model == ##
    PMI.boot <- lavaan::cfa(model.PMI, data = data.source,
                meanstructure = TRUE,
                auto.fix.first = FALSE,
                marker.int.zero = TRUE,
                ordered = FALSE,
                missing = "fiml", se = "none", test = "none", check.start = FALSE, check.post = FALSE, check.gradient = FALSE)

    ## == Bootstrapping == ##
    bootcoef <- lavaan::bootstrapLavaan(PMI.boot, R = b.no, ncpus = max(1L, parallel::detectCores() - 1L), parallel="snow")

    # Remove error and nonadmissible bootstrap samples #
    B.na <- attr(bootcoef,"nonadmissible")
    B.er <- attr(bootcoef,"error.idx")
    B.de <- c(B.na,B.er)
    if (length(B.de) != 0) {
      B.re <- bootcoef[-c(B.de),]
      bootcoef <- B.re
    }
    bootno <- nrow(bootcoef)  # number of successful bootstrap samples
    cat(paste0("Number of Successful Bootstrapped Samples = ", bootno, "\n"))
    bootcoef <- bootcoef[,ext]

  } ## end MonteCarlo or Bootstrap


  ## == Start the factor.no loop for CompareIntercepts == ##
  for (factor.no in 1: no.factor) {

    alpha <- ERate[factor.no]

    txY <<- matrix(" ",1, (no.group+2))
    txYY <<- matrix(" ",1, (no.group+2))
    EP <<- 1  # estimated parameter number

    no.k <- no.items[factor.no]*no.group  # number of factor loadings in FL.PMI
    FL.PMI <- matrix(1:no.k, nrow = no.group)  # factor loading matrix
    colnames(FL.PMI) <- paste0(' Item ', 1:no.items[factor.no])
    rownames(FL.PMI) <- c(paste0("Group ", 1:no.group))

    TX.PMI <- matrix(1:no.k, nrow = no.group)  # Intercept matrix
    colnames(TX.PMI) <- paste0(' Item ', 1:no.items[factor.no])
    rownames(TX.PMI) <- c(paste0("Group ", 1:no.group))

    for (i in 1:no.group) {
      for (j in 1:no.items[factor.no]) {
        if (j == no.markers[factor.no]) {
          FL.PMI[i, j] <- 1
          TX.PMI[i, j] <- 0
        } else {
          FL.PMI[i, j] <- format(round(par.est[paste0(names.lv[factor.no],"_G",i,"=~", names.item[[factor.no]][[j]],"_G",i)], digits=4), scientific = FALSE)
          TX.PMI[i, j] <- format(round(par.est[paste0(names.item[[factor.no]][[j]],"_G", i, "~1")], digits=4), scientific = FALSE)
        }  ## end if j
      }  ## end for j
    }  ## end for i

    class(FL.PMI) <- "numeric"
    class(TX.PMI) <- "numeric"

    for (Referent in 1:no.items[factor.no]) {  ## Loop referent item number for comparison

      if(length(unique(FL.PMI[,Referent])) != 1) {
        txZ <<- matrix(" ", 1, (no.group+2))
        txX <<- matrix(1:(no.group+2), 1, (no.group+2))
        txX[1, 1:no.group] <- 0
        uni.fl <- unique(FL.PMI[, Referent])
        Luni.fl <- length(uni.fl)
        for (uni in 1: Luni.fl) {
          flfl <- uni.fl[uni]
          for (temp.txZ in 1: no.group) {
            if (txZ[1, temp.txZ] == flfl) {
              txZ[1,temp.txZ] <- paste0("F", factor.no, "I", EP)
            }
          }
          EP <<- EP + 1
        }
        txX[1, (no.group+1)] <- Referent
        txZ[1, (no.group+1)] <- Referent

        Arg <- 1

        if (Referent > Arg) {
          txX[1, (no.group+2)] <- Arg
          txZ[1, (no.group+2)] <- Arg
        } else {
          txX[1, (no.group+2)] <- Arg + 1
          txZ[1, (no.group+2)] <- Arg + 1
        }  ## end if Referent
#$      print(txX, quote=FALSE)
#        txY <<- rbind(txY,txX)
#$      print(txZ, quote=FALSE)
#        txYY <<- rbind(txYY,txZ)

        next
       } # skip if unequal factor loadings


      ## ==  Print Factor Loadings FL of all groups  == ##
      FL <- matrix(1:no.k, nrow = no.group)  # factor loading matrix
      colnames(FL) <- paste0(' Item ', 1:no.items[factor.no])
      rownames(FL) <- c(paste0("Group ", 1:no.group))

      if (Referent == no.markers[factor.no]) {  ## if referent is the marker item
        for (i in 1:no.group) {
          for (j in 1:no.items[factor.no]) {
            if (j == Referent) {
              FL[i, j] <- 1
            } else {
              FL[i, j] <- FL.PMI[i, j]
            }  ## end if j
          }  ## end loop j
        }  ## end loop i
      } else {  ## referent is not the marker item
        for (i in 1:no.group) {
          for (j in 1:no.items[factor.no]) {
            if (j == Referent) {
              FL[i, j] <- 1
            } else {
              FL[i, j] <- FL.PMI[i, j]/FL.PMI[i, Referent]
            }  ## end if j
          }  ## end loop j
        }  ## end loop i
      }  ## end if referent

#$      cat(rep("\n", 3), "Factor Loadings ", rep("\n", 2))
#$      print(round(FL[], digits=4))

      ## == End (Print Factor Loadings)  == ##

      ## ==  Print Intercepts of all groups  == ##
      no.k <- no.items[factor.no]*no.group  # number of intercepts in TX
      TX <- matrix(1:no.k, nrow = no.group)  # intercept matrix
      colnames(TX) <- paste0(' Item ', 1:no.items[factor.no])
      rownames(TX) <- c(paste0("Group ", 1:no.group))

      if (Referent == no.markers[factor.no]) {  ## if referent is the marker item
        for (i in 1:no.group) {
          for (j in 1:no.items[factor.no]) {
            if (j == Referent) {
              TX[i, j] <- 0
            } else {
              TX[i, j] <- TX.PMI[i, j]
            }  ## end if j
          }  ## end loop j
        }  ## end loop i
      } else {  ## referent is not the marker item
        for (i in 1:no.group) {
          for (j in 1:no.items[factor.no]) {
            if (j == Referent) {
              TX[i, j] <- 0
            } else {
              TX[i, j] <- TX.PMI[i, j] - FL.PMI[i,j]/FL.PMI[i, Referent]*TX.PMI[i, Referent]
            }  ## end if j
          }  ## end loop j
        }  ## end loop i
      }  ## end if referent

#$      cat(rep("\n", 3), "Intercepts ", rep("\n", 2))
#$      print(round(TX[], digits=4))

      ## == End (Print Intercepts)  == ##


      for (nArg in 1: no.items[factor.no]) {  ## Argument item number for comparison
        no.dif <- factorial(no.group) / factorial(no.group - 2) / 2  # No. of pairwise comparisons for each item
        boot.dif.tx <- matrix(0, nrow(bootcoef), no.dif)  # Create bootstrap difference matrix
        samp.dif.tx <- matrix(0, no.group, no.group) # Create sample difference matrix
        if (nArg == Referent) {
          next
        } else {
          Arg <- nArg
        } ## end if nArg

        ## == Calculate bootstrap difference and sample estimate difference == ##
        comp = 0
        for (r in 1:(no.group-1)) {  ## r is referent group
          for (a in (r+1):no.group) {  ## a is argument group
            comp = comp + 1
            if (Referent == no.markers[factor.no]) {
                boot.dif.tx[,comp] <- bootcoef[, paste0(names.item[[factor.no]][[Arg]],"_G1", "~1")] -
                                      bootcoef[, paste0(names.item[[factor.no]][[Arg]],"_G2", "~1")]
                samp.dif.tx[r,a] <- par.est[paste0(names.item[[factor.no]][[Arg]],"_G1", "~1")] -
                                    par.est[paste0(names.item[[factor.no]][[Arg]],"_G2", "~1")]
#$ print(paste0("Referent = ", Referent, " & Arg = ", Arg))
#$ print(par.est[paste0(names.item[[factor.no]][[Arg]],"_G1", "~1")])
#$ print(par.est[paste0(names.item[[factor.no]][[Arg]],"_G2", "~1")])
            } else { # if (Referent > no.markers[factor.no]) {
              if (Arg == no.markers[factor.no]) {
                boot.dif.tx[, comp] <- (-1*bootcoef[, paste0(names.item[[factor.no]][[Referent]],"_G1", "~1")]/
                                           bootcoef[, paste0(names.lv[factor.no],"_G1=~", names.item[[factor.no]][[Referent]],"_G1")]) -
                                       (-1*bootcoef[, paste0(names.item[[factor.no]][[Referent]],"_G2", "~1")]/
                                           bootcoef[, paste0(names.lv[factor.no],"_G2=~", names.item[[factor.no]][[Referent]],"_G2")])
                samp.dif.tx[r,a] <- (-1*par.est[paste0(names.item[[factor.no]][[Referent]],"_G1", "~1")]/
                                        par.est[paste0(names.lv[factor.no],"_G1=~", names.item[[factor.no]][[Referent]],"_G1")]) -
                                    (-1*par.est[paste0(names.item[[factor.no]][[Referent]],"_G2", "~1")]/
                                        par.est[paste0(names.lv[factor.no],"_G2=~", names.item[[factor.no]][[Referent]],"_G2")])
#$ print(paste0("Referent = ", Referent, " & Arg = ", Arg))
#$ print(par.est[paste0(names.item[[factor.no]][[Referent]],"_G1", "~1")])
#$ print(par.est[paste0(names.lv[factor.no],"_G1=~", names.item[[factor.no]][[Referent]],"_G1")])
#$ print(par.est[paste0(names.item[[factor.no]][[Referent]],"_G2", "~1")])
#$ print(par.est[paste0(names.lv[factor.no],"_G2=~", names.item[[factor.no]][[Referent]],"_G2")])
              } else { ## (Arg != no.markers[factor.no])
                boot.dif.tx[,comp] <- (bootcoef[, paste0(names.item[[factor.no]][[Arg]],"_G1", "~1")]-
                                       bootcoef[, paste0(names.lv[factor.no],"_G1=~", names.item[[factor.no]][[Arg]],"_G1")]*
                                       bootcoef[, paste0(names.item[[factor.no]][[Referent]],"_G1", "~1")]/
                                       bootcoef[, paste0(names.lv[factor.no],"_G1=~", names.item[[factor.no]][[Referent]],"_G1")]) -
                                      (bootcoef[, paste0(names.item[[factor.no]][[Arg]],"_G2", "~1")]-
                                       bootcoef[, paste0(names.lv[factor.no],"_G2=~", names.item[[factor.no]][[Arg]],"_G2")]*
                                       bootcoef[, paste0(names.item[[factor.no]][[Referent]],"_G2", "~1")]/
                                       bootcoef[, paste0(names.lv[factor.no],"_G2=~", names.item[[factor.no]][[Referent]],"_G2")])
                samp.dif.tx[r,a] <- (par.est[paste0(names.item[[factor.no]][[Arg]],"_G1", "~1")]-
                                     par.est[paste0(names.lv[factor.no],"_G1=~", names.item[[factor.no]][[Arg]],"_G1")]*
                                     par.est[paste0(names.item[[factor.no]][[Referent]],"_G1", "~1")]/
                                     par.est[paste0(names.lv[factor.no],"_G1=~", names.item[[factor.no]][[Referent]],"_G1")]) -
                                    (par.est[paste0(names.item[[factor.no]][[Arg]],"_G2", "~1")]-
                                     par.est[paste0(names.lv[factor.no],"_G2=~", names.item[[factor.no]][[Arg]],"_G2")]*
                                     par.est[paste0(names.item[[factor.no]][[Referent]],"_G2", "~1")]/
                                     par.est[paste0(names.lv[factor.no],"_G2=~", names.item[[factor.no]][[Referent]],"_G2")])
#$ print(paste0("Referent = ", Referent, " & Arg = ", Arg))
#$ print(par.est[paste0(names.item[[factor.no]][[Arg]],"_G2", "~1")])
#$ print(par.est[paste0(names.lv[factor.no],"_G2=~", names.item[[factor.no]][[Arg]],"_G2")])
#$ print(par.est[paste0(names.item[[factor.no]][[Referent]],"_G2", "~1")])
#$ print(par.est[paste0(names.lv[factor.no],"_G2=~", names.item[[factor.no]][[Referent]],"_G2")])
              }  ## end if Arg
            }  ## end if Referent
          }  ## end loop a
        }  ## end loop r

#$ print(paste0("Referent = ", Referent, " Argument = ", nArg))
#$ print(samp.dif.tx)


        ## == Calculate Percentile Probability == ##

        F1.comp.pp <- matrix(0, no.group, no.group)  ## Percentile Probability
        pno_nipair <- 0  ## Set number of non-invariant pairs to zero

        comp = 0
        for (r in 1:(no.group-1)) {  ## r is the referent group
          for (a in (r+1):no.group) {  ## a is the argument
            comp = comp + 1
            if (quantile(boot.dif.tx[, comp], probs = 0.5, na.rm = TRUE) > 0) {
              F1.comp.pp[r,a] = 2*(sum(boot.dif.tx[, comp] < 0, na.rm=TRUE)/bootno)
            } else {
              F1.comp.pp[r,a] = 2*(sum(boot.dif.tx[, comp] > 0, na.rm=TRUE)/bootno)
            }  ## end if
          }  ## end loop a
        }  ## end loop r


        ## == Percentile Confidence Intervals of Pairwise Comparisons == ##

        PCI <- matrix(1:(no.dif*10), nrow = no.dif)
        colnames(PCI) <- c("Group 1","Group 2","0.5%","2.5%","5%","Estimate","95%","97.5%","99.5%","p-value")

        comp = 0
        for (r in 1:(no.group-1)) {
          for (a in (r+1):no.group) {
            comp = comp + 1
            PCI[comp, 1] <- r
            PCI[comp, 2] <- a
            PCI[comp, 3] <- quantile(boot.dif.tx[, comp],c(0.005), na.rm = TRUE)
            PCI[comp, 4] <- quantile(boot.dif.tx[, comp],c(0.025), na.rm = TRUE)
            PCI[comp, 5] <- quantile(boot.dif.tx[, comp],c(0.05), na.rm = TRUE)
            PCI[comp, 6] <- samp.dif.tx[r,a]
            PCI[comp, 7] <- quantile(boot.dif.tx[, comp],c(0.95), na.rm = TRUE)
            PCI[comp, 8] <- quantile(boot.dif.tx[, comp],c(0.975), na.rm = TRUE)
            PCI[comp, 9] <- quantile(boot.dif.tx[, comp],c(0.995), na.rm = TRUE)
            PCI[comp,10] <- F1.comp.pp[r,a]
            if (FL[r,Arg] != FL[a,Arg]) {
              if (pno_nipair == 0)	{
                pnipair <- c(r,a)
                pno_nipair <- pno_nipair + 1
              } else {
                pnipair <- c(pnipair,r,a)
                pno_nipair <- pno_nipair + 1
              }  ## end if pno_nipair
            } else {
              if (F1.comp.pp[r,a] < alpha) {
                if (pno_nipair == 0)	{
                  pnipair <- c(r,a)
                  pno_nipair <- pno_nipair + 1
                } else {
                  pnipair <- c(pnipair,r,a)
                  pno_nipair <- pno_nipair + 1
                }  ## end if pno_nipair
              }  ## end if F1.comp.pp
            }  ## end if FL[, Arg]
          }  ## end loop a
        }  ## end loop r


#$      cat("\n")
#$      if (Referent > Arg) {
#$        cat(paste0("Referent Item: ", Referent, "  Argument Item: ", Arg, "\n"))
#$      } else {
#$        cat(paste0("Referent Item: ", Referent, "  Argument Item: ", (Arg+1), "\n"))
#$      }  ## end if Referent
#$      cat("\n")
#$      print(round(PCI[], digits=4))

        count.tx <- count.tx + 1

        ## == Run List and Delete == ##
#$      cat(rep("\n", 2), "Sets of Groups with Invariant Intercepts", "\n")
        if (pno_nipair == no.dif) {

          txZ <<- matrix(" ", 1, (no.group+2))
          txX <<- matrix(1:(no.group+2), 1, (no.group+2))
          txX[1, 1:no.group] <- 0
          for (temp.txZ in 1: no.group) {
            txZ[1,temp.txZ] <- paste0("F", factor.no, "I", EP)
            EP <<- EP + 1
          }
          txX[1, (no.group+1)] <- Referent
          txZ[1, (no.group+1)] <- Referent
###          if (Referent > Arg) {
            txX[1, (no.group+2)] <- Arg
            txZ[1, (no.group+2)] <- Arg
###          } else {
###            txX[1, (no.group+2)] <- Arg + 1
###            txZ[1, (no.group+2)] <- Arg + 1
###          }  ## end if Referent
#$        print(txX, quote=FALSE)
          txY <<- rbind(txY,txX)
#$        print(txZ, quote=FALSE)
          txYY <<- rbind(txYY,txZ)

        } else if (pno_nipair > 0) {
#          txX <<- matrix(" ", 1, (no.group+2))
#          txZ <<- matrix(" ", 1, (no.group+2))

          listanddelete.tx(factor.no, no.group, no_nipair=pno_nipair, nipair=pnipair, Referent, Arg, txY, txYY, EP) ## Run list and delete tx

#$        print(txX, quote=FALSE)
#$        print(txZ, quote=FALSE)

        } else {  # (pno_nipair = 0)
          txX <<- matrix(1:(no.group+2), 1, (no.group+2))
          txZ <<- matrix(" ", 1, (no.group+2))
          txZ <- matrix(paste0("F", factor.no, "I", EP), 1, (no.group+2))
          EP <<- EP + 1
          txX[1, (no.group+1)] <- Referent
          txZ[1, (no.group+1)] <- Referent
###          if (Referent > Arg) {
            txX[1, (no.group+2)] <- Arg
            txZ[1, (no.group+2)] <- Arg
###          } else {
###            txX[1, (no.group+2)] <- Arg + 1
###            txZ[1, (no.group+2)] <- Arg + 1
###          }  ## end if Referent
#$        print(txX, quote=FALSE)
          txY <<- rbind(txY,txX)
#$        print(txZ, quote=FALSE)
          txYY <<- rbind(txYY,txZ)

        }  ## end if pno_nipair
      } ## end loop Arg
    } ## end loop Referent

    if (length(unique(txY[, no.group+1])) == 2 & no.items[factor.no] == 2) { # Only 2 groups with more than 2 items
      txX <<- matrix(1:(no.group+2), 1, (no.group+2))
      txZ <<- matrix(" ", 1, (no.group+2))
      txYa <<- txY[nrow(txY), no.group+1]
      txYb <<- txY[nrow(txY), no.group+2]
      txX[1, no.group + 1] <- txYb
      txX[1, no.group + 2] <- txYa
      txY <<- rbind(txY,txX)

      txZ[1, no.group + 1] <- txYb
      txZ[1, no.group + 2] <- txYa
      for (temp.txYY in 1: no.group) {
        txZ[1, temp.txYY] <- paste0("F", factor.no, "I", EP)
        EP <<- EP + 1
      }
      txYY <<- rbind(txYY,txZ)
    }

  txY <<- txY[-c(1),]
#$  cat("\n")
#$  print(txY, quote = FALSE)

  txYY <<- txYY[-c(1),]
#$  cat("\n")
#$  print(txYY, quote = FALSE)

  if (is.matrix(txY) == FALSE) { txY <<- matrix(txY, nrow=1) }
  if (is.matrix(txYY) == FALSE) { txYY <<- matrix(txYY, nrow=1) }

    NIcombine <- table(txY[, no.group+1], txY[, no.group+2])
#    NI.level <- unique(txY[, no.group+1])
#    NIcombine <- NIcombine[NI.level, NI.level]


    tempR <- unique(txYY[, no.group+1])

    for (aR in 1:length(tempR)) {
      R <- tempR[aR]
      for (r in 1:(no.items[factor.no])) {
        if (r == 1) {
          assign(paste0("Rset", r), c(length(which(txYY[,no.group+1] < R))+1:(NIcombine[R, r])))
        } else if (r == R) {
          assign(paste0("Rset", r), 0)
        } else {
          assign(paste0("Rset", r), c(length(which(txYY[,no.group+1] < R))+sum(NIcombine[R, 1:(r-1)])+1:(NIcombine[R, r])))
        }
        if (R == 1) { assign("Rset1", 0) }
      }  ## end loop r

      ## expand.grid(Rset1, Rset2, Rset3, Rset4, Rset5 ...)

      temp.row <- eval(parse(text = paste0("expand.grid(", paste0(c("Rset"), 1:no.items[factor.no], collapse=","),")")))
      if (aR == 1) {
        Set.row <- temp.row
      } else {
        Set.row <- rbind(Set.row, temp.row)
      } ## end if R
    }  ## end loop aR

    Model.load <- array(" ", dim = c(no.items[factor.no], no.group, nrow(Set.row)))  # all possible invariance models

    for (Xset in 1: nrow(Set.row)) {
      a <- matrix(" ", no.items[factor.no], no.group)
      for (r in 1:no.items[factor.no]) {
        if (Set.row[Xset, r] == 0) {
          a[r, ] <- 0
        } else {
          a[r, ] <- txYY[Set.row[Xset, r], 1:no.group]
        } ## end if Set.row
###@@@        if (txYY[r, no.group+1] > 1) { break }
      }  ## end loop r
      Model.load[, , Xset] <- a

    }  ## end loop Xset

#$ print(Model.load)

##########################

  ## == Model.summary == ##
  # Column 1 - number of estimated intercepts
  # Column 2 - number of measurement invariance items
  # Column 3 - p-value of Model chi-square
  # Column 4 - CFI
  # Column 5 - RMSEA
  # Column 6 - SRMR

    Model.summary <- matrix(0, nrow(Set.row), 7)
    for (R in 1:nrow(Set.row)) {
      Model.summary[R, 1] <- length(unique(c(Model.load[,,R]))) - 1
      no.MI <- no.items[factor.no] - 1
      for (Ra in 1: no.items[factor.no]) {
        if (length(unique(c(Model.load[Ra,,R]))) > 1) {no.MI <- no.MI - 1}
      }  ## end loop Ra
      Model.summary[R, 2] <- no.MI
    }  ## end loop R


    ## == Recommended Models (Rec.Model) - Largest no. of MI items, then smallest number of estimated loadings == ##
    Rec.Model <- which(Model.summary[,2] == max(Model.summary[,2]))
    Rec.Model <- Rec.Model[which(Model.summary[Rec.Model,1] == min(Model.summary[Rec.Model,1]))]


    ## == Find recommended model if no.items > 2 == ##
    if (no.items[factor.no] > 2) {

      Rec.Model.load <- matrix(nrow = length(Rec.Model), ncol = length(lavaan::coef(PMI.Model.fit)))
      MODFIT <- matrix(0, nrow = length(Rec.Model), ncol = 9)
      colnames(MODFIT) <-
        c("  chisq.scaled", "  df.scaled", "  pvalue.scaled", "  rmsea.scaled", "  cfi.scaled", "  tli.scaled", "  srmr_bentler", "  aic", "  bic")
      Rec.Model.FL <- array(0, dim = c(no.group, no.items[factor.no], length(Rec.Model)))
      colnames(Rec.Model.FL) <- paste0(" Item ", 1:no.items[factor.no])
      rownames(Rec.Model.FL) <- c(paste0("Group ", 1:no.group))

      Rec.Model.TX <- array(0, dim = c(no.group, no.items[factor.no], length(Rec.Model)))
      colnames(Rec.Model.TX) <- paste0(" Item ", 1:no.items[factor.no])
      rownames(Rec.Model.TX) <- c(paste0("Group ", 1:no.group))

      PSI.XX <- array(" ", length(Rec.Model))
      for (Model.R in 1: length(Rec.Model)) {

        ## == Specify the factor loadings for partial scalar invariance models (PSI.X) == ##

        PSI.marker <- which(Model.load[,1,Rec.Model[Model.R]] == 0)
        if (length(PSI.marker) > 1) {PSI.marker <<- PSI.marker[1]}
        PSI.FL <- matrix(1:no.k, nrow = no.group)  # factor loading matrix
        PSI.load <- matrix(" ", no.items[factor.no], no.group)  # factor loading matrix

        for (FL.item.g in 1:no.group) {
          for (FL.item in 1:no.items[factor.no]) {
            PSI.FL[FL.item.g, FL.item] <- FL.PMI[FL.item.g,FL.item]/FL.PMI[FL.item.g,PSI.marker]
          }  ## end loop FL.item
        }  ## end loop FL.item.g

        FL.count <- 1
        for (FL.item in 1:no.items[factor.no]) {
          if (PSI.FL[FL.item.g, FL.item] == 1) {
            PSI.load[FL.item,] <- 1
          } else {
            temp.fl <- unique(PSI.FL[,FL.item])
            for (temp.x in 1:length(temp.fl)) {
              for (FL.item.g in 1:no.group) {
                if (PSI.FL[FL.item.g,FL.item] == temp.fl[temp.x]) {
                  PSI.load[FL.item,FL.item.g] <- paste0("F",factor.no,"L",FL.count)
                }  ## end if
              }  ## end loop FL.item.g
              FL.count <- FL.count + 1
            } ## end loop temp.x
          } ## end if PSI.FL
        }  ## end loop FL.item


        ## == Specify the partial metric invariance models (PMI.X) == ##
        PMI.X.1 <- paste0("PSI.X <- '", "\n")
        for (j in 1:2) {
          PMI.X.1 <- paste0(PMI.X.1, "  ", names.lv[factor.no], "_G",j," =~ ", PSI.load[1, j], "*", names.item[[factor.no]][[1]], "_G",j," + ")
          for (i in 2: (no.items[factor.no]-1)) {
            PMI.X.1 <- paste0(PMI.X.1, PSI.load[i, j], "*", names.item[[factor.no]][[i]], "_G",j," + ")
          } ## end loop i
          PMI.X.1 <- paste0(PMI.X.1, PSI.load[(i+1), j], "*", names.item[[factor.no]][[i+1]], "_G",j, "\n")
        } # end for j

        ## -- Create Residual Covariance -- ##
        PMI.X.1 <- paste0(PMI.X.1, "\n")
        for (i in 1:no.items[factor.no]) {
          PMI.X.1 <- paste0(PMI.X.1, "  ", names.item[[factor.no]][[i]], "_G1 ~~ ", names.item[[factor.no]][[i]], "_G2", "\n")
        }

        ## == Specify the partial scalar invariance models (PSI.X) == ##

        final.intercept <- unique(as.vector(Model.load[,, Rec.Model[Model.R]]))
        for (j in 1:2) {
          for (i in 1:no.items[factor.no]) {
            if (Model.load[i,j,Rec.Model[Model.R]] != 0) { Model.load[i,j, Rec.Model[Model.R]] <- paste0("F", factor.no, "I", which(Model.load[i,j,Rec.Model[Model.R]] == final.intercept)) }
          } # end for i
        } # end for j

        PSI.XX[Model.R] <- paste0(PMI.X.1, "\n")
        for (j in 1:2) {
          for (i in 1:no.items[factor.no]) {
            PSI.XX[Model.R] <- paste0(PSI.XX[Model.R], "  ", names.item[[factor.no]][[i]], "_G", j, " ~ ", Model.load[i, j, Rec.Model[Model.R]], "*1", "\n")
          } ## end for i
        } # end for j

        PSI.XX[Model.R] <- paste0(PSI.XX[Model.R], "'", "\n")
        eval(parse(text = PSI.XX[Model.R]))

#$ cat(PSI.XX[Model.R])

        ## == Run PSI.X == ##
        if (Cluster == "NULL") {
          eval(parse(text=   "PSI.X.fit <- lavaan::sem(PSI.X,
                           data.source,
                           missing = 'fiml',
                           auto.fix.first = FALSE,
                           marker.int.zero = TRUE,
                           meanstructure = T,
                           information = 'observed',
                           estimator = 'MLR')"))
        } else {
          eval(parse(text=   "PSI.X.fit <- lavaan::sem(PSI.X,
                           data.source,
                           missing = 'fiml',
                           auto.fix.first = FALSE,
                           marker.int.zero = TRUE,
                           meanstructure = T,
                           information = 'observed',
                           cluster = Cluster,
                           estimator = 'MLR')"))
        }

        ## == Request summary outputs == ##
#$      eval(parse(text = "print(lavaan::summary(PSI.X.fit, fit.measure = T, standardized = T, rsq = T))"))


        ## == Save fit indices == ##
        if (correct.cfi == TRUE) {
          fitmeasures <- PAIRcorrected.cfi(PSI.object=PSI.X.fit, data=data.source, cluster=Cluster)
          XX <- fitmeasures[c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr_bentler", "aic", "bic")]
        } else {
          XX <- lavaan::fitMeasures(PSI.X.fit,
                                  c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr_bentler", "aic", "bic"))
        } # end (if correct.cfi)

        MODFIT[Model.R,] <- format(round((XX), digits = 4), nsmall = 4, scientific = FALSE)

        ## == Save Intercepts == ##
        temp.tx <- lavaan::lavInspect(PSI.X.fit, what='est')$nu
        for (G in 1: no.group) {
          Rec.Model.TX[G, , Model.R] <- temp.tx[(no.items[factor.no]*(G-1)+1):(no.items[factor.no]*G)]
        }  ## end loop G
      } ## end loop Model.R


      ## == Print model fit indices == ##
      rownames(MODFIT) <- c(paste0("Model ", 1: length(Rec.Model)))
#$    cat(rep("\n", 3))
#$    cat("####    MODEL FIT INDICES   ####", rep("\n", 2))
#$    print(MODFIT, quote=FALSE, right=TRUE)

      ## == Print factor loadings == ##
#$    for (Model.R in 1: length(Rec.Model)) {
#$      cat(rep("\n", 2), "Factor Loadings of Model", Model.R, rep("\n", 2))
#$      print(round(Rec.Model.FL[,,Model.R], digits=4))
#$    }  ## end loop Model.R


      ## == Save Recommended Model (Recommend.Model) -- First model with best model fit (R.Model) == ##
      class(MODFIT) <- "numeric"
      if (BMSC == "Chisquare") {
        Model.R <- which.min(MODFIT[,1])
      } else if (BMSC == "CFI") {
        Model.R <- which.max(MODFIT[,5])
      } else if (BMSC == "RMSEA") {
        Model.R <- which.min(MODFIT[,4])
      } else if (BMSC == "SRMR") {
        Model.R <- which.min(MODFIT[,7])
      }
      R.Model <- Model.R[1] # Select the first model if 2 or more models have the same fit
      if (factor.no == 1) {Recommend.Model <- matrix("PSI.Model.R <- '", 1)}  ## Reset Recommend.Model

      Recommend.Model <- rbind(Recommend.Model, PSI.XX[R.Model])

    }  ## end if no.items[factor.no] > 2


    if (no.items[factor.no] < 3) {   ## Begin Models with 2 items ##
        Rec.Model <- Rec.Model[1] # Select the first model if 2 or more models have the same fit
        PSI.marker <- which(Model.load[,1,Rec.Model] == 0)
        if (length(PSI.marker) > 1) {PSI.marker <<- PSI.marker[1]}
        PSI.FL <- matrix(1:no.k, nrow = no.group)  # factor loading matrix
        PSI.load <- matrix(" ", no.items[factor.no], no.group)  # factor loading matrix

        for (FL.item.g in 1:no.group) {
          for (FL.item in 1:no.items[factor.no]) {
            PSI.FL[FL.item.g, FL.item] <- FL.PMI[FL.item.g,FL.item]/FL.PMI[FL.item.g,PSI.marker]
          }  ## end loop FL.item
        }  ## end loop FL.item.g

        FL.count <- 1
        for (FL.item in 1:no.items[factor.no]) {
          if (PSI.FL[FL.item.g, FL.item] == 1) {
            PSI.load[FL.item,] <- 1
          } else {
            temp.fl <- unique(PSI.FL[,FL.item])
            for (temp.x in 1:length(temp.fl)) {
              for (FL.item.g in 1:no.group) {
                if (PSI.FL[FL.item.g,FL.item] == temp.fl[temp.x]) {
                  PSI.load[FL.item,FL.item.g] <- paste0("F",factor.no,"L",FL.count)
                }  ## end if
              }  ## end loop FL.item.g
              FL.count <- FL.count + 1
            } ## end loop temp.x
          } ## end if PSI.FL
        }  ## end loop FL.item

      R.Model <- Model.load[,,Rec.Model]  # Recommended Model

      Recommend.Model <- rbind(Recommend.Model, "\n")

      ## == Save Recommended Model (Recommend.Model) == ##
        if (factor.no == 1) {Recommend.Model <- matrix("PSI.Model.R <- '", 1)}

        for (j in 1:2) {
          Recommend.Model <- rbind(Recommend.Model, paste0("  ", names.lv[factor.no], "_G", j, " =~ ", PSI.load[1,j], "*", names.item[[factor.no]][[1]], "_G", j, " + "))
          Recommend.Model <- rbind(Recommend.Model, paste0(PSI.load[(2), j],"*", names.item[[factor.no]][[2]], "_G",j, "\n"))
        } # end for j

        ## -- Create Residual Covariance -- ##
        Recommend.Model <- rbind(Recommend.Model, "\n")
        for (i in 1:no.items[factor.no]) {
          Recommend.Model <- rbind(Recommend.Model, paste0("  ", names.item[[factor.no]][[i]], "_G1 ~~ ", names.item[[factor.no]][[i]], "_G2", "\n"))
        }

        ## == Specify the partial scalar invariance models (PSI.X) == ##

        final.intercept <- unique(as.vector(Model.load[,, Rec.Model[Model.R]]))
        for (j in 1:2) {
          for (i in 1:no.items[factor.no]) {
            if (Model.load[i,j,Rec.Model[Model.R]] != 0) { Model.load[i,j, Rec.Model[Model.R]] <- paste0("F", factor.no, "I", which(Model.load[i,j,Rec.Model[Model.R]] == final.intercept)) }
          } # end for i
        } # end for j

        Recommend.Model <- rbind(Recommend.Model, "\n")
        for (j in 1:(2)) {
          for (i in 1:no.items[factor.no]) {
            Recommend.Model <- rbind(Recommend.Model, paste0("  ", names.item[[factor.no]][[i]], "_G", j, " ~ ", Model.load[i, j, Rec.Model[Model.R]], "*1", "\n"))
          } # end for i
        } # end for j

    }  ## End (Models with 2 items)

  } ## End (loop factor.no)


  Recommend.Model[1] <- sub("''","'", Recommend.Model[1])
  for (i in 1:no.factor) {
    Recommend.Model[i+1] <-  sub("PSI.X <- '", " ", Recommend.Model[i+1])
    Recommend.Model[i+1] <- gsub("\n'", " ", Recommend.Model[i+1])
  } # end loop i
  Recommend.Model <- rbind(Recommend.Model, "  '", "\n")  ## last line of Recommended Model


  cat(rep("\n", 2), "## =====  Recommended Model  ===== ##", rep("\n", 2))
  for (i in 1: nrow(Recommend.Model)) {
    if (i == 1) {
      cat(Recommend.Model[i], "\n")
    } else {
      cat(Recommend.Model[i])
    }
  } ## end loop i

  cat(rep("\n", 2))
  cat("## ===== lavaan Outputs of PSI.Model.R ===== ##")
  cat(rep("\n", 2))

  ## == Start print recommended model to PSI.txt (Partial Metric Invariance Model) == ##
  sink('PSI.txt')

  cat("## =====  Recommended Model  ===== ##")
  cat("\n")
    for (i in 1: nrow(Recommend.Model)) {
      if (i == 1) {
        cat(Recommend.Model[i], "\n")
      } else {
        cat(Recommend.Model[i])
      }
    } ## end loop i


  cat(rep("\n",2), "## Run PSI.Model.R")
  cat(rep("\n",2), paste0("  PSI.Model.fit <- lavaan::sem(PSI.Model.R,"))
  cat("\n", paste0("    ", arg2_char, ","))
  cat("\n", "    missing = 'fiml',")
  cat("\n", "    auto.fix.first = FALSE,")
  cat("\n", "    marker.int.zero = TRUE,")
  cat("\n", "    meanstructure = T,")
  cat("\n", "    information = 'observed',")
  if (Cluster != "NULL") { cat("\n", paste0("    cluster = ", arg4_char,",")) }
  cat("\n", "    estimator = 'MLR')")

  cat(rep("\n", 3), "## Request summary outputs")
  cat(rep("\n", 2), paste0("  print(lavaan::summary(PSI.Model.fit, fit.measure = T, standardized = T, rsq = T))"), "\n")

  sink()  ## Stop writing to file

  PSI.Model.R <<- eval(parse(text=Recommend.Model))

  ## == Run PSI.Model.R == ##
  if (Cluster == "NULL") {
    eval(parse(text=   "PSI.Model.fit <- lavaan::sem(PSI.Model.R,
                       data.source,
                       missing = 'fiml',
                       auto.fix.first = FALSE,
                       marker.int.zero = TRUE,
                       meanstructure = T,
                       information = 'observed',
                       estimator = 'MLR')"))
  } else {
    eval(parse(text=   "PSI.Model.fit <- lavaan::sem(PSI.Model.R,
                       data.source,
                       missing = 'fiml',
                       auto.fix.first = FALSE,
                       marker.int.zero = TRUE,
                       meanstructure = T,
                       information = 'observed',
                       cluster = Cluster,
                       estimator = 'MLR')"))
  }

  print(lavaan::summary(PSI.Model.fit, fit.measure = T, standardized = T, rsq = T))

  cat(rep("\n", 2), "## Compare fit indices across models ##", "\n")
  DMODFIT <- matrix(0, nrow = 3, ncol = 9)
  colnames(DMODFIT) <-
        c("  chisq.scaled", "  df.scaled", "  pvalue.scaled", "  rmsea.scaled", "  cfi.scaled", "  tli.scaled", "  srmr", "  aic", "  bic")
  rownames(DMODFIT) <- c("PMI.Model.fit", "PSI.Model.fit", "PSI.Model.fit - PMI.Model.fit")
  XX1 <- lavaan::fitMeasures(PMI.Model.fit, fit.measures=
      c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr", "aic", "bic"))
  DMODFIT[1,] <- XX1

  if (correct.cfi == TRUE) {
    fitmeasures <- PAIRcorrected.cfi(PSI.object=PSI.Model.fit, data=data.source, cluster=Cluster)
    XX2 <- fitmeasures[c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr_bentler", "aic", "bic")]
  } else {
    XX2 <- lavaan::fitMeasures(PSI.Model.fit, fit.measures=
                             c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr", "aic", "bic"))
  } # end (if correct.cfi)

  DMODFIT[2,] <- XX2

  DMODFIT[3,1] <- DMODFIT[2,1] -DMODFIT[1,1]
  DMODFIT[3,2] <- DMODFIT[2,2] -DMODFIT[1,2]
  DMODFIT[3,3] <- pchisq(q=DMODFIT[3,1], df=DMODFIT[3,2], lower.tail=FALSE)
  DMODFIT[3,4] <- DMODFIT[2,4] -DMODFIT[1,4]
  DMODFIT[3,5] <- DMODFIT[2,5] -DMODFIT[1,5]
  DMODFIT[3,6] <- DMODFIT[2,6] -DMODFIT[1,6]
  DMODFIT[3,7] <- DMODFIT[2,7] -DMODFIT[1,7]
  DMODFIT[3,8] <- DMODFIT[2,8] -DMODFIT[1,8]
  DMODFIT[3,9] <- DMODFIT[2,9] -DMODFIT[1,9]
  print(round(DMODFIT, digits=4))


  ## == Print Factor Loadings == ##
  Final.Model.FL <- matrix(0, sum(no.items), 2)
  rownames(Final.Model.FL) <- paste0(" Item ", names.ind[1:sum(no.items)])
  colnames(Final.Model.FL) <- c(paste0("Group ", 1:no.group))

  Final.par.est <- lavaan::parameterEstimates(PSI.Model.fit)
  ext <- c(which(Final.par.est[,"op"] == "=~"))
  aa <- Final.par.est[ext, "est"]

  ## Extract factor loadings ##
  n <- 1
  for (j in 1: 2) {
    Final.Model.FL[1:no.items[1], j] <- aa[(no.items[1]*(j-1)+1):(no.items[1]*j)]
  }
  n <- n + no.items[1]*2

  if(no.factor > 1){
    for (i in 2:no.factor) {
      for (j in 1: 2) {
        Final.Model.FL[(sum(no.items[1:(i-1)])+1):sum(no.items[1:i]), j] <- aa[(n+no.items[i]*(j-1)):(n+(no.items[i]*j-1))]
      }
      n <- n + no.items[i]*2
    }
  }

  cat("\n")
  cat("## === Factor Loadings in Final Model === ##")
  cat("\n")
  print(round(Final.Model.FL, digits=4))
  ## == End (Print Factor Loadings) == ##



  ## == Print Intercepts == ##
  Final.Model.TX <- matrix(0, sum(no.items), no.group)
  rownames(Final.Model.TX) <- paste0(" Item ", names.ind[1:sum(no.items)])
  colnames(Final.Model.TX) <- c(paste0("Group ", 1:no.group))

  temp.tx <- lavaan::lavInspect(PSI.Model.fit, what='est')$nu
  for (G in 1: no.group) {
    Final.Model.TX[1:no.items[1], G] <- temp.tx[(no.items[1]*(G-1)+1):(no.items[1]*G)]
  }
  if (no.factor>1) {
    for (i in 2: no.factor) {
      for (G in 1: no.group) {
        Final.Model.TX[(sum(no.items[1:(i-1)])+1):(sum(no.items[1:(i-1)])+no.items[i]),G] <-
         temp.tx[(sum(no.items[1:(i-1)])*2+no.items[i]*(G-1)+1):(sum(no.items[1:(i-1)])*2+no.items[i]*G)]
      }
    }
  }
  cat("\n")
  cat("## === Intercepts in Final Model === ##")
  cat("\n")
  print(round(Final.Model.TX, digits=4))

  ## == End (Print Intercepts) == ##


  ## == Print Latent Means == ##
  Final.Model.LM <- matrix(0, no.factor, no.group)
  rownames(Final.Model.LM) <- paste0(" ", names.lv)
  colnames(Final.Model.LM) <- c(paste0("Group ", 1:no.group))

  temp.LM <- lavaan::lavInspect(PSI.Model.fit, what='est')$alpha

  for (G in 1: factor.no) { Final.Model.LM[G,] <- temp.LM[(no.group*(G-1)+1):(no.group*G)]}
  cat("\n")
  cat("## === Latent Means in Final Model === ##")
  cat("\n")
  print(round(Final.Model.LM, digits=4))
  ## == End (Print Latent Means) == ##


  ## ========== Compare Latent Means ========== ##

  ## Extraxt latent means ##
  par.est <- lavaan::coef(PSI.Model.fit)  # sample parameters
  temp <- lavaan::parameterEstimates(PSI.Model.fit, remove.nonfree=TRUE)
  simvcov <- lavInspect(PSI.Model.fit, what="vcov")

  parsed <- temp[temp$op == "=~",]
  names.lm <- unique(parsed[,"lhs"])

  ext <- c(which(temp$op %in% "~1" & temp$lhs %in% names.lm))

  par.est <- par.est[ext]
  simvcov <- simvcov[ext,ext]
  no.lm.g <- length(names.lv)  # number of estimated latent means per group #

  cat(rep("\n", 3), "## ======= COMPARISON OF LATENT MEANS ======= ##", rep("\n", 2))  ## print heading

  if (TYPE == "MonteCarlo") {

    mcmc <- MASS::mvrnorm(n=2000000, mu=par.est, Sigma=simvcov, tol = 1e-6)  # Run 2,000,000 simulations
    bootcoef <- mcmc
    bootno <- nrow(mcmc)  # No. of successful simulated samples
    cat(paste0("Number of Successful Simulated Samples = ", bootno, "\n"))

  } else { # TYPE=Bootstrap

    ## == Simplified bootstrapping model == ##
    PSI.boot <- lavaan::cfa(PSI.Model.R, data = data.source,
                  meanstructure = TRUE,
                  auto.fix.first = FALSE,
                  marker.int.zero = TRUE,
                  ordered = FALSE,
                  missing = "fiml", se = "none", test = "none", check.start = FALSE, check.post = FALSE, check.gradient = FALSE)

    ## == Bootstrapping == ##
    bootcoef <- lavaan::bootstrapLavaan(PSI.boot, R = b.no, ncpus = max(1L, parallel::detectCores() - 1L), parallel="snow")

    # Remove error and nonadmissible bootstrap samples #
    B.na <- attr(bootcoef,"nonadmissible")
    B.er <- attr(bootcoef,"error.idx")
    B.de <- c(B.na,B.er)
    if (length(B.de) != 0) {
      B.re <- bootcoef[-c(B.de),]
      bootcoef <- B.re
    }
    bootno <- nrow(bootcoef)  # number of successful bootstrap samples
    cat(paste0("Number of Successful Bootstrapped Samples = ", bootno, "\n"))
    bootcoef <- bootcoef[,ext]  ## Extract the bootstrapped intercepts

  } ## end MonteCarlo or Bootstrap


  ## == Start the factor.no loop for CompareMeans == ##
  for (factor.no in 1: no.factor) {

    boot.dif.lm <- matrix(0, nrow(bootcoef), no.dif)  # Create bootstrap difference matrix
    samp.dif.lm <- matrix(0, no.group, no.group) # Create sample difference matrix

    ## == Calculate bootstrap difference and sample estimate difference == ##
    comp = 0
    for (r in 1:(no.group-1)) {  ## r is referent group
      for (a in (r+1):no.group) {  ## a is argument group
        comp = comp + 1
        boot.dif.lm[,comp] <- bootcoef[, paste0(names.lv[factor.no], "_G", r, "~1")] - bootcoef[,paste0(names.lv[factor.no], "_G", a, "~1")]
        samp.dif.lm[r,a] <- par.est[paste0(names.lv[factor.no], "_G", a, "~1")] - par.est[paste0(names.lv[factor.no], "_G", r, "~1")]
        samp.dif.lm[a,r] <- par.est[paste0(names.lv[factor.no], "_G", r, "~1")] - par.est[paste0(names.lv[factor.no], "_G", a, "~1")]
      }  ## end loop a
    }  ## end loop r

    colnames(samp.dif.lm) <- c(paste0("Group ", 1:no.group))
    rownames(samp.dif.lm) <- c(paste0("Group ", 1:no.group))


    ## == Calculate Number of Invariant Intercepts == ##
    LM.inv.tau <- matrix(0, no.group, no.group)
    for (r in 1:(no.group-1)) {  ## r is the referent group
      for (a in (r+1):no.group) {  ## a is the argument
        if (factor.no == 1) {
          LM.inv.tau[r,a] <- sum(round(Final.Model.TX[1:no.items[factor.no],r],6) == round(Final.Model.TX[1:no.items[factor.no],a],6)) 
        } else {
          LM.inv.tau[r,a] <- sum(round(Final.Model.TX[(sum(no.items[1:(factor.no-1)])+1):sum(no.items[1:factor.no]),r], 6) == 
                                 round(Final.Model.TX[(sum(no.items[1:(factor.no-1)])+1):sum(no.items[1:factor.no]),a], 6)) 
        }
        LM.inv.tau[a,r] <- LM.inv.tau[r,a]
      } # end for a
    } # end for r 


    ## == Calculate Percentile Probability == ##
    LM.comp.pp <- matrix(0, no.group, no.group)  ## Percentile Probability
    colnames(LM.comp.pp) <- c(paste0("Group ", 1:no.group))
    rownames(LM.comp.pp) <- c(paste0("Group ", 1:no.group))

    comp = 0
    for (r in 1:(no.group-1)) {  ## r is the referent group
      for (a in (r+1):no.group) {  ## a is the argument
        comp = comp + 1
        if (quantile(boot.dif.lm[, comp], probs = 0.5, na.rm = TRUE) > 0) {
          LM.comp.pp[r,a] = 2*(sum(boot.dif.lm[, comp] < 0, na.rm=TRUE)/bootno)
          LM.comp.pp[a,r] = 2*(sum(boot.dif.lm[, comp] < 0, na.rm=TRUE)/bootno)
        } else {
          LM.comp.pp[r,a] = 2*(sum(boot.dif.lm[, comp] > 0, na.rm=TRUE)/bootno)
          LM.comp.pp[a,r] = 2*(sum(boot.dif.lm[, comp] > 0, na.rm=TRUE)/bootno)
        }  ## end if
      }  ## end loop a
    }  ## end loop r

    ## == Percentile Confidence Intervals of Pairwise Comparisons == ##

    PCI <- matrix(1:(no.dif*10), nrow = no.dif)
    colnames(PCI) <- c("Group 1","Group 2","0.5%","2.5%","5%","Estimate","95%","97.5%","99.5%","p-value")

    comp = 0
    for (r in 1:(no.group-1)) {
      for (a in (r+1):no.group) {
        comp = comp + 1
        PCI[comp, 1] <- paste0("Group ", r)
        PCI[comp, 2] <- paste0("Group ", a)
        PCI[comp, 3] <- round(quantile(boot.dif.lm[, comp],c(0.005), na.rm = TRUE), digits=4)
        PCI[comp, 4] <- round(quantile(boot.dif.lm[, comp],c(0.025), na.rm = TRUE), digits=4)
        PCI[comp, 5] <- round(quantile(boot.dif.lm[, comp],c(0.05), na.rm = TRUE), digits=4)
        PCI[comp, 6] <- round(samp.dif.lm[r,a], digits=4)
        PCI[comp, 7] <- round(quantile(boot.dif.lm[, comp],c(0.95), na.rm = TRUE), digits=4)
        PCI[comp, 8] <- round(quantile(boot.dif.lm[, comp],c(0.975), na.rm = TRUE), digits=4)
        PCI[comp, 9] <- round(quantile(boot.dif.lm[, comp],c(0.995), na.rm = TRUE), digits=4)
        PCI[comp,10] <- round(LM.comp.pp[r,a], digits=4)
      }  ## end loop a
    }  ## end loop r

    cat("\n")
#$  print(PCI[],quote=F, nsmall=4, scientific=FALSE)


    # == Combine LM.comp.pp with LM.inv.tau == #
    LM.comp.pp <- matrix(sprintf("%.4f", LM.comp.pp), nrow = nrow(LM.comp.pp))
    LM.inv.tau <- matrix(sprintf("%.0f", LM.inv.tau), nrow = nrow(LM.inv.tau))
    LM.comp.pp <- matrix(paste0(LM.comp.pp, " (", LM.inv.tau, ")  "), nrow = nrow(LM.comp.pp), ncol = ncol(LM.comp.pp))
    colnames(LM.comp.pp) <- c(paste0("Group ", 1:no.group))
    rownames(LM.comp.pp) <- c(paste0("Group ", 1:no.group))
    diag(LM.comp.pp) <- ""

    samp.dif.lm <- matrix(sprintf("%.4f", samp.dif.lm), nrow = nrow(samp.dif.lm))
    diag(samp.dif.lm) <- ""
    colnames(samp.dif.lm) <- c(paste0("Group ", 1:no.group))
    rownames(samp.dif.lm) <- c(paste0("Group ", 1:no.group))


    cat("\n")
    cat("## ===== Latent Variable: ", names.lv[factor.no], " ===== ##")
    cat("\n")
    cat("== Factor Loadings ==")
    cat("\n")
    if (factor.no == 1) {
      print(format(round(Final.Model.FL[1:no.items[factor.no],],4)), quote=F)
    } else {
      print(format(round(Final.Model.FL[(sum(no.items[1:(factor.no-1)])+1):sum(no.items[1:factor.no]),],4)), quote=F)
    }
    cat("\n")
    cat("== Intercepts ==")
    cat("\n")
    if (factor.no == 1) {
      print(format(round(Final.Model.TX[1:no.items[factor.no],],4)), quote=F)
    } else {
      print(format(round(Final.Model.TX[(sum(no.items[1:(factor.no-1)])+1):sum(no.items[1:factor.no]),],4)), quote=F)
    }
    cat("\n")
    cat("== Latent Means ==")
    cat("\n")
    print(Final.Model.LM[factor.no,], quote=F, digits=5, nsmall=4, scientific=FALSE)
    cat("\n")
    cat("== Pairwise Difference in Latent Means ==")
    cat("\n")
    print(format(samp.dif.lm, justify="right"), quote=F)
    cat("\n")
    cat("== p-values for pairwise comparisons ==")
    cat("\n")
    print(format(LM.comp.pp, justify="right"), quote=F)
    cat("\n")
    cat("Note:", "\n")
    cat("Numbers in parentheses are numbers of items with invariant intercepts.", "\n")
    cat("There are ", no.items[factor.no], " items for this factor. Latent means should only compare with at least ", ceiling((no.items[factor.no]+1)/2), " items with invariant intercepts.", "\n")
  } # End (factor.no loop)

  cat(rep("\n",2),"The recommended model PSI.Model.R is saved in the file 'PSI.txt'", "\n")

} ## End (Function PAIRCompareMeans)

# ==================== Finish Function "PAIRCompareMeans" ==================== #





## ==========  Create Function List and Delete for LX ========== ##

listanddelete.lx <- function(factor.no = 1, no.group = 3, no_nipair = 2, nipair = c(1,3,2,3), Referent=1, Arg=2) {

#  no.group <- 3
#  no_nipair <- 1
#  nipair <- c(1,3,2,3,1,4,2,4,3,4,1,5,1,6,5,6,4,6,2,7,5,7,6,7,2,8,7,8,3,8,4,8)
#  nipair <- c(1,3)

#  flY <<- matrix(" ",1, (no.group+2))
#  flYY <<- matrix(" ",1, (no.group+2))
  nEP <- EP
#  factor.no <- 1
#  Referent <- 1
#  Arg <- 2

  MI <- t(matrix(nipair,no_nipair,nrow=2))
  MIset <- no.group-2
  mm3 <- matrix(0, 1, no.group+2)  ## Initialize mm3 (Matrix with all invariant sets) ##

  for (k in 1:MIset) {
    Noset <- no.group - k
    NIset <- factorial(no.group)/(factorial(no.group-k)*factorial(k))
    mm <- t(combn(no.group, Noset))  ## All combinations with Noset items
    for (j in 1:no_nipair) {
      a <- MI[j,1]
      b <- MI[j,2]
      for (i in 1:NIset) {
        count4 <- length(which(mm[i,] == a | mm[i,] == b))
        if (count4 > 1) { mm[i,] <- 0 }
      }  ## end for i
    }  ## end for j

    mm2 <- matrix(0, nrow(mm), no.group+2)  ## Temporary mm2
    for (i in 1:nrow(mm)) {
      for (j in 1:Noset) {
        if (mm[i, j] > 0) {mm2[i, mm[i,j]] = mm[i,j] }
      } ## end j
    }  ## end i
    mm3 <- rbind(mm3, mm2)
  }  ## end for k

  for (i in 1: nrow(mm3)) { mm3[i, no.group+2] <- length(which(mm3[i,] != 0)) }
  mm3 <- subset(mm3, mm3[, no.group+2] > 0)

  mm8 <- matrix(0, 1, no.group+2)  ## Initialize mm8 (Matrix with all combinations) ##

  if (nrow(mm3) == 0) {  ## no invariant set ##

    mmy <- matrix(0, 1, no.group+2)
    mm8 <- rbind(mm8, mmy)

  } else if (nrow(mm3) == 1) {  ## only 1 invariant set ##

    mmy <- matrix(0, 1, no.group+2)
    mmy[1] <- 1
    mm8 <- rbind(mm8, mmy)

  } else {  ## more than 1 invariant sets ##

    mm4 <- mm3  ## Non-redundant mm3 ##
    for (i in 1: (nrow(mm3)-1)) {
      for (j in (i+1):nrow(mm3)) {
        if (all(mm3[j, 1:no.group] %in% mm4[i, 1:no.group]) == TRUE) { mm4[j, no.group+1] <- 1 }
      }  ## end for j
    }  ## end for i
    for (i in 1:nrow(mm4)) { mm4[i, no.group+2] <- length(which(mm4[i, 1:no.group] != 0)) } ## end for i
    for (i in 1:nrow(mm4)) { if (mm4[i,no.group+2] < max(mm4[, no.group+2])) { mm4[i, no.group+1] <- 1 } }  ## end for i
    mm3[, no.group+2] <- 0
    mm4[, no.group+2] <- 0
    mm4 <- subset(mm4, mm4[, no.group+1] == 0)

    ## more than 1 empty cell ##

    for (Set.a in 1:nrow(mm4)) {
      mm6 <- mm4[Set.a,]  ## Set for comparison ##

      if (length(which(mm6[1:no.group] == 0)) == 1) {
        mmy <- matrix(0, 1, no.group+2)
        mmy[1] <- Set.a
        mm8 <- rbind(mm8, mmy)
        next
      } ## end if only 1 empty cell

      for (Set.b in (Set.a+1):nrow(mm3)) {
        if (Set.b > nrow(mm3)) { break }
        if (all(mm3[Set.b, mm6!=0]==0) == FALSE) {  ## without overlap
          next
        } else {
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          if (length(which(mm6[1:no.group] == 0)) > 1) {
            mm6[which(mm3[Set.b, 1:no.group] != 0)] <- which(mm3[Set.b,1:no.group] != 0)

      for (Set.c in (Set.b+1):nrow(mm3)) {
        if (Set.c > nrow(mm3)) { break }
        if (all(mm3[Set.c, mm6!=0]==0) == FALSE) { ## without overlap
          next
        } else { # find overlap
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          mmy[3] <- Set.c
          if (length(which(mm6[1:no.group] == 0)) > 1) {
            mm6[which(mm3[Set.c, 1:no.group] != 0)] <- which(mm3[Set.c,1:no.group] != 0)

      for (Set.d in (Set.c+1):nrow(mm3)) {
        if (Set.d > nrow(mm3)) { break }
        if (all(mm3[Set.d, mm6!=0]==0) == FALSE) { ## without overlap
          next
        } else { # find overlap
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          mmy[3] <- Set.c
          mmy[4] <- Set.d
          if (length(which(mm6[1:no.group] == 0)) > 1) {
            mm6[which(mm3[Set.d, 1:no.group] != 0)] <- which(mm3[Set.d,1:no.group] != 0)

      for (Set.e in (Set.d+1):nrow(mm3)) {
        if (Set.e > nrow(mm3)) { break }
        if (all(mm3[Set.e, mm6!=0]==0) == FALSE) { ## without overlap
          next
        } else { # find overlap
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          mmy[3] <- Set.c
          mmy[4] <- Set.d
          mmy[5] <- Set.e
          if (length(which(mm6[1:no.group] == 0)) > 1) {
            mm6[which(mm3[Set.e, 1:no.group] != 0)] <- which(mm3[Set.e,1:no.group] != 0)

      for (Set.f in (Set.e+1):nrow(mm3)) {
        if (Set.f > nrow(mm3)) { break }
        if (all(mm3[Set.f, mm6!=0]==0) == FALSE) { ## without overlap
          next
        } else { # find overlap
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          mmy[3] <- Set.c
          mmy[4] <- Set.d
          mmy[5] <- Set.e
          mmy[6] <- Set.f
          if (length(which(mm6[1:no.group] == 0)) > 1) {
            mm6[which(mm3[Set.f, 1:no.group] != 0)] <- which(mm3[Set.f,1:no.group] != 0)

      for (Set.g in (Set.f+1):nrow(mm3)) {
        if (Set.g > nrow(mm3)) { break }
        if (all(mm3[Set.g, mm6!=0]==0) == FALSE) { ## without overlap
          next
        } else { # find overlap
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          mmy[3] <- Set.c
          mmy[4] <- Set.d
          mmy[5] <- Set.e
          mmy[6] <- Set.f
          mmy[7] <- Set.g
          if (length(which(mm6[1:no.group] == 0)) > 1) {
            mm6[which(mm3[Set.g, 1:no.group] != 0)] <- which(mm3[Set.g,1:no.group] != 0)

      for (Set.h in (Set.g+1):nrow(mm3)) {
        if (Set.h > nrow(mm3)) { break }
        if (all(mm3[Set.h, mm6!=0]==0) == FALSE) { ## without overlap
          next
        } else { # find overlap
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          mmy[3] <- Set.c
          mmy[4] <- Set.d
          mmy[5] <- Set.e
          mmy[6] <- Set.f
          mmy[7] <- Set.g
          mmy[8] <- Set.h
          if (length(which(mm6[1:no.group] == 0)) > 1) {
            mm6[which(mm3[Set.h, 1:no.group] != 0)] <- which(mm3[Set.h,1:no.group] != 0)

      for (Set.i in (Set.h+1):nrow(mm3)) {
        if (Set.i > nrow(mm3)) { break }
        if (all(mm3[Set.i, mm6!=0]==0) == FALSE) { ## without overlap
          next
        } else { # find overlap
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          mmy[3] <- Set.c
          mmy[4] <- Set.d
          mmy[5] <- Set.e
          mmy[6] <- Set.f
          mmy[7] <- Set.g
          mmy[8] <- Set.h
          mmy[9] <- Set.i
          if (length(which(mm6[1:no.group] == 0)) > 1) {
            mm6[which(mm3[Set.i, 1:no.group] != 0)] <- which(mm3[Set.i,1:no.group] != 0)

      for (Set.j in (Set.i+1):nrow(mm3)) {
        if (Set.j > nrow(mm3)) { break }
        if (all(mm3[Set.j, mm6!=0]==0) == FALSE) { ## without overlap
          next
        } else { # find overlap
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          mmy[3] <- Set.c
          mmy[4] <- Set.d
          mmy[5] <- Set.e
          mmy[6] <- Set.f
          mmy[7] <- Set.g
          mmy[8] <- Set.h
          mmy[9] <- Set.i
          mmy[10] <- Set.j
          if (length(which(mm6[1:no.group] == 0)) > 1) {
            mm6[which(mm3[Set.j, 1:no.group] != 0)] <- which(mm3[Set.j,1:no.group] != 0)

      for (Set.k in (Set.j+1):nrow(mm3)) {
        if (Set.k > nrow(mm3)) { break }
        if (all(mm3[Set.k, mm6!=0]==0) == FALSE) { ## without overlap
          next
        } else { # find overlap
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          mmy[3] <- Set.c
          mmy[4] <- Set.d
          mmy[5] <- Set.e
          mmy[6] <- Set.f
          mmy[7] <- Set.g
          mmy[8] <- Set.h
          mmy[9] <- Set.i
          mmy[10] <- Set.j
          mmy[11] <- Set.k
          if (length(which(mm6[1:no.group] == 0)) > 1) {
            mm6[which(mm3[Set.k, 1:no.group] != 0)] <- which(mm3[Set.k,1:no.group] != 0)

      for (Set.l in (Set.k+1):nrow(mm3)) {
        if (Set.l > nrow(mm3)) { break }
        if (all(mm3[Set.l, mm6!=0]==0) == FALSE) { ## without overlap
          next
        } else { # find overlap
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          mmy[3] <- Set.c
          mmy[4] <- Set.d
          mmy[5] <- Set.e
          mmy[6] <- Set.f
          mmy[7] <- Set.g
          mmy[8] <- Set.h
          mmy[9] <- Set.i
          mmy[10] <- Set.j
          mmy[11] <- Set.k
          mmy[12] <- Set.l
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.l
    }  # if length Set.k
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.k
    }  # if length Set.j
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.j
    }  # if length Set.i
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.i
    }  # if length Set.h
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.h
    }  # if length Set.g
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.g
    }  # if length Set.f
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.f
    }  # if length Set.e
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.e
    }  # if length Set.d
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.d
    }  # if length Set.c
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.c
    } # if length Set.b
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.b

      if (Set.b == nrow(mm3)) {
        if (all(mm3[Set.b, mm6!=0]==0) == FALSE) {
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mm8 <- rbind(mm8, mmy)
        }
      }
    } # for Set.a

    for (i in 1: (nrow(mm8)-1)) {
      for (j in (i+1):nrow(mm8)) {
        if (all(mm8[j, 1:no.group] %in% mm8[i, 1:no.group]) == TRUE) { mm8[j, no.group+2] <- 1 }
      }  ## end for j
    }  ## end for i

  } ## end if no invariant groups


  flZ <- matrix(0, 1, no.group+2)

  mm8 <- subset(mm8, mm8[, no.group+2] != 1)
  mm8 <- mm8[2:nrow(mm8),]

  if (is.matrix(mm8) == FALSE) {  ## mm8 has only 1 row
    mm7 <- matrix(0, 1, no.group+2)
    for (i in 1: length(which(mm8[1:no.group] > 0))) {
      xx <- which(mm3[mm8[i], 1:no.group] != 0)
      for (j in 1:length(xx)) { mm7[xx[j]] <- paste0("F", factor.no, "L", nEP) } ## end for j
      if (length(xx != 0)) { nEP <- nEP + 1 }
      xx <- which(mm7[1:no.group] == 0)
      for (k in 1:length(xx)) { mm7[xx[k]] <- paste0("F", factor.no, "L", (nEP-1+k)) } ## end for k
      nEP <- nEP + k
    }  ## end for i
    flZ <- mm7
    flX <- mm7
    flX[which(flZ[1:no.group] != 0)] <- which(flZ[1:no.group] != 0)

  } else {
    for (i in 1: nrow(mm8)) {
      mm7 <- matrix(0, 1, no.group+2)
      for (j in 1: length(which(mm8[1, 1:no.group] > 0))) {
        xx <- which(mm3[mm8[i,j], 1:no.group] != 0)
        for (k in 1:length(xx)) { mm7[xx[k]] <- paste0("F", factor.no, "L", nEP) } ## end for k
        nEP <- nEP + 1
      } ## end for j
      xx <- which(mm7[1:no.group] == 0)
      for (j in 1:length(xx)) { mm7[xx[j]] <- paste0("F", factor.no, "L", (nEP-1+j)) } ## end for j
      nEP <- nEP + j
      flZ <- rbind(flZ, mm7)
    } ## end for i

    flZ <- flZ[2:nrow(flZ),]
    flX <- flZ
    for (i in 1:nrow(flZ)) { flX[i, which(flZ[i,1:no.group] != 0)] <- which(flZ[i,1:no.group] != 0) }
    for (i in 1: nrow(flZ)) { flZ[i, no.group+2] <- length(unique(flZ[i, 1:no.group])) }
    flX <- subset(flX, flZ[, no.group+2] == min(flZ[, no.group+2]))
    flZ <- subset(flZ, flZ[, no.group+2] == min(flZ[, no.group+2]))
  }


  flX[1:nrow(flX), (no.group+1)] <- Referent
  flZ[1:nrow(flZ), (no.group+1)] <- Referent

  flX[1:nrow(flX), (no.group+2)] <- Arg
  flZ[1:nrow(flZ), (no.group+2)] <- Arg

  flY <<- rbind(flY,flX)
  flYY <<- rbind(flYY,flZ)

  #$ print(flY)
  #$ print(flYY)
  EP <<- nEP

}  ## end function listanddelete.lx

## ==========  Finish Function List and Delete for LX ========== ##




## ==========  Create Function List and Delete for TX ========== ##

listanddelete.tx <- function(factor.no = 1, no.group = 3, no_nipair = 2, nipair = c(1,3,2,3), Referent=1, Arg=2, txY, txYY, EP) {


#  no.group <- 3
#  no_nipair <- 1
#  nipair <- c(1,3,2,3,1,4,2,4,3,4,1,5,1,6,5,6,4,6,2,7,5,7,6,7,2,8,7,8,3,8,4,8)
#  nipair <- c(1,3)
#  txY <<- matrix(" ",1, (no.group+2))
#  txYY <<- matrix(" ",1, (no.group+2))
#  factor.no <- 1
#  Referent <- 1
#  Arg <- 2

  txX <<- matrix(" ", 1, (no.group+2))
  txZ <<- matrix(" ", 1, (no.group+2))
  nEP <- EP

  MI <- t(matrix(nipair,no_nipair,nrow=2))
  MIset <- no.group-2
  mm3 <- matrix(0, 1, no.group+2)  ## Initialize mm3 (Matrix with all invariant sets) ##

  for (k in 1:MIset) {
    Noset <- no.group - k
    NIset <- factorial(no.group)/(factorial(no.group-k)*factorial(k))
    mm <- t(combn(no.group, Noset))  ## All combinations with Noset items
    for (j in 1:no_nipair) {
      a <- MI[j,1]
      b <- MI[j,2]
      for (i in 1:NIset) {
        count4 <- length(which(mm[i,] == a | mm[i,] == b))
        if (count4 > 1) { mm[i,] <- 0 }
      }  ## end for i
    }  ## end for j

    mm2 <- matrix(0, nrow(mm), no.group+2)  ## Temporary mm2
    for (i in 1:nrow(mm)) {
      for (j in 1:Noset) {
        if (mm[i, j] > 0) {mm2[i, mm[i,j]] = mm[i,j] }
      } ## end j
    }  ## end i
    mm3 <- rbind(mm3, mm2)
  }  ## end for k

  for (i in 1: nrow(mm3)) { mm3[i, no.group+2] <- length(which(mm3[i,] != 0)) }
  mm3 <- subset(mm3, mm3[, no.group+2] > 0)

  mm8 <- matrix(0, 1, no.group+2)  ## Initialize mm8 (Matrix with all combinations) ##

  if (nrow(mm3) == 0) {  ## no invariant set ##

    mmy <- matrix(0, 1, no.group+2)
    mm8 <- rbind(mm8, mmy)

  } else if (nrow(mm3) == 1) {  ## only 1 invariant set ##

    mmy <- matrix(0, 1, no.group+2)
    mmy[1] <- 1
    mm8 <- rbind(mm8, mmy)

  } else {  ## more than 1 invariant sets ##

    mm4 <- mm3  ## Non-redundant mm3 ##
    for (i in 1: (nrow(mm3)-1)) {
      for (j in (i+1):nrow(mm3)) {
        if (all(mm3[j, 1:no.group] %in% mm4[i, 1:no.group]) == TRUE) { mm4[j, no.group+1] <- 1 }
      }  ## end for j
    }  ## end for i
    for (i in 1:nrow(mm4)) { mm4[i, no.group+2] <- length(which(mm4[i, 1:no.group] != 0)) } ## end for i
    for (i in 1:nrow(mm4)) { if (mm4[i,no.group+2] < max(mm4[, no.group+2])) { mm4[i, no.group+1] <- 1 } }  ## end for i
    mm3[, no.group+2] <- 0
    mm4[, no.group+2] <- 0
    mm4 <- subset(mm4, mm4[, no.group+1] == 0)

    ## more than 1 empty cell ##

    for (Set.a in 1:nrow(mm4)) {
      mm6 <- mm4[Set.a,]  ## Set for comparison ##

      if (length(which(mm6[1:no.group] == 0)) == 1) {
        mmy <- matrix(0, 1, no.group+2)
        mmy[1] <- Set.a
        mm8 <- rbind(mm8, mmy)
        next
      } ## end if only 1 empty cell

      for (Set.b in (Set.a+1):nrow(mm3)) {
        if (Set.b > nrow(mm3)) { break }
        if (all(mm3[Set.b, mm6!=0]==0) == FALSE) {  ## without overlap
          next
        } else {
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          if (length(which(mm6[1:no.group] == 0)) > 1) {
            mm6[which(mm3[Set.b, 1:no.group] != 0)] <- which(mm3[Set.b,1:no.group] != 0)

      for (Set.c in (Set.b+1):nrow(mm3)) {
        if (Set.c > nrow(mm3)) { break }
        if (all(mm3[Set.c, mm6!=0]==0) == FALSE) { ## without overlap
          next
        } else { # find overlap
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          mmy[3] <- Set.c
          if (length(which(mm6[1:no.group] == 0)) > 1) {
            mm6[which(mm3[Set.c, 1:no.group] != 0)] <- which(mm3[Set.c,1:no.group] != 0)

      for (Set.d in (Set.c+1):nrow(mm3)) {
        if (Set.d > nrow(mm3)) { break }
        if (all(mm3[Set.d, mm6!=0]==0) == FALSE) { ## without overlap
          next
        } else { # find overlap
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          mmy[3] <- Set.c
          mmy[4] <- Set.d
          if (length(which(mm6[1:no.group] == 0)) > 1) {
            mm6[which(mm3[Set.d, 1:no.group] != 0)] <- which(mm3[Set.d,1:no.group] != 0)

      for (Set.e in (Set.d+1):nrow(mm3)) {
        if (Set.e > nrow(mm3)) { break }
        if (all(mm3[Set.e, mm6!=0]==0) == FALSE) { ## without overlap
          next
        } else { # find overlap
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          mmy[3] <- Set.c
          mmy[4] <- Set.d
          mmy[5] <- Set.e
          if (length(which(mm6[1:no.group] == 0)) > 1) {
            mm6[which(mm3[Set.e, 1:no.group] != 0)] <- which(mm3[Set.e,1:no.group] != 0)

      for (Set.f in (Set.e+1):nrow(mm3)) {
        if (Set.f > nrow(mm3)) { break }
        if (all(mm3[Set.f, mm6!=0]==0) == FALSE) { ## without overlap
          next
        } else { # find overlap
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          mmy[3] <- Set.c
          mmy[4] <- Set.d
          mmy[5] <- Set.e
          mmy[6] <- Set.f
          if (length(which(mm6[1:no.group] == 0)) > 1) {
            mm6[which(mm3[Set.f, 1:no.group] != 0)] <- which(mm3[Set.f,1:no.group] != 0)

      for (Set.g in (Set.f+1):nrow(mm3)) {
        if (Set.g > nrow(mm3)) { break }
        if (all(mm3[Set.g, mm6!=0]==0) == FALSE) { ## without overlap
          next
        } else { # find overlap
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          mmy[3] <- Set.c
          mmy[4] <- Set.d
          mmy[5] <- Set.e
          mmy[6] <- Set.f
          mmy[7] <- Set.g
          if (length(which(mm6[1:no.group] == 0)) > 1) {
            mm6[which(mm3[Set.g, 1:no.group] != 0)] <- which(mm3[Set.g,1:no.group] != 0)

      for (Set.h in (Set.g+1):nrow(mm3)) {
        if (Set.h > nrow(mm3)) { break }
        if (all(mm3[Set.h, mm6!=0]==0) == FALSE) { ## without overlap
          next
        } else { # find overlap
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          mmy[3] <- Set.c
          mmy[4] <- Set.d
          mmy[5] <- Set.e
          mmy[6] <- Set.f
          mmy[7] <- Set.g
          mmy[8] <- Set.h
          if (length(which(mm6[1:no.group] == 0)) > 1) {
            mm6[which(mm3[Set.h, 1:no.group] != 0)] <- which(mm3[Set.h,1:no.group] != 0)

      for (Set.i in (Set.h+1):nrow(mm3)) {
        if (Set.i > nrow(mm3)) { break }
        if (all(mm3[Set.i, mm6!=0]==0) == FALSE) { ## without overlap
          next
        } else { # find overlap
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          mmy[3] <- Set.c
          mmy[4] <- Set.d
          mmy[5] <- Set.e
          mmy[6] <- Set.f
          mmy[7] <- Set.g
          mmy[8] <- Set.h
          mmy[9] <- Set.i
          if (length(which(mm6[1:no.group] == 0)) > 1) {
            mm6[which(mm3[Set.i, 1:no.group] != 0)] <- which(mm3[Set.i,1:no.group] != 0)

      for (Set.j in (Set.i+1):nrow(mm3)) {
        if (Set.j > nrow(mm3)) { break }
        if (all(mm3[Set.j, mm6!=0]==0) == FALSE) { ## without overlap
          next
        } else { # find overlap
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          mmy[3] <- Set.c
          mmy[4] <- Set.d
          mmy[5] <- Set.e
          mmy[6] <- Set.f
          mmy[7] <- Set.g
          mmy[8] <- Set.h
          mmy[9] <- Set.i
          mmy[10] <- Set.j
          if (length(which(mm6[1:no.group] == 0)) > 1) {
            mm6[which(mm3[Set.j, 1:no.group] != 0)] <- which(mm3[Set.j,1:no.group] != 0)

      for (Set.k in (Set.j+1):nrow(mm3)) {
        if (Set.k > nrow(mm3)) { break }
        if (all(mm3[Set.k, mm6!=0]==0) == FALSE) { ## without overlap
          next
        } else { # find overlap
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          mmy[3] <- Set.c
          mmy[4] <- Set.d
          mmy[5] <- Set.e
          mmy[6] <- Set.f
          mmy[7] <- Set.g
          mmy[8] <- Set.h
          mmy[9] <- Set.i
          mmy[10] <- Set.j
          mmy[11] <- Set.k
          if (length(which(mm6[1:no.group] == 0)) > 1) {
            mm6[which(mm3[Set.k, 1:no.group] != 0)] <- which(mm3[Set.k,1:no.group] != 0)

      for (Set.l in (Set.k+1):nrow(mm3)) {
        if (Set.l > nrow(mm3)) { break }
        if (all(mm3[Set.l, mm6!=0]==0) == FALSE) { ## without overlap
          next
        } else { # find overlap
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mmy[2] <- Set.b
          mmy[3] <- Set.c
          mmy[4] <- Set.d
          mmy[5] <- Set.e
          mmy[6] <- Set.f
          mmy[7] <- Set.g
          mmy[8] <- Set.h
          mmy[9] <- Set.i
          mmy[10] <- Set.j
          mmy[11] <- Set.k
          mmy[12] <- Set.l
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.l
    }  # if length Set.k
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.k
    }  # if length Set.j
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.j
    }  # if length Set.i
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.i
    }  # if length Set.h
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.h
    }  # if length Set.g
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.g
    }  # if length Set.f
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.f
    }  # if length Set.e
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.e
    }  # if length Set.d
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.d
    }  # if length Set.c
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.c
    } # if length Set.b
          mm8 <- rbind(mm8, mmy)
        } ## without overlap
      } # for Set.b

      if (Set.b == nrow(mm3)) {
        if (all(mm3[Set.b, mm6!=0]==0) == FALSE) {
          mmy <- matrix(0, 1, no.group+2)
          mmy[1] <- Set.a
          mm8 <- rbind(mm8, mmy)
        }
      }
    } # for Set.a

    for (i in 1: (nrow(mm8)-1)) {
      for (j in (i+1):nrow(mm8)) {
        if (all(mm8[j, 1:no.group] %in% mm8[i, 1:no.group]) == TRUE) { mm8[j, no.group+2] <- 1 }
      }  ## end for j
    }  ## end for i

  } ## end if no invariant groups


  txZ <- matrix(0, 1, no.group+2)

  mm8 <- subset(mm8, mm8[, no.group+2] != 1)
  mm8 <- mm8[2:nrow(mm8),]

  if (is.matrix(mm8) == FALSE) {  ## mm8 has only 1 row
    mm7 <- matrix(0, 1, no.group+2)
    for (i in 1: length(which(mm8[1:no.group] > 0))) {
      xx <- which(mm3[mm8[i], 1:no.group] != 0)
      for (j in 1:length(xx)) { mm7[xx[j]] <- paste0("F", factor.no, "I", nEP) } ## end for j
      if (length(xx != 0)) { nEP <- nEP + 1 }
      xx <- which(mm7[1:no.group] == 0)
      for (k in 1:length(xx)) { mm7[xx[k]] <- paste0("F", factor.no, "I", (nEP-1+k)) } ## end for k
      nEP <- nEP + k
    }  ## end for i
    txZ <- mm7
    txX <- mm7
    txX[which(txZ[1:no.group] != 0)] <- which(txZ[1:no.group] != 0)

  } else {
    for (i in 1: nrow(mm8)) {
      mm7 <- matrix(0, 1, no.group+2)
      for (j in 1: length(which(mm8[1, 1:no.group] > 0))) {
        xx <- which(mm3[mm8[i,j], 1:no.group] != 0)
        for (k in 1:length(xx)) { mm7[xx[k]] <- paste0("F", factor.no, "I", nEP) } ## end for k
        nEP <- nEP + 1
      } ## end for j
      xx <- which(mm7[1:no.group] == 0)
      for (j in 1:length(xx)) { mm7[xx[j]] <- paste0("F", factor.no, "I", (nEP-1+j)) } ## end for j
      nEP <- nEP + j
      txZ <- rbind(txZ, mm7)
    } ## end for i

    txZ <- txZ[2:nrow(txZ),]
    txX <- txZ
    for (i in 1:nrow(txZ)) { txX[i, which(txZ[i,1:no.group] != 0)] <- which(txZ[i,1:no.group] != 0) }
    for (i in 1: nrow(txZ)) { txZ[i, no.group+2] <- length(unique(txZ[i, 1:no.group])) }
    txX <- subset(txX, txZ[, no.group+2] == min(txZ[, no.group+2]))
    txZ <- subset(txZ, txZ[, no.group+2] == min(txZ[, no.group+2]))
  }


  txX[1:nrow(txX), (no.group+1)] <- Referent
  txZ[1:nrow(txZ), (no.group+1)] <- Referent

  txX[1:nrow(txX), (no.group+2)] <- Arg
  txZ[1:nrow(txZ), (no.group+2)] <- Arg

  txY <<- rbind(txY,txX)
  txYY <<- rbind(txYY,txZ)

 #$  print(txY)
#$  print(txYY)
  EP <<- nEP + 1

}  ## end function listanddelete.tx

## ==========  Finish Function List and Delete for TX ========== ##




## ==========  Create Function List and Delete ========== ##

listanddelete <- function(no_item = 5, no_nipair = 2, nipair = c(1,3,1,5)) {

  MI <- t(matrix(nipair,no_nipair,nrow=2))
  MI

  MIset <- no_item-2


  for (k in 1:MIset) {
    Noset <- no_item - k
    NIset <- factorial(no_item)/(factorial(no_item-k)*factorial(k))

    mm <- t(combn(no_item,Noset))
    for (j in 1:no_nipair) {
      a <- MI[j,1]
      b <- MI[j,2]
      for (i in 1:NIset) {
        count4 <- length(which(mm[i,] == a | mm[i,] == b))
        if (count4 > 1) {
          mm[i,] <- 0
        }
      }
    }
    for (ii in 1:NIset) {
      count4 <- sum(mm[ii,])
      if (count4 > 0) {
        print(mm[ii,])
      }
    }
  }

}  ## End (Function listanddelete)

## ==========  Finish Function List and Delete ========== ##





## ==========  Create Function correct.cfi ========== ##

corrected.cfi <- function(PSI.object, data, Groups, cluster="NULL") {

  ind.names <- lavNames(PSI.object) # PSI.object = PSI.Model.fit
  no.ind <- length(ind.names)
  no.group <- lavInspect(PSI.object, "ngroups")

  PSI.parTable <- parTable(PSI.object)
  psi.2 <- PSI.parTable[, c("lhs", "op", "group", "label")]
  psi.2 <- psi.2[psi.2$op == "~1", ]

  baseline_parTable <- lav_partable_independence(PSI.object)
  if (cluster == "NULL") {
    BASE2 <- lavaan::sem(baseline_parTable, data, group=Groups, # data = Example.A, group=Groups
                         meanstructure=TRUE,
                         marker.int.zero = TRUE,
                         estimator="MLR",
                         missing = "ML",
                         information="observed")
  } else {
    BASE2 <- lavaan::sem(baseline_parTable, data, group=Groups, # data = Example.A, group=Groups
                         meanstructure=TRUE,
                         marker.int.zero = TRUE,
                         estimator="MLR",
                         missing = "ML",
                         information="observed",
                         cluster = Cluster)
  }

  BASE.parTable <- parTable(BASE2)

  null.2 <- BASE.parTable[, c("lhs", "op", "group")]
  null.2 <- null.2[null.2$op == "~1", ]

  BASE3 <- merge(null.2, psi.2, by = c("lhs", "op", "group"))
  BASE3 <- BASE3[order(BASE3$lhs, BASE3$group), ]


  correct.baseline <- paste0("PSI.Null <- '", "\n")
  for (i in 1:no.ind) {
    correct.baseline <- paste0(correct.baseline, ind.names[i], " ~~ ", ind.names[i], "\n")
    FL.base <- paste0(ind.names[i], " ~ c(")
    for (j in 1:(no.group-1) ){
      if(BASE3[which(BASE3[,"lhs"] == ind.names[i] & BASE3[,"group"] == j), "label"] == "") {
        FL.base <- paste0(FL.base, "XFI", i,",")
      } else {
        FL.base <- paste0(FL.base, BASE3[which(BASE3[,"lhs"] == ind.names[i] & BASE3[,"group"] == j), "label"],",")
      } # end (if BASE3)
    } # end (for j)
    if(BASE3[which(BASE3[,"lhs"] == ind.names[i] & BASE3[,"group"] == no.group), "label"] == "") {
      FL.base <- paste0(FL.base, "XFI", i,")*1")
    } else {
      FL.base <- paste0(FL.base, BASE3[which(BASE3[,"lhs"] == ind.names[i] & BASE3[,"group"] == no.group), "label"], ")*1")
    } # end (if BASE3)
    #   print(FL.base)
    correct.baseline <- paste0(correct.baseline, FL.base, "\n")
  } # end (for i)
  correct.baseline <- paste0(correct.baseline, "'", "\n")
  eval(parse(text=correct.baseline))

  if (cluster == "NULL") {
    BASE.fit <- lavaan::sem(PSI.Null, data, group=Groups, # data = Example.A, Groups="Region"
                            meanstructure=TRUE,
                            marker.int.zero = TRUE,
                            estimator="MLR",
                            missing = "ML",
                            information="observed")
  } else {
    BASE.fit <- lavaan::sem(PSI.Null, data, group=Groups, # data = Example.A, Groups="Region"
                            meanstructure=TRUE,
                            marker.int.zero = TRUE,
                            estimator="MLR",
                            missing = "ML",
                            information="observed",
                            cluster = Cluster)
  } # end (if cluster)

  corrected.fit <- fitmeasures(PSI.object, baseline.model=BASE.fit)
  return(corrected.fit)
}  # end (function corrected.cfi)


# fitmeasures <- corrected.cfi(PSI.Model.fit, data=Example.A, Groups="Region")
# XX <- fitmeasures[c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr_bentler", "aic", "bic")]

## ==========  Finish Function correct.cfi ========== ##



## ==========  Create Function LGcorrect.cfi ========== ##

LGcorrected.cfi <- function(PSI.object, data, cluster="NULL") {

  PSI.parTable <- parTable(PSI.object)
  psi.2 <- PSI.parTable[, c("lhs", "op", "label")]
  psi.2 <- psi.2[psi.2$op == "~1", ]

  baseline_parTable <- lav_partable_independence(PSI.object)
  if (cluster == "NULL") {
    BASE2 <- lavaan::sem(baseline_parTable, data,
                         meanstructure=TRUE,
                         marker.int.zero = TRUE,
                         estimator="MLR",
                         missing = "ML",
                         information="observed")
  } else {
    BASE2 <- lavaan::sem(baseline_parTable, data,
                         meanstructure=TRUE,
                         marker.int.zero = TRUE,
                         estimator="MLR",
                         missing = "ML",
                         information="observed",
                         cluster = Cluster)
  }

  BASE.parTable <- parTable(BASE2)

  null.2 <- BASE.parTable[, c("lhs", "op")]
  null.2 <- null.2[null.2$op == "~1", ]

  BASE3 <- merge(null.2, psi.2, by = c("lhs", "op"))
  BASE3 <- BASE3[order(BASE3$lhs), ]


  delimiter <- "_T"
  BASE4 <- BASE3[which(BASE3[, "label"] == ""),]
  marker <- unique(sub(paste0("\\", delimiter, ".*"), "", BASE4[,"lhs"]))

  for (i in 1: length(marker)) {
    BASE3[which(sub(paste0("\\", delimiter, ".*"), "", BASE3[,"lhs"]) == marker[i]), "label"] <- paste0("XFI", i)
  }

  correct.baseline <- paste0("PSI.Null <- '", "\n")
  for (i in 1:nrow(BASE3)) {
    correct.baseline <- paste0(correct.baseline, BASE3[i, "lhs"], " ~~ ", BASE3[i, "lhs"], "\n")
    FL.base <- paste0(BASE3[i, "lhs"], " ~ ", BASE3[i, "label"], "*1")
    #    print(FL.base)
    correct.baseline <- paste0(correct.baseline, FL.base, "\n")
  } # end (for i)
  correct.baseline <- paste0(correct.baseline, "'", "\n")
  eval(parse(text=correct.baseline))

  if (cluster == "NULL") {
    BASE.fit <- lavaan::sem(PSI.Null, data,
                            meanstructure=TRUE,
                            marker.int.zero = TRUE,
                            estimator="MLR",
                            missing = "ML",
                            information="observed")
  } else {
    BASE.fit <- lavaan::sem(PSI.Null, data,
                            meanstructure=TRUE,
                            marker.int.zero = TRUE,
                            estimator="MLR",
                            missing = "ML",
                            information="observed",
                            cluster = Cluster)
  } # end (if cluster)

  corrected.fit <- fitmeasures(PSI.object, baseline.model=BASE.fit)
  return(corrected.fit)
}  # end (function LGcorrected.cfi)


# fitmeasures <- LGcorrected.cfi(PSI.Model.fit, data=Example.C)
# XX <- fitmeasures[c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr_bentler", "aic", "bic")]

## ==========  Finish Function LGcorrect.cfi ========== ##



## ==========  Create Function PAIRcorrect.cfi ========== ##

PAIRcorrected.cfi <- function(PSI.object, data, cluster="NULL") {

  PSI.parTable <- parTable(PSI.object)
  psi.2 <- PSI.parTable[, c("lhs", "op", "label")]
  psi.2 <- psi.2[psi.2$op == "~1", ]

  baseline_parTable <- lav_partable_independence(PSI.object)
  if (cluster == "NULL") {
    BASE2 <- lavaan::sem(baseline_parTable, data,
                         meanstructure=TRUE,
                         marker.int.zero = TRUE,
                         estimator="MLR",
                         missing = "ML",
                         information="observed")
  } else {
    BASE2 <- lavaan::sem(baseline_parTable, data,
                         meanstructure=TRUE,
                         marker.int.zero = TRUE,
                         estimator="MLR",
                         missing = "ML",
                         information="observed",
                         cluster = Cluster)
  }

  BASE.parTable <- parTable(BASE2)

  null.2 <- BASE.parTable[, c("lhs", "op")]
  null.2 <- null.2[null.2$op == "~1", ]

  BASE3 <- merge(null.2, psi.2, by = c("lhs", "op"))
  BASE3 <- BASE3[order(BASE3$lhs), ]


  delimiter <- "_G"
  BASE4 <- BASE3[which(BASE3[, "label"] == ""),]
  marker <- unique(sub(paste0("\\", delimiter, ".*"), "", BASE4[,"lhs"]))

  for (i in 1: length(marker)) {
    BASE3[which(sub(paste0("\\", delimiter, ".*"), "", BASE3[,"lhs"]) == marker[i]), "label"] <- paste0("XFI", i)
  }

  correct.baseline <- paste0("PSI.Null <- '", "\n")
  for (i in 1:nrow(BASE3)) {
    correct.baseline <- paste0(correct.baseline, BASE3[i, "lhs"], " ~~ ", BASE3[i, "lhs"], "\n")
    FL.base <- paste0(BASE3[i, "lhs"], " ~ ", BASE3[i, "label"], "*1")
    #    print(FL.base)
    correct.baseline <- paste0(correct.baseline, FL.base, "\n")
  } # end (for i)
  correct.baseline <- paste0(correct.baseline, "'", "\n")
  eval(parse(text=correct.baseline))

  if (cluster == "NULL") {
    BASE.fit <- lavaan::sem(PSI.Null, data,
                            meanstructure=TRUE,
                            marker.int.zero = TRUE,
                            estimator="MLR",
                            missing = "ML",
                            information="observed")
  } else {
    BASE.fit <- lavaan::sem(PSI.Null, data,
                            meanstructure=TRUE,
                            marker.int.zero = TRUE,
                            estimator="MLR",
                            missing = "ML",
                            information="observed",
                            cluster = Cluster)
  } # end (if cluster)

  corrected.fit <- fitmeasures(PSI.object, baseline.model=BASE.fit)
  return(corrected.fit)
}  # end (function PAIRcorrected.cfi)


# fitmeasures <- LGcorrected.cfi(PSI.Model.fit, data=Example.C)
# XX <- fitmeasures[c("chisq.scaled", "df.scaled", "pvalue.scaled", "rmsea.scaled", "cfi.scaled", "tli.scaled", "srmr_bentler", "aic", "bic")]

## ==========  Finish Function PAIRcorrect.cfi ========== ##


