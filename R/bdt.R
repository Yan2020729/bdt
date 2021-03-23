

############################################################################################################
#------------------------------------------------ Stage 2 --------------------------------------------------
############################################################################################################

#---------------------- produce Q1, Q0, Ystar in observed data -------------------------
# input: Y, A, W, seed, Qform, gform, outcome type(continuous or binary)
# output: an augmented data with new outcome Ystar and Q1, Q0
# (for continuous outcome, Y*, Q1 and Q0 are all bounded, so also output the range of the three quantities)

.generateQnewY <- function(Y, A, W, seed, Qform, gform, outcome_type){

  Qform <- stringr::str_replace_all(string = Qform, pattern = "A+", repl = "1")
  ObsData <- data.frame(Y, A, W)
  ab <- NULL

  glm_g <- glm(gform, data = ObsData, family = "binomial")

  # split the observed dataset to two treatment groups
  base1 <- ObsData[ObsData$A ==1,]
  base0 <- ObsData[ObsData$A ==0,]
  designM <- model.matrix (as.formula(Qform), ObsData)

  # regress the outcome separately in two groups
  set.seed(seed)
  if (identical(outcome_type,"continuous")){
    glm_Q1 <- glm(Qform, data = base1)
    glm_Q0 <- glm(Qform, data = base0)

    Q1pred <- designM[, !is.na(coef(glm_Q1))] %*% glm_Q1$coeff[!is.na(coef(glm_Q1))]
    Q0pred <- designM[, !is.na(coef(glm_Q0))] %*% glm_Q0$coeff[!is.na(coef(glm_Q0))]
    Ypred1 <- rnorm(dim(ObsData)[1], mean = Q1pred, sd= sigma(glm_Q1))
    Ypred0 <- rnorm(dim(ObsData)[1], mean = Q0pred, sd= sigma(glm_Q0))
    Ystar <- ifelse(ObsData$A==1, Ypred1, Ypred0)

    ab <- range(Q1pred, Q0pred, Ystar)
    Q1pred = (Q1pred-ab[1])/diff(ab); Q0pred = (Q0pred-ab[1])/diff(ab); Ystar = (Ystar-ab[1])/diff(ab)
  }

  if (identical(outcome_type,"binary")){
    glm_Q1 <- glm(Qform, data = base1, family = "binomial")
    glm_Q0 <- glm(Qform, data = base0, family = "binomial")

    Q1pred = plogis( designM[, !is.na(coef(glm_Q1))] %*% coef(glm_Q1)[!is.na(coef(glm_Q1))] )
    Q0pred = plogis( designM[, !is.na(coef(glm_Q0))] %*% coef(glm_Q0)[!is.na(coef(glm_Q0))] )

    Ypred1 <- rbinom(dim(ObsData)[1], 1, Q1pred)
    Ypred0 <- rbinom(dim(ObsData)[1], 1, Q0pred)
    Ystar <- ifelse(ObsData$A==1, Ypred1, Ypred0)
  }

  augObsData <- data.frame(ObsData[,!names(ObsData)%in% c("Y")], Q1pred, Q0pred, Ystar)
  return(list(augObsData = augObsData, glm_g = glm_g, ab = ab))
}




############################################################################################################
#------------------------------------------------ Stage 3 --------------------------------------------------
############################################################################################################

.bound <- function(x, bounds){
  x[x<min(bounds)] <- min(bounds); x[x>max(bounds)] <- max(bounds)
  return(x)
}

# use the bound function to compute the weight
.weight <- function (x, gbound){
  if (gbound ==0){weigh <- 1/x}
  if (gbound!=0){weigh <- 1/.bound(x, c(gbound, 1-gbound))}
  return(weigh)
}

#----------------------------- produce an augmented bootstrap data -------------------------------

# This function is using the coeffs of the treatment regression and SL to predict weights
#        in the bootstrap data then create au augmented bootstrap dataset and the true value as well
# input: obsData, seed, gGLM, the coeffs of outcome regression, SL.library, gbound and range ab
# output: one augmented bootstrap dataset, the "true" ate and p(A=1|W)

.sim_data <- function(augObsData, seed, gGLM, gform, glm_g, SL.library, gbound, ab)
{

  # sample n subject with replacement
  set.seed(seed)
  n <- dim(augObsData)[1]
  resamp <- sample(1:n, n, replace = TRUE)
  bootData <- augObsData[resamp,]

  # obtain the true value
  if (!is.null(ab)){tru <- mean(bootData$Q1pred-bootData$Q0pred)*diff(ab)
  } else {tru <- mean(bootData$Q1pred-bootData$Q0pred)}

  # estimate ps by GLM or SL
  if (gGLM){
    designM <- model.matrix (as.formula(gform), bootData)
    ps <- plogis( designM[, !is.na(coef(glm_g))] %*% coef(glm_g)[!is.na(coef(glm_g))] )
  } else {
    # library(SuperLearner)
    ps <- suppressWarnings(SuperLearner(Y = augObsData$A, X= augObsData[, !names(augObsData) %in% c("A", "Q0pred", "Q1pred", "Ystar")],
                                           newX = bootData[, !names(bootData) %in% c("A", "Q0pred", "Q1pred", "Ystar")],
                                           family = binomial(),SL.library = SL.library)$SL.predict)
  }
  weight <- ifelse(bootData$A == 0, 1- ps, ps)
  weigBound <- .weight(weight, gbound)
  data.sim <- data.frame(bootData, w = weigBound)

  return(list(data.sim = data.sim, true = tru, ps = ps))
}


############################################################################################################
#------------------------------------------------ Stage 5 --------------------------------------------------
############################################################################################################

# this function is to covert categorical variables to corresponding dichotomous variables.
# depend on the data, there are two choices: remove first dummy or remove most frequent dummy
.convertCategor <- function(W, gform, Qform, remove_first_dummy, remove_most_frequent_dummy){
  if (any(lapply(W, is.factor) == TRUE)){

    cat_names <- names(W)[lapply(W, is.factor) == TRUE]
    W_dums = fastDummies::dummy_cols(W, remove_first_dummy = remove_first_dummy,
                                     remove_most_frequent_dummy=remove_most_frequent_dummy)
    W_dum = W_dums[, !names(W_dums) %in% c(cat_names)]

    noSpace.gform <- stringr::str_replace_all(string=gform, pattern=" ", repl="")
    g_val <- chartr(old = "+", new = " ", substring(noSpace.gform, 3))
    g_val <- unlist(strsplit( g_val, split = " "))

    gform_cat <- stringr::str_c(c(unlist(sapply(1:length(g_val), function(x)
      names(W_dum)[stringr::str_detect(names(W_dum), g_val[x])]))), collapse = "+")
    gform <- paste0("A~", gform_cat)

    noSpace.Qform <- stringr::str_replace_all(string=Qform, pattern=" ", repl="")
    Q_val <- chartr(old = "+", new = " ", substring(noSpace.Qform, 5))
    Q_val <- unlist(strsplit( Q_val, split = " "))

    Qform_cat <- stringr::str_c(c(unlist(sapply(1:length(Q_val), function(x)
        names(W_dum)[stringr::str_detect(names(W_dum), Q_val[x])]))), collapse = "+")
    Qform <- paste0("Y~A+", Qform_cat)

    return(list(W=W_dum, gform=gform, Qform=Qform))
  } else {return(list(W=W, gform=gform, Qform=Qform))}
}


# output the message of method used based on the input of gform and SL.library
.check_gMeth <- function(gGLM, gform, SL.library){
  if (gGLM){
    if (!is.null(gform)){
      if (!is.null(SL.library)) {# if user set gGLM is TRUE, define both gform and SL function simultaneously
        cat("Estimating ate using user-supplied regression formula and ignoring supplied SuperLearner functions for g \n\n")}
      else (cat("Estimating ate using user-supplied regression formula for g \n\n"))
    }
    else {
      if (!is.null(SL.library)) {# if user set gGLM is TRUE but define SL function instead of gform
        cat("Estimating ate ignoring supplied SuperLearner functions and using default regression formula 'A~ W' for g \n\n")}
      else (cat("Estimating ate using default regression formula 'A~ W' for g \n\n"))
    }
  }
  else {
    if (!is.null(SL.library)){
      if (!is.null(gform)){
        cat("Estimating ate using user-supplied SuperLearner functions and ignoring supplied regression formula for g \n\n")}
      else {cat("Estimating ate using user-supplied SuperLearner functions for g \n\n")}
    }
    else (stop("Please specify SuperLearner library for g. \n\n"))
  }
}


#------------------------ main function of TMLE to estimate ATE with flexible g (SL or GLM) ----------------------

# input: one augmented bootstrap dataset, and range ab if outcome is continuous
# returns: the estimated ate, sd and 95% confidence intervals

TMLE_ate <- function (bootdata, ab){

  n = dim(bootdata)[1]; Ystar <- bootdata$Ystar; A <- bootdata$A
  Q1m = bootdata$Q1pred; Q0m = bootdata$Q0pred; weight = bootdata$w

  update1 <- glm(Ystar~1, subset=(A==1), offset = qlogis(Q1m), weights=weight, family=quasibinomial(), data = bootdata)
  update0 <- glm(Ystar~1, subset=(A==0), offset = qlogis(Q0m), weights=weight, family=quasibinomial(), data = bootdata)
  Q1m_updated <- plogis(qlogis(Q1m) + as.numeric(coef(update1)))
  Q0m_updated <- plogis(qlogis(Q0m) + as.numeric(coef(update0)))

  if (!is.null(ab)){
    ate_tmle <- mean(Q1m_updated- Q0m_updated, na.rm = TRUE)*diff(ab)
    IF <- ((Ystar-Q1m_updated)*A*weight + Q1m_updated)-((Ystar - Q0m_updated)*(1-A)*weight + Q0m_updated)-mean(Q1m_updated- Q0m_updated, na.rm = TRUE) # !!!! last term is not ate_tmle
    sd <- sqrt( sum(IF^2)/(n^2))*diff(ab)
  } else {
    ate_tmle <- mean(Q1m_updated- Q0m_updated, na.rm = TRUE)
    IF <- ((Ystar-Q1m_updated)*A*weight + Q1m_updated)-((Ystar - Q0m_updated)*(1-A)*weight + Q0m_updated)-mean(Q1m_updated- Q0m_updated, na.rm = TRUE) # !!!! last term is not ate_tmle
    sd <- sqrt( sum(IF^2)/(n^2))
  }
  ci <- ate_tmle + c(-1,1)*1.96*(sd)

  return(list(ATE_tmle = ate_tmle, SE_tmle =sd, CI_tmle=ci))
}



#-------------------------- main function of AIPTW to estimate ATE with flexible g (SL or GLM) -----------------------

# input: one augmented bootstrap dataset, and range ab if outcome is continuous
# returns: the estimated ate, se and 95% confidence intervals

AIPTW_ate <- function (bootdata, ab){

  n = dim(bootdata)[1]; Ystar <- bootdata$Ystar; A <- bootdata$A
  Q1m = bootdata$Q1pred; Q0m = bootdata$Q0pred; weight = bootdata$w

  if (!is.null(ab)){
    mu1 <- sum((Ystar - Q1m)*A*weight)/n + mean(Q1m); mu0 <- sum((Ystar - Q0m)*(1-A)*weight)/n + mean(Q0m)
    ate_aiptw <- (mu1-mu0)*diff(ab)
    IF <- ((Ystar-Q1m)*A*weight + Q1m)-((Ystar - Q0m)*(1-A)*weight + Q0m)-(mu1-mu0)
    sd <- sqrt( sum(IF^2)/(n^2))*diff(ab)
  } else {
    mu1 <- sum((Ystar - Q1m)*A*weight)/n + mean(Q1m); mu0 <- sum((Ystar - Q0m)*(1-A)*weight)/n + mean(Q0m)
    ate_aiptw <- (mu1-mu0)
    IF <- ((Ystar-Q1m)*A*weight + Q1m)-((Ystar - Q0m)*(1-A)*weight + Q0m)-(mu1-mu0)
    sd <- sqrt( sum(IF^2)/(n^2))
  }
  ci <- ate_aiptw + c(-1,1)*1.96*(sd)

  return(list(ATE_aiptw = ate_aiptw, SE_aiptw = sd, CI_aiptw = ci))
}


############################################################################################################
#------------------------------------------------ Stage 6 --------------------------------------------------
############################################################################################################

#-------------- helper function to check all input variables --------------------
.check_var2 <- function(Y, A, W, outcome_type, gGLM, M, gbound, gform, SL.library, seed){
  # check the dimension of W and Y, and A
  if(nrow(W)!= length(Y) | length(A)!= length(Y) ){stop("The data dimensions do not correspond.")}

  # check missing values
  miss <- sapply(data.frame(Y, A, W), function(x)sum(is.na(x)))
  if (sum(miss) > 0) {stop("There are missing values in the dataset. ")}

  # check gbound
  if (gbound > 1 || gbound < 0){stop("gbound must be between 0 and 1.")}

  # check M
  if (class(M)!= "numeric" || M < 1){stop("M must be a number more than 0.")}
  else if (M <= 50) {cat("Small number of replications may lead to unreliable estimates. \n\n")}

  # check binary treatment A
  if (is.factor(A)){if (length(levels(A)) != 2L) stop("A must be a binary variable")}
  if (is.numeric(A)) {
    uniq_a <- unique(A)
    if ( length(unique(uniq_a)) > 2L || any(as.integer(uniq_a) != uniq_a) ) {stop("A must be a binary variable")}
    if (any(! uniq_a %in% c(0,1))) stop("A must be defined as 1 for treatment and 0 for control.")
  }
  if (!is.numeric(A) & !is.factor(A) & !is.logical(A)) {stop("A must be a binary variable. ")}

  # check outcome Y
  y_type <- NULL
  if (is.factor(Y)){
    if (length(levels(Y)) != 2L) stop("Only continuous or binary outcomes are accepted.")
    else (y_type <-"binary")
  }
  else if (is.logical(Y)) {y_type <-"binary"}
  else if (is.numeric(Y)) {
    uniq_y <- unique(Y)
    if ( length(unique(uniq_y)) > 2L || any(as.integer(uniq_y) != uniq_y) ) {y_type <-"continuous"}
    else (y_type <-"binary")
  }
  else (stop("Only continuous or binary outcomes are accepted."))

  # check outcome type
  if (class(outcome_type)!="character"|!outcome_type %in% c("continuous", "binary"))
  {stop("outcome_type must be 'continuous' or 'binary'.")}

  # check the correspondence of outcome_type and Y value
  if ((outcome_type =="continuous"& y_type!="continuous") || (outcome_type =="binary"& y_type!="binary")) {
    stop("The outcome_type does not correspond to type of Y value.")}

  # check gGLM and verbose
  if (is.null(gGLM) || class(gGLM)!="logical"){stop("You should define 'gGLM' to be TRUE or FALSE.")}

  # remove the space in the string of input gform
  if (!is.null(gform)){
    noSpace.gform <- stringr::str_replace_all(string = gform, pattern = " ", repl= "")
    if (substr(noSpace.gform, 1, 2)!= "A~"){stop("Error inputs for gform occur. ")}}

  # check seed
  if (!is.null(seed)){if(length(seed)!=M+1){stop("Number of Seed must be M+1.")}}
}



#----------------- main function to compute the bias of the estimates on AIPTW,TMLE ------------------

# input: outcome- Y(binary and/or continuous); binary treatment A; a set of potential confounders- W
# outcome_type -- "continuous" or "binary"
# M is the number of replications
# seed is needed if user want to compare estimators by different estimate methods of ps; length should equal to M+1
# gform -- optional regression formula for estimating P(A=1|W), default is A~W
# Qform -- regression formula for estimating P(Y=1|A, W), default is Y~A,W
# gbound -- bound on P(A=1|W), defaults to 0.025
# gGLM -- logical: if TRUE, use GLM; otherwise use SL
# SL.library -- prediction algorithms for data adaptive estimation of g
# remove_first_dummy -- for categorical covariates, if true remove the first dummy of each covariates
# remove_most_frequent_dummy -- for categorical covariates, if true remove the most frequently observed category
# returns: the estimated ate and 95% confidence intervals for AIPTW, TMLE and pooled ps by glm or SL.

bdt <- function (Y, A, W, outcome_type, M, seed = NULL, gGLM = NULL, gform = NULL, gbound = 0.025, SL.library = NULL,
                 remove_first_dummy = FALSE, remove_most_frequent_dummy = FALSE){

  bias_TMLE <- NULL; bias_AIPTW <- NULL; sd_TMLE <- NULL; sd_AIPTW <- NULL; ps_M <- NULL; trueM <- NULL
  cov_TMLE <- cov_AIPTW <- 0

  # check all inputs
  .check_var2(Y, A, W, outcome_type=outcome_type, gGLM=gGLM, M=M, gbound=gbound, gform=gform, SL.library=SL.library, seed=seed)
  .check_gMeth(gGLM=gGLM, gform=gform, SL.library=SL.library)

  # create the outcome and default treatment regression formula
  Qform <- paste ("Y~A+", paste(colnames(W), collapse = "+"))
  if(is.null(gform)){gform=paste("A~", paste(colnames(W), collapse = "+"))}

  # test factors in W and change to dummy variables
  coverCate <- .convertCategor(W=W, gform=gform, Qform=Qform, remove_first_dummy=remove_first_dummy, remove_most_frequent_dummy=remove_most_frequent_dummy)
  W_dumm = coverCate$W; Qform_dumm = coverCate$Qform; gform_dumm = coverCate$gform

  # generate augmented observed data containing Q1, Q0 and Ystar
  augObs <- .generateQnewY(Y, A, W=W_dumm, seed = seed[1], Qform = Qform_dumm, gform = gform_dumm, outcome_type = outcome_type)
  augObsData = augObs$augObsData; glm_g = augObs$glm_g; ab = augObs$ab

  # loops to create M bootstrap data sets
  for (i in 1:M){
    setseed = seed[i+1]
    augBoot <- .sim_data(augObsData=augObsData, seed= setseed, gGLM = gGLM, gform=gform_dumm, glm_g=glm_g,
                         SL.library = SL.library, gbound=gbound, ab = ab)
    true_eff <- augBoot$true
    trueM[i] <- true_eff

    tmle_res <- TMLE_ate(bootdata = augBoot$data.sim, ab = ab)
    aiptw_res <- AIPTW_ate(bootdata = augBoot$data.sim, ab = ab)

    # compute the bias, se and ci of tmle
    bias_TMLE[i] <- tmle_res$ATE_tmle-true_eff
    sd_TMLE[i] <- tmle_res$SE_tmle
    cov_TMLE <- cov_TMLE + ifelse(true_eff >= tmle_res$CI_tmle[1] & true_eff <= tmle_res$CI_tmle[2], 1, 0 )

    # compute the bias, se and ci of aiptw
    bias_AIPTW[i] <- aiptw_res$ATE_aiptw-true_eff
    sd_AIPTW[i] <- aiptw_res$SE_aiptw
    cov_AIPTW <- cov_AIPTW + ifelse(true_eff >= aiptw_res$CI_aiptw[1] & true_eff <= aiptw_res$CI_aiptw[2], 1, 0 )

    #---------- pooled ps ---------
    ps <- data.frame(A= augBoot$data.sim$A, ps = augBoot$ps)
    ps_M <- rbind(ps_M, ps)
  }

  results <- list(true_Effect = trueM, gbound =  gbound, ps_M = ps_M, seeds = seed, gGLM= gGLM,
                  bias_AIPTW = bias_AIPTW, bias_TMLE = bias_TMLE,
                  sd_AIPTW = sd_AIPTW, sd_TMLE = sd_TMLE,
                  cov_AIPTW = cov_AIPTW/M, cov_TMLE = cov_TMLE/M)
  class(results) <- "bdt"
  return(results)
}



#----------------- print functions to output the estimated ate bias of AIPTW, TMLE -------------------
# input the results of main bias_boot function
# return the bias of two estimated ate
bias.bdt <- function(object){
  if(identical(class(object), "bdt")){
    cat("Eestimated bias of ATE with AIPTW and TMLE:\n\n")
    bias.ate <- data.frame(cbind(object$bias_AIPTW, object$bias_TMLE))
    colnames(bias.ate) <- c("bias(AIPTW)", "bias(TMLE)" )
  } else {stop("Object must have class 'bdt'")}

  return(bias.ate)
}
bias <- function(object){ UseMethod("bias",object)}




#----------------- summary of BDT -------------------
# create BDT summary object

summary.bdt <- function(object,...){
  if(identical(class(object), "bdt")){

    smda = data.frame(object$ps_M)

    ps1.sum = summary(smda$ps)
    ps0.sum = summary(1-smda$ps)
    ps11.sum = summary(smda[smda$A==1,]$ps)
    ps00.sum = summary(1-smda[smda$A==0,]$ps)

    ps <- data.frame(rbind(ps1.sum, ps0.sum, ps11.sum, ps00.sum))
    dimnames(ps) <- list(c("P(A=1|X) for all", "P(A=0|X) for all",
                           "P(A=1|X) for A=1", "P(A=0|X) for A=0"),
                         c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."))

    bias <- data.frame(rbind(summary(object$bias_AIPTW), summary(object$bias_TMLE)))
    dimnames(bias) <- list(c("Bias of AIPTW", "Bias of TMLE"),
                           c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."))

    sd <- paste0("SE of AIPTW: ",round(mean(object$sd_AIPTW),4),
                 ";  SE of TMLE: ",round(mean(object$sd_TMLE),4))

    coverage <- paste0("Coverage rate of AIPTW: ",round(object$cov_AIPTW,4),
                       ";  Coverage rate of TMLE: ",round(object$cov_TMLE,4))

    summary.bdt <- list(true.effect = mean(object$true_Effect), gbound = object$gbound, gGLM = object$gGLM,
                      bias.aiptw = summary(object$bias_AIPTW),
                      bias.tmle = summary(object$bias_TMLE),
                      sd.aiptw = mean(object$sd_AIPTW),
                      sd.tmle = mean(object$sd_TMLE),
                      cov.aiptw = round(object$cov_AIPTW, 4),
                      cov.tmle = round(object$cov_TMLE, 4),
                      bias = bias,  coverage = coverage, sd = sd, ps = ps)

    class(summary.bdt) <- "summary.bdt"
  } else {
    stop("Object must have class 'bdt'")
    summary.bdt <- NULL}
  return(summary.bdt)
}



#-------------print.summary.BDT------------------
# print BDT summary object

print.summary.bdt <- function(x,...){
  if(identical(class(x), "summary.bdt")){
    specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall = k))

    cat("\n\nTrue simulated data effect: ", specify_decimal(x$true.effect, 5), " \n")
    if (x$gGLM){
      cat("\nMean and median bias, SE and coverage rates of AIPTW and TMLE with the ps estimation by GLM:\n\n")
    } else {
      cat("\nMean and median bias, SE and coverage rates of AIPTW and TMLE with the ps estimation by SL:\n\n")
    }
    a <- specify_decimal(x$bias.aiptw[4], 5); t <- specify_decimal(x$bias.tmle[4], 5)
    u <- specify_decimal(x$bias.aiptw[3], 5); v <- specify_decimal(x$bias.tmle[3], 5)
    m <- specify_decimal(x$sd.aiptw, 5); n <- specify_decimal(x$sd.tmle, 5)

    res_data <- data.frame(c("AIPTW"," TMLE"), c(a, t), c(u, v), c(m, n), c(x$cov.aiptw, x$cov.tmle))
    names(res_data) <- c("     ","MeanBias",  "Med.Bias", "   SE ", " Cov.")
    gdata::write.fwf(res_data)

    cat("\n\nSummary of propensity scores truncated at gbound", x$gbound, " \n\n")
    print(x$ps)
  }
}



#--------------------------- plots of the results -------------------------------
# input the results of main bias_boot function
# return two plots:
#   1. the boxplots and points of bias of two estimated ate and coverage rates
#   2. the density plots of pooled propensity scores
# For the boxplot: y-axis represents the bias of estimates; x-axis represents estimates: AIPTW, TMLE
#                 coverage rates are shown in red boxes
# For the density plot: return 2 x 2 density plots (A, B, C, D) of the log of the estimated weights for
#         A--treatment A=1 in subset of subjects with A=1
#         B--treatment A=0 in subset of subjects with A=0
#         C--treatment A=1 for all subjects
#         D--treatment A=0 for all subjects

plot.bdt <- function(x, xlab = NULL, ylab = "Bias", outlierSize= 0.4, notch = FALSE,
                     pointsize = 0.7, pointcol = "blue", labelsize = 0.35, labelcol = "red",
                     linetype = "dotted", linecol = "darkred", ...){
  if(identical(class(x), "bdt")){

    #--------------------- boxplot of bias ate ------------------
    N <- length(x$bias_TMLE)
    type_est <- c(rep("AIPTW",N),rep("TMLE",N) )
    Bias_ATE <- c(x$bias_AIPTW, x$bias_TMLE)
    est_data <- data.frame(type_est, Bias_ATE)
    colnames(est_data) <- c("type_est", "Bias_ATE")
    mt = ifelse(x$gGLM, "GLM for g", "SL for g")
    title <- paste0("Boxplots of the ATE bias with ", x$gbound, " gbound under ", mt)
    type = c("AIPTW", "TMLE")
    high = min(est_data$Bias_ATE)-0.03

    coverage = round(c(x$cov_AIPTW, x$cov_TMLE), 3)
    anno.data <- data.frame(type, high, coverage)
    # library(ggplot2)
    est_plot <- ggplot(est_data, aes(x = as.factor(type_est), y = Bias_ATE))+
      geom_boxplot(outlier.size = outlierSize, fill= "cornflowerblue", notch = notch)+
      geom_rug(color = "black") + ggtitle(title)+
      geom_point(position = "jitter", size=pointsize, color=pointcol, alpha=0.5)+
      geom_label(data = anno.data, aes(x = type,  y = high, label = paste("Cov: ", coverage,sep="")),
                 label.padding = unit(0.40, "lines"), # Rectangle size around label
                 label.size = labelsize, colour = labelcol, size = 3.5, fill = "white") +
      geom_hline(yintercept = 0,col = linecol, linetype = linetype) +
      ylim(min(Bias_ATE)-0.03, max(Bias_ATE) + 0.03) +
      ylab(ylab) + xlab(xlab) +
      theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=0))


    #------------------- density plot of ps ------------------
    object <- x
    n <- length(object$ps_M$ps)
    smda <- data.frame(object$ps_M)
    if (object$gGLM){ gmethod <- "GLM"} else {gmethod <- "SL"}

    ps1 <- data.frame(as.factor(c(rep(gmethod, n))),c(log(1/smda$ps)))
    ps0 <- data.frame(as.factor(c(rep(gmethod, n))),c(log(1/(1-smda$ps))))
    ps.A1 <- data.frame(as.factor(c(rep(gmethod, length(smda[smda$A ==1,]$ps)))),c(log(1/smda[smda$A ==1,]$ps)))
    ps.A0 <- data.frame(as.factor(c(rep(gmethod, length(smda[smda$A ==0,]$ps)))),c(log(1/(1-smda[smda$A ==0,]$ps))))
    colnames(ps1) <- colnames(ps0) <- colnames(ps.A1) <- colnames(ps.A0) <- c(gmethod, "Propensity")

    ps.datas <- list(ps.A1, ps.A0, ps1, ps0)
    xlab_name <- paste0(c("log(1/g_1)[A=1]", "log(1/g_1)[A=0]", "log(1/g_1)", "log(1/g_0)"), " on ", gmethod)

    ps_plot <- sapply(1:4, function(y) list(ggplot(ps.datas[[y]],
                                          aes(fill = eval(paste0("ps.datas[[y]]$", gmethod)), x = ps.datas[[y]]$Propensity))+
                                          geom_density(alpha = 0.3) + theme(axis.text.x = element_text(angle = 0))+
                                          xlab(xlab_name[y]) + theme(legend.position = "none")))

    ps_plots <- ggpubr::ggarrange(ps_plot[[1]], ps_plot[[2]], ps_plot[[3]], ps_plot[[4]], labels = c("A","B","C","D"), ncol = 2, nrow = 2)
    fig_box_density <- list(boxplot_bias_ATE = est_plot, densityplot_ps = ps_plots)
  }
  else {
    stop("Object must have class 'bdt'")
    fig_box_density <- NULL
  }
  par(ask=TRUE)
  return(fig_box_density)
}

