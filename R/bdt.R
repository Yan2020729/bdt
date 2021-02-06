

############################################################################################################
#------------------------------------------------ Stage 2 --------------------------------------------------
############################################################################################################


#---------------------- function to fit glm on Y in the observed data -------------------------
# input: Y, A, W, Qform, outcome type(continuous or binary)
# output: the coefficients of the glm and valide variable names for two treatment groups

.fit_glm_Y <- function(Y, A, W, Qform, outcome_type){

  Qform <- stringr::str_replace_all(string = Qform, pattern = "A+", repl = "1")
  ObsData <- data.frame(Y, A, W)

  # split the observed dataset to two treatment groups
  base1 <- ObsData[ObsData$A ==1,]
  base0 <- ObsData[ObsData$A ==0,]

  # regress the outcome separately in two groups
  if (identical(outcome_type,"continuous")){
    regY1 <- lm(Qform, data = base1)
    regY0 <- lm(Qform, data = base0)
  }
  if (identical(outcome_type,"binary")){
    regY0 <- glm(Qform, data = base0, family = "binomial")
    regY1 <- glm(Qform, data = base1, family = "binomial")
  }
  return(list(fit_y0=regY0, fit_y1=regY1))
}



############################################################################################################
#------------------------------------------------ Stage 3 --------------------------------------------------
############################################################################################################


#----------------------------- function to produce one bootstrap data -------------------------------

# This function is using the coeffs of the outcome regression in the original dataset to predict
#       the outcome in the bootstrap data then create the bootstrap dataset
# input: A, W, outcome type(continuous or binary), and the coeffs of outcome regression and valide
#        variable names for two treatment groups
# output: one bootstrap dataset

.sim_data <- function(A, W, Qform, outcome_type, fit_y1, fit_y0)
{

  ObsData <- data.frame(A, W)
  base1 <- ObsData[ObsData$A == 1,]
  base0 <- ObsData[ObsData$A == 0,]

  # sample n subject with replacement
  n <- dim(ObsData)[1]
  resamp <- sample(1:n, n, replace = TRUE)
  bootData <- ObsData[resamp,]

  # split the sample data into treated and untreated
  BaseBData1 <- bootData[bootData$A == 1,]
  BaseBData0 <- bootData[bootData$A == 0,]

  # create design matrix on simulated data
  Qform <- substring(stringr::str_replace_all(string = Qform, pattern = "A+", repl = "1"),2)
  designM1 <- model.matrix (as.formula(Qform), BaseBData1)
  designM0 <- model.matrix (as.formula(Qform), BaseBData0)

  # get the new outcome based on the coeffs of glm on original data
  if (identical(outcome_type,"continuous")){
    mu1 = designM1[, !is.na(coef(fit_y1))] %*% fit_y1$coeff[!is.na(coef(fit_y1))]
    mu0 = designM0[, !is.na(coef(fit_y0))] %*% fit_y0$coeff[!is.na(coef(fit_y0))]
    Ypred1 <- rnorm(n=dim(BaseBData1)[1], mean= mu1, sd= summary(fit_y1)$sigma)
    Ypred0 <- rnorm(n=dim(BaseBData0)[1], mean= mu0, sd= summary(fit_y0)$sigma)
  }
  if (identical(outcome_type,"binary")){
    Ypred1 <- rbinom(dim(BaseBData1)[1],1, plogis( designM1[, !is.na(coef(fit_y1))] %*% coef(fit_y1)[!is.na(coef(fit_y1))] ))
    Ypred0 <- rbinom(dim(BaseBData0)[1],1, plogis( designM0[, !is.na(coef(fit_y0))] %*% coef(fit_y0)[!is.na(coef(fit_y0))] ))
  }

  # combine the two group data to obtain a dataframe
  data1 <- data.frame(Ypred1,BaseBData1)
  data0 <- data.frame(Ypred0,BaseBData0)
  colnames(data1) <- colnames(data0) <- c("Y","A",colnames(W))

  data.sim <- data.frame(rbind(data1,data0))
  return(data.sim)
}




############################################################################################################
#------------------------------------------------ Stage 4 --------------------------------------------------
############################################################################################################

#-------------------------- main function to compute the true effect of bootstrap data ---------------------

# This function is to compute the true effect of one bootstrap dataset.
# input: one bootstrap dataset, outcome type (binary or continuous) and the coeffs of outcome regression
#        and valid variable names for two treatment groups
# output: one numerical true value

.trueVal <- function (Y, A, W, Qform, outcome_type, fit_y1, fit_y0){

  # create design matrix on Qform
  Qform <- stringr::str_replace_all(string = Qform, pattern = "A+", repl = "1")
  designM <- model.matrix (as.formula(Qform), data.frame(Y, A, W))

  # compute two potential outcomes of original data
  if (identical(outcome_type,"continuous")){
    yy1 <- designM[, !is.na(coef(fit_y1))] %*% fit_y1$coeff[!is.na(coef(fit_y1))]
    yy0 <- designM[, !is.na(coef(fit_y0))] %*% fit_y0$coeff[!is.na(coef(fit_y0))]
  }
  if (identical(outcome_type,"binary")){
    yy1 <- plogis( designM[, !is.na(coef(fit_y1))] %*% coef(fit_y1)[!is.na(coef(fit_y1))] )
    yy0 <- plogis( designM[, !is.na(coef(fit_y0))] %*% coef(fit_y0)[!is.na(coef(fit_y0))] )
  }

  # compute the difference of the mean of two potential outcomes
  trueVal <- mean(yy1, na.rm = TRUE) - mean(yy0, na.rm = TRUE)
  return(trueVal)
}




############################################################################################################
#------------------------------------------------ Stage 5 --------------------------------------------------
############################################################################################################
# bound x with min and max values of bounds
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

# this function is to covert categorical variables to corresponding dichotomous variables.
# depend on the data, there are two choices: remove first dummy or remove most frequent dummy
.convertCategor <- function(W, gform, Qform = NULL, remove_first_dummy = FALSE, remove_most_frequent_dummy = FALSE){
  if (any(lapply(W, is.factor) == TRUE)){

    cat_names <- names(W)[lapply(W, is.factor) == TRUE]
    W_dums = fastDummies::dummy_cols(W, remove_first_dummy = remove_first_dummy,
                                     remove_most_frequent_dummy=remove_most_frequent_dummy)
    W_dum = W_dums[, !names(W_dums) %in% c(cat_names)]

    # extract baseline variable names from gform
    noSpace.gform <- stringr::str_replace_all(string=gform, pattern=" ", repl="")
    g_val <- chartr(old = "+", new = " ", substring(noSpace.gform, 3))
    g_val <- unlist(strsplit( g_val, split = " "))

    gform_cat <- stringr::str_c(c(unlist(sapply(1:length(g_val), function(x)
      names(W_dum)[stringr::str_detect(names(W_dum), g_val[x])]))), collapse = "+")
    gform <- paste0("A~", gform_cat)

    if (!is.null(Qform)){
      noSpace.Qform <- stringr::str_replace_all(string=Qform, pattern=" ", repl="")
      Q_val <- chartr(old = "+", new = " ", substring(noSpace.Qform, 5))
      Q_val <- unlist(strsplit( Q_val, split = " "))

      Qform_cat <- stringr::str_c(c(unlist(sapply(1:length(Q_val), function(x)
        names(W_dum)[stringr::str_detect(names(W_dum), Q_val[x])]))), collapse = "+")
      Qform <- paste0("Y~A+", Qform_cat)
    }
    return(list(W=W_dum, gform=gform, Qform=Qform))
  } else {
    return(list(W=W, gform=gform, Qform=Qform))
  }
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


# this function is to produce the estimates of ps by glm (gform) or SL depending on value of gGLM
.weight_glm_sl <- function(Y, A, W, gform, gbound, gGLM, SL.library){
  if (gGLM){
    # w <- .weight.glm(data.frame(Y, A, W), gform, gbound)
    ObsData <- data.frame(Y, A, W)
    ps <- glm(gform, data = ObsData, family = "binomial")
    designM <- model.matrix (as.formula(gform), ObsData)
    ps.pred <- plogis( designM[, !is.na(coef(ps))] %*% coef(ps)[!is.na(coef(ps))] )
    ps.obs <- ifelse(ObsData$A == 0, 1- ps.pred, ps.pred)
  } else {
  # library(SuperLearner)
    ps.SL <- SuperLearner(Y = A, X = W, family = binomial(),SL.library = SL.library)$SL.predict
    ps.obs <- ifelse(data.frame(Y, A, W)$A == 0, 1- ps.SL, ps.SL)
  }
  w <- .weight(ps.obs, gbound)# use the weight function
  return(w)
}


#------------------------ main function of TMLE to estimate ATE with flexible g (SL or GLM) ----------------------

# input: one observed dataset with outcome- Y(binary and/or continuous),binary treatment A and potential confounders W
# gform -- optional regression formula for estimating P(A=1|W)
# Qform -- regression formula for estimating P(Y=1|A, W)
# gbound -- bound on P(A=1|W), defaults to 0.025
# gGLM -- logical: if TRUE, use GLM; otherwise use SL
# User could define gform in GLM, default is A~W
# SL.library -- prediction algorithms for data adaptive estimation of g
# returns: the estimated ate, sd and 95% confidence interval

TMLE_ate <- function (ObsData, Qform=NULL, outcome_type, gGLM, gbound = 0.025, gform = NULL, SL.library = NULL,
                      gMeth = FALSE, remove_first_dummy = FALSE, remove_most_frequent_dummy = FALSE){

  Y <- ObsData$Y
  A <- ObsData$A
  W <- ObsData[, !names(ObsData) %in% c("Y", "A")]
  n <- length(Y)

  if (gMeth){.check_gMeth(gGLM, gform, SL.library)} # if gMeth=TRUE, show the status message
  if (is.null(Qform)){Qform <- paste ("Y~A+", paste(colnames(W), collapse = "+"))}
  if(is.null(gform)){gform=paste("A~", paste(colnames(W), collapse = "+"))}

  coverCate <- .convertCategor(W, gform, Qform, remove_first_dummy, remove_most_frequent_dummy)
  W = coverCate$W; Qform = coverCate$Qform; gform = coverCate$gform

  # ----------------------------------------------------------- estimate ps
  w <- .weight_glm_sl(Y, A, W, gform, gbound, gGLM, SL.library)

  # ----------------------------------------------------------- estimate Q
  if (identical(outcome_type,"continuous")){
    ab <- range(Y); Y.t <- (Y-ab[1])/diff(ab)
    Ybound <- 0.001 # make value of ystar bounded away from 0 and 1
    Y = .bound(Y.t, c(Ybound, 1-Ybound))
  }
  reg <- suppressWarnings(glm(Qform, data = data.frame(Y, A, W), family = "binomial"))

  designM1 <- model.matrix (as.formula(Qform), data.frame(A=1, W))
  designM0 <- model.matrix (as.formula(Qform), data.frame(A=0, W))
  Q1m <- plogis(designM1[, !is.na(coef(reg))] %*% reg$coeff[!is.na(coef(reg))])
  Q0m <- plogis(designM0[, !is.na(coef(reg))] %*% reg$coeff[!is.na(coef(reg))])

  # TMLE algorithm
  update1 <- suppressWarnings(glm(Y~-1+w, subset=(A==1), offset = qlogis(Q1m), family="binomial", data = data.frame(Y, A, W)))
  update0 <- suppressWarnings(glm(Y~-1+w, subset=(A==0), offset = qlogis(Q0m), family="binomial", data = data.frame(Y, A, W)))
  Q1m_updated <- plogis(qlogis(Q1m) + as.numeric(coef(update1))*w)
  Q0m_updated <- plogis(qlogis(Q0m) + as.numeric(coef(update0))*w)

  if (identical(outcome_type,"continuous")){
    res_ate <- mean((Q1m_updated - Q0m_updated)*diff(ab), na.rm = TRUE)
    IF <- ((Y-Q1m_updated)*A*w + Q1m_updated)-((Y - Q0m_updated)*(1-A)*w + Q0m_updated)-mean((Q1m_updated - Q0m_updated), na.rm = TRUE) # compute the influence function
    sd <- sqrt( sum(IF^2)/(n^2))*diff(ab)
  }
  if (identical(outcome_type,"binary")){
    res_ate <- mean((Q1m_updated - Q0m_updated), na.rm = TRUE)
    IF <- ((Y-Q1m_updated)*A*w + Q1m_updated)-((Y - Q0m_updated)*(1-A)*w + Q0m_updated)-res_ate # compute the influence function
    sd <- sqrt( sum(IF^2)/(n^2))
  }

  ci <- res_ate + c(-1,1)*1.96*(sd) # compute the standard error then the 95% CI
  eff_tmle <- list(ATE = res_ate, SD = sd, CI = ci)
  return(eff_tmle)# return list type values: ate and ci
}



#-------------------------- main function of AIPTW to estimate ATE with flexible g (SL or GLM) -----------------------

# input: one observed dataset with outcome- Y(binary and/or continuous),binary treatment A and potential confounders W
# gform -- optional regression formula for estimating the  P(A=1|W)
# gbound -- bound on P(A=1|W), defaults to 0.025
# gGLM -- logical: if TRUE, use GLM; otherwise use SL
# User could define gform in GLM, default is A~W
# SL.library -- prediction algorithms for data adaptive estimation of g
# returns: the estimated ate, sd and 95% confidence interval
#
AIPTW_ate <- function (ObsData, Qform = NULL, outcome_type, gGLM, gbound = 0.025, gform = NULL, SL.library  = NULL,
                       gMeth = FALSE, remove_first_dummy = FALSE, remove_most_frequent_dummy = FALSE){

  Y <- ObsData$Y
  A <- ObsData$A
  W <- ObsData[, !names(ObsData) %in% c("Y", "A")]
  n <- length(Y)

  if (gMeth){.check_gMeth(gGLM, gform, SL.library)} # if gMeth=TRUE, show the status message
  if (is.null(Qform)){Qform <- paste ("Y~A+", paste(colnames(W), collapse = "+"))}
  if(is.null(gform)){gform=paste("A~", paste(colnames(W), collapse = "+"))}

  coverCate <- .convertCategor(W, gform, Qform, remove_first_dummy, remove_most_frequent_dummy)
  W = coverCate$W; Qform = coverCate$Qform; gform = coverCate$gform

  # get the potential outcome based on the coeffs of glm on original data
  designM1 <- model.matrix (as.formula(Qform), data.frame(A=1, W))
  designM0 <- model.matrix (as.formula(Qform), data.frame(A=0, W))

  if (identical(outcome_type,"continuous")){
    reg <- glm(Qform, data = data.frame(Y, A, W), family = "gaussian")
    Q1m <- designM1[, !is.na(coef(reg))] %*% reg$coeff[!is.na(coef(reg))]
    Q0m <- designM0[, !is.na(coef(reg))] %*% reg$coeff[!is.na(coef(reg))]
  }

  if (identical(outcome_type,"binary")){
    reg <- glm(Qform, data = data.frame(Y, A, W), family = "binomial")
    Q1m <- plogis(designM1[, !is.na(coef(reg))] %*% reg$coeff[!is.na(coef(reg))])
    Q0m <- plogis(designM0[, !is.na(coef(reg))] %*% reg$coeff[!is.na(coef(reg))])
  }

  # estimate propensity score
  w <- .weight_glm_sl(Y, A, W, gform, gbound, gGLM, SL.library)

  # estimate ate
  mu1 <- sum((Y - Q1m)*A*w)/n + mean(Q1m); mu0 <- sum((Y - Q0m)*(1-A)*w)/n + mean(Q0m)
  res_ate <- mu1 - mu0
  IF <- ((Y-Q1m)*A*w + Q1m)-((Y - Q0m)*(1-A)*w + Q0m) - res_ate # compute the influence function
  sd <- sqrt( sum(IF^2)/(n^2))
  ci <- res_ate + c(-1,1)*1.96*sd # compute the standard error then the 95% CI

  eff_aiptw <- list(ATE = res_ate, SD = sd, CI = ci)
  return(eff_aiptw)# return list type values: ate and ci
}



#-------------------------- main function of IPTW to estimate ATE with flexible g (SL or GLM) -----------------------

# input: one observed dataset with outcome- Y(binary and/or continuous),binary treatment A and potential confounders W
# outcome_type -- "continuous" or "binary"
# gform -- optional regression formula for estimating the  P(A=1|W)
# gbound -- bound on P(A=1|W), defaults to 0.025
# gGLM -- logical: if TRUE, use GLM; otherwise use SL
# User could define gform in GLM, default is A~W
# SL.library -- prediction algorithms for data adaptive estimation of g
# returns: the estimated ate, sd and 95% confidence interval

# IPTW_ate <- function (ObsData, outcome_type, gGLM, gbound = 0.025, gform= NULL, SL.library = NULL,
#                       gMeth = FALSE, remove_first_dummy = FALSE, remove_most_frequent_dummy = FALSE){
#   Y <- ObsData$Y
#   A <- ObsData$A
#   W <- ObsData[, !names(ObsData) %in% c("Y", "A")]
#
#   if (gMeth){.check_gMeth(gGLM, gform, SL.library)} # if gMeth=TRUE, show the status message
#   if(is.null(gform)){gform=paste("A~", paste(colnames(W), collapse = "+"))}
#   coverCate <- .convertCategor(W, gform, remove_first_dummy, remove_most_frequent_dummy)
#   W = coverCate$W; gform = coverCate$gform
#
#   # estimate propensity score
#   w <- .weight_glm_sl(Y, A, W, gform, gbound, gGLM, SL.library)
#
#   mod <- lm(Y~ A, data = data.frame(Y, A, W), weights = w)
#   res_ate <-summary(mod)$coef[2,1]
#   # compute 95% CI
#   sd <- sqrt(sandwich::vcovHC(mod)[2,2])
#   ci <- res_ate + 1.96*c(-1,1)*sd
#   eff_iptw <- list(ATE = res_ate, SD = sd, CI = ci)
#   return(eff_iptw)# return list type values: ate, CI
# }



############################################################################################################
#------------------------------------------------ Stage 6 --------------------------------------------------
############################################################################################################

#-------------- helper function to check all input variables --------------------
.check_var2 <- function(Y, A, W, outcome_type, M, gbound, gGLM, gform, SL.library){
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
  if (class(gGLM)!="logical"){stop("gGLM must be TRUE or FALSE")}

  ######### check gform ########
  # remove the space in the string of input gform
  if (!is.null(gform)){
    noSpace.gform <- stringr::str_replace_all(string = gform, pattern = " ", repl= "")
    if (substr(noSpace.gform, 1, 2)!= "A~"){stop("Error inputs for gform occur. ")}
  }
}



#----------------- main function to compute the bias of the estimates on IPTW,AIPTW,TMLE ------------------


# input: outcome- Y(binary and/or continuous); binary treatment A; a set of potential confounders- W
# outcome_type -- "continuous" or "binary"
# M is the number of replications
# gform -- optional regression formula for estimating P(A=1|W), default is A~W
# Qform -- regression formula for estimating P(Y=1|A, W), default is Y~A,W
# gbound -- bound on P(A=1|W), defaults to 0.025
# gGLM -- logical: if TRUE, use GLM; otherwise use SL
# SL.library -- prediction algorithms for data adaptive estimation of g
# returns: the estimated ate and 95% confidence intervals for three methods: IPTW, AIPTW, TMLE and pooled ps.

bdt <- function (Y, A, W, outcome_type, M, gbound = 0.025, gGLM, Qform = NULL, gform = NULL, SL.library = NULL,
                remove_first_dummy = FALSE, remove_most_frequent_dummy = FALSE){

  bias_TMLE <- NULL; bias_AIPTW <- NULL # bias_IPTW <- NULL
  sd_TMLE <- NULL; sd_AIPTW <- NULL
  cov_TMLE <- cov_AIPTW <- 0 # cov_IPTW <- 0
  ps.all.GLM <- NULL; ps.all.SL <- NULL

  # check all inputs
  .check_var2(Y, A, W, outcome_type, M, gbound, gGLM, gform, SL.library)
  .check_gMeth(gGLM, gform, SL.library)

  # create the outcome regression formula
  if (is.null(Qform)){Qform <- paste ("Y~A+", paste(colnames(W), collapse = "+"))}
  Qform <- stringr::str_replace_all(string = Qform, pattern = " ", repl = "")
  if(is.null(gform)){gform=paste("A~", paste(colnames(W), collapse = "+"))}

  # test factors in W and change to dummy variables
  coverCate <- .convertCategor(W, gform, Qform, remove_first_dummy, remove_most_frequent_dummy)
  W_dumm = coverCate$W; Qform_dumm = coverCate$Qform; gform_dumm = coverCate$gform

  # compute the coeffs derived from the glm of outcome for original dataset
  fit_y0 <- .fit_glm_Y(Y, A, W = W_dumm, Qform = Qform_dumm, outcome_type = outcome_type)$fit_y0
  fit_y1 <- .fit_glm_Y(Y, A, W = W_dumm, Qform = Qform_dumm, outcome_type = outcome_type)$fit_y1

  true_eff <- .trueVal(Y, A, W = W_dumm, Qform = Qform_dumm, outcome_type = outcome_type, fit_y1= fit_y1, fit_y0=fit_y0)

  # loops to create N bootstrap datasets
  for (i in 1:M){
    bootdata <- .sim_data( A, W = W_dumm, Qform = Qform_dumm, outcome_type = outcome_type, fit_y1= fit_y1, fit_y0=fit_y0)

    if (gGLM){
      ps.GLM <- data.frame(.ps_GLM(A = bootdata$A, W = bootdata[, !names(bootdata) %in% c("Y", "A")],
                                   gform = gform_dumm, gbound = gbound,
                                   remove_first_dummy, remove_most_frequent_dummy)$ps_glm)
      names(ps.GLM) <- c("A", "ps")
      ps.all.GLM <- rbind(ps.all.GLM, ps.GLM)
    }
    if (!is.null(SL.library) & !gGLM){
      ps.SL <- data.frame(.ps_SL(A= bootdata$A, W = bootdata[, !names(bootdata) %in% c("Y", "A")], SL.library = SL.library, gbound = gbound)$ps_sl)
      names(ps.SL) <- c("A", "ps")
      ps.all.SL <- rbind(ps.all.SL, ps.SL)
    }

    # compute estimates of tmle, iptw, and aiptw of the bootstrap dataset
    tmle_res <- TMLE_ate(bootdata, Qform = Qform_dumm, outcome_type = outcome_type, gGLM, gbound = gbound,
                         gform = gform_dumm, SL.library = SL.library, gMeth = FALSE,
                         remove_first_dummy, remove_most_frequent_dummy)

    aiptw_res <- AIPTW_ate(bootdata, Qform = Qform_dumm, outcome_type = outcome_type, gGLM, gbound = gbound,
                           gform = gform_dumm, SL.library = SL.library, gMeth = FALSE,
                           remove_first_dummy, remove_most_frequent_dummy)

    # iptw_res <- IPTW_ate(bootdata, outcome_type = outcome_type, gGLM, gbound = gbound, gform = gform_dumm,
    #                      SL.library = SL.library, gMeth = FALSE,
    #                      remove_first_dummy, remove_most_frequent_dummy)

    # compute the bias
    bias_TMLE[i] <- tmle_res$ATE-true_eff
    # bias_IPTW[i] <- iptw_res$ATE-true_eff
    bias_AIPTW[i] <- aiptw_res$ATE-true_eff

    sd_TMLE[i] <- tmle_res$SD
    sd_AIPTW[i] <- aiptw_res$SD

    # compute the times that true value is within the intervals of three methods
    cov_TMLE <- cov_TMLE + ifelse(true_eff >= tmle_res$CI[1] & true_eff <= tmle_res$CI[2], 1, 0 )
    # cov_IPTW <- cov_IPTW + ifelse(true_eff >= iptw_res$CI[1] & true_eff <= iptw_res$CI[2], 1, 0 )
    cov_AIPTW <- cov_AIPTW + ifelse(true_eff >= aiptw_res$CI[1] & true_eff <= aiptw_res$CI[2], 1, 0 )
  }

  results <- list(true_Effect = true_eff, gbound =  gbound, #bias_IPTW = bias_IPTW,
                  bias_AIPTW = bias_AIPTW, bias_TMLE = bias_TMLE,#cov_IPTW = cov_IPTW/M,
                  sd_AIPTW = sd_AIPTW, sd_TMLE = sd_TMLE,#cov_IPTW = cov_IPTW/M,
                  cov_AIPTW = cov_AIPTW/M, cov_TMLE = cov_TMLE/M,
                  ps_GLM = ps.all.GLM, ps_SL = ps.all.SL)
  class(results) <- "bdt"
  return(results)
}



#----------------- print functions to output the estimated ate bias of AIPTW, IPTW, TMLE -------------------
# input the results of main bias_boot function
# return the bias of the three estimated ate
bias.bdt <- function(object){
  if(identical(class(object), "bdt")){
    # cat("Eestimated bias of ATE with IPTW, AIPTW, TMLE:\n\n")
    cat("Eestimated bias of ATE with AIPTW and TMLE:\n\n")

    bias.ate <- data.frame(cbind(#object$bias_IPTW,
                        object$bias_AIPTW, object$bias_TMLE))
    # colnames(bias.ate) <- c("bias_IPTW", "bias_AIPTW", "bias_TMLE")
    colnames(bias.ate) <- c("bias_AIPTW", "bias_TMLE")
  } else {
    stop("Object must have class 'bdt'")}
  return(bias.ate)
}
bias <- function(object){ UseMethod("bias",object)}




#----------------- summary of BDT -------------------
# create BDT summary object
# input the object of main bias_boot function
# return the min, max, mean, median of bias of the three estimated ate, coverage rates and ps

summary.bdt <- function(object,...){
  if(identical(class(object), "bdt")){

    if (!is.null(object$ps_GLM)){
      ps1.sum = summary(object$ps_GLM$ps)
      ps0.sum = summary(1-object$ps_GLM$ps)
      smda = data.frame(object$ps_GLM)
      ps11.sum = summary(smda[smda$A==1,]$ps)
      ps00.sum = summary(1-smda[smda$A==0,]$ps)}
    else {
      ps1.sum = summary(object$ps_SL$ps)
      ps0.sum = summary(1-object$ps_SL$ps)
      smda = data.frame(object$ps_SL)
      ps11.sum = summary(smda[smda$A==1,]$ps)
      ps00.sum = summary(1-smda[smda$A==0,]$ps)}

    bias <- data.frame(rbind(#summary(object$bias_IPTW),
                          summary(object$bias_AIPTW), summary(object$bias_TMLE)))
    dimnames(bias) <- list(c(#"Bias of IPTW",
                              "Bias of AIPTW", "Bias of TMLE"),
                           c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."))
    # coverage <- paste0("Coverage rate of IPTW: ", round(object$cov_IPTW,4),
    #                   ";  Coverage rate of AIPTW: ",round(object$cov_AIPTW,4),
    #                   ";  Coverage rate of TMLE: ",round(object$cov_TMLE,4))
    sd <- paste0("Standard error of AIPTW: ",round(mean(object$sd_AIPTW),4),
                 ";  standard error of TMLE: ",round(mean(object$sd_TMLE),4))
    coverage <- paste0("Coverage rate of AIPTW: ",round(object$cov_AIPTW,4),
                       ";  Coverage rate of TMLE: ",round(object$cov_TMLE,4))
    ps <- data.frame(rbind(ps1.sum, ps0.sum, ps11.sum, ps00.sum))
    dimnames(ps) <- list(c("P(A=1|X) for all subjects", "P(A=0|X) for all subjects",
                      "P(A=1|X) for subgroups A=1", "P(A=0|X) for subgroups A=0"),
                      c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max."))
    summary.bdt <- list(true.effect = object$true_Effect, gbound = object$gbound,
                      #bias.iptw = summary(object$bias_IPTW),
                      bias.aiptw = summary(object$bias_AIPTW), bias.tmle = summary(object$bias_TMLE),
                      #cov.iptw = round(object$cov_IPTW,4),
                      sd.aiptw = mean(object$sd_AIPTW), sd.tmle = mean(object$sd_TMLE),
                      cov.aiptw = round(object$cov_AIPTW,4), cov.tmle = round(object$cov_TMLE,4),
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

    # cat("\nAverage bias and coverage rates of IPTW, AIPTW and TMLE estimators:\n\n")
    cat("\nMean and median bias, standard errors and coverage rates of AIPTW and TMLE estimators:\n\n")
    #i <- specify_decimal(x$bias.iptw[4], 5);
    a <- specify_decimal(x$bias.aiptw[4], 5); t <- specify_decimal(x$bias.tmle[4], 5)
    #s <- specify_decimal(x$bias.iptw[3], 5);
    u <- specify_decimal(x$bias.aiptw[3], 5); v <- specify_decimal(x$bias.tmle[3], 5)
    m <- specify_decimal(x$sd.aiptw, 5); n <- specify_decimal(x$sd.tmle, 5)
    # res_data <- data.frame(c(" IPTW","AIPTW"," TMLE"), c(i, a, t), c(s, u, v), c(x$cov.iptw, x$cov.aiptw, x$cov.tmle))
    res_data <- data.frame(c("AIPTW"," TMLE"), c(a, t), c(u, v), c(m, n), c(x$cov.aiptw, x$cov.tmle))
    names(res_data) <- c("     ","MeanBias",  "Med.Bias", "   SD ", " Cov.")
    gdata::write.fwf(res_data)

    cat("\n\nSummary of propensity scores truncated at gbound", x$gbound, " \n\n")
    print(x$ps)
  }
}



#--------------------------- plots of the results -------------------------------
# input the results of main bias_boot function
# return two plots:
#   1. the boxplots and points of bias of the three estimated ate and coverage rates
#   2. the density plots of pooled propensity scores
# For the boxplot: y-axis represents the bias of estimates; x-axis represents estimates: IPTW, AIPTW, TMLE
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
    # type_est <- c(rep("IPTW",N), rep("AIPTW",N),rep("TMLE",N) )
    # Bias_ATE <- c(x$bias_IPTW, x$bias_AIPTW, x$bias_TMLE)
    type_est <- c(rep("AIPTW",N),rep("TMLE",N) )
    Bias_ATE <- c(x$bias_AIPTW, x$bias_TMLE)
    est_data <- data.frame(type_est, Bias_ATE)
    colnames(est_data) <- c("type_est", "Bias_ATE")
    mt = ifelse(is.null(x$ps_GLM), "SL for g", "GLM for g")
    title <- paste0("Boxplots of the ATE bias with ", x$gbound, " gbound under ", mt)
    # type = c("IPTW", "AIPTW", "TMLE")
    type = c("AIPTW", "TMLE")
    high = min(est_data$Bias_ATE)-0.03
    # coverage = round(c(x$cov_IPTW, x$cov_AIPTW, x$cov_TMLE), 2)
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
    if (!is.null(object$ps_GLM) ){
      n <- length(object$ps_GLM$ps)
      smda <- data.frame(object$ps_GLM)
      ps1 <- data.frame(as.factor(c(rep("GLM", n))),c(log(1/smda$ps)))
      ps0 <- data.frame(as.factor(c(rep("GLM", n))),c(log(1/(1-smda$ps))))
      ps.A1 <- data.frame(as.factor(c(rep("GLM", length(smda[smda$A ==1,]$ps)))),c(log(1/smda[smda$A ==1,]$ps)))
      ps.A0 <- data.frame(as.factor(c(rep("GLM", length(smda[smda$A ==0,]$ps)))),c(log(1/(1-smda[smda$A ==0,]$ps))))
      colnames(ps1) <- colnames(ps0) <- colnames(ps.A1) <- colnames(ps.A0) <- c("GLM", "Propensity")

      ps.datas <- list(ps.A1, ps.A0, ps1, ps0)
      xlab_name <- paste0(c("log(1/g_1)[A=1]", "log(1/g_1)[A=0]", "log(1/g_1)", "log(1/g_0)"), " on GLM")

      ps_plot <- sapply(1:4, function(y) list(ggplot(ps.datas[[y]], aes(fill = ps.datas[[y]]$GLM, x = ps.datas[[y]]$Propensity))+
                                                geom_density(alpha = 0.3) + theme(axis.text.x = element_text(angle = 0))+
                                                xlab(xlab_name[y]) + theme(legend.position = "none")))

      ps_plots <- ggpubr::ggarrange(ps_plot[[1]], ps_plot[[2]], ps_plot[[3]], ps_plot[[4]], labels = c("A","B","C","D"), ncol = 2, nrow = 2)
    }

    if (!is.null(object$ps_SL) ){
      n <- length(object$ps_SL$ps)
      smda <- data.frame(object$ps_SL)
      ps1 <- data.frame(as.factor(c(rep("SL", n))), c(log(1/smda$ps)))
      ps0 <- data.frame(as.factor(c(rep("SL", n))), c(log(1/(1-smda$ps))))
      ps.A1 <- data.frame(as.factor(c(rep("SL", length(smda[smda$A ==1,]$ps)))), c(log(1/smda[smda$A ==1,]$ps)))
      ps.A0 <- data.frame(as.factor(c(rep("SL", length(smda[smda$A ==0,]$ps)))), c(log(1/(1-smda[smda$A ==0,]$ps))))
      colnames(ps1) <- colnames(ps0) <- colnames(ps.A1) <- colnames(ps.A0) <- c("SL", "Propensity")

      ps.datas <- list(ps.A1, ps.A0, ps1, ps0)
      xlab_name <- paste0(c("log(1/g_1)[A=1]", "log(1/g_1)[A=0]", "log(1/g_1)", "log(1/g_0)"), " on SL")

      ps_plot <- sapply(1:4, function(y) list(ggplot(ps.datas[[y]], aes(fill = ps.datas[[y]]$SL, x = ps.datas[[y]]$Propensity)) +
                                                geom_density(alpha = 0.3) + theme(axis.text.x = element_text(angle = 0)) +
                                                xlab(xlab_name[y]) + theme(legend.position = "none")))

      ps_plots <- ggpubr::ggarrange(ps_plot[[1]], ps_plot[[2]], ps_plot[[3]], ps_plot[[4]], labels = c("A","B","C","D"), ncol = 2, nrow = 2)
    }
    fig_box_density <- list(boxplot_bias_ATE = est_plot, densityplot_ps = ps_plots)
  }
  else {
    stop("Object must have class 'bdt'")
    fig_box_density <- NULL
  }
  par(ask=TRUE)
  return(fig_box_density)
}

