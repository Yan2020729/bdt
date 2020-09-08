
#---------------------------------------------- Description ------------------------------------------
# This R code implements adapted bootstrap diagnostic tool (BDT), originally described in Bahamyirou A, Blais L,
# Forget A, Schnitzer ME. "Understanding and diagnosing the potential for bias when using machine learning
# methods with doubly robust causal estimators."Statistical Methods in Medical Research. 2019 Jun;28(6):1637-50.

# This bootstrap diagnostic tool (BDT) is to identify the instability in the estimation of marginal causal effects
# in the presence of near practical positivity violations. Data-adaptive methods can produce a separation of the
# two exposure groups in terms of propensity score densities which can lead to biased finite-sample estimates
# of the treatment effect. BDT is based on a bootstrap resampling of the subjects and simulation of the outcome
# data conditional on observed treatment and covariates in order to diagnose whether the estimation using data-adaptive
# methods for the propensity score may lead to large bias and poor coverage in the finite sample.

# The parameter of interest in this program is the additive effect of a binary point treatment- A, on an outcome- Y
# (binary and/or continuous), adjusting for a set of potential confounders- W based on observing n i.i.d. observations
# on (W,A,Y). Missing data are not allowed.
#
# Super Learner (SL) [1] is used for the data-adaptive estimation of the propensity score. This package calls the function
# SuperLearner in the package SuperLearner [2]. This package also depends on the implementation of Targeted maximum likelihood
# estimation (TMLE) [3] in the package tmle [4]. Inverse probability of treatment weighting (IPTW) and augmented inverse
# probability of treatment weighting (AIPTW) [5] are also implemented in this package.
#
# There are six stages:
# Stage 1:  For a given observed dataset (A, W), use GLM and/or SL to estimate the propensity score and plot the distribution for treatment A=1 and A=0 for all subjects and for subgroups A=1, A=0.
# Stage 2:  For an observed dataset with n subjects, fit a GLM of Y on W for each of the two subgroups with A=1,0. The outcome regression only includes the main terms of W.
# Stage 3:  Create bootstrap datasets by resampling the n subjects with replacement and removing the observed outcome values. The new outcome realizations are predicted using the estimated coefficients from the 2nd Stage.
# Stage 4:  Compute the "true" effect of the bootstrap datasets. We achieve this by computing the difference of the two potential outcomes using the estimated coefficients from the 2nd stage.
# Stage 5:  Use three different approaches to estimate the additive effect where the outcome model is "true" GLM (specified in Stage 2) and the propensity score is modeled by SL and/or logistic regression of A on W as specified by the user:
#           1) TMLE, 2) IPTW, and 3) AIPTW.
#           The user can define the regression formula for the propensity in GLM or the SL functions in SL.
#           There is no default SL function.
#           Default GLM formula: A~W.
#           Default bound on the propensity score is 0.025.
# Stage 6:  Compute the bias of the estimates and the coverage rates then give a summary of the results over the bootstraps (numerical and boxplot).
#           The output of the main function includes each estimated propensity score and  the corresponding density plots for the combined treatment groups and separately.
#
# References:
# [1] Van der Laan, M.J., Polley, E.C. and Hubbard, A.E., 2007. Super learner. Statistical applications in genetics and molecular biology, 6(1).
# [2] Polley E, LeDell E, Kennedy C, Lendle S, van der Laan M. 2019. Package 'SuperLearner'.
# [3] Van Der Laan MJ, Rubin D. 2006. Targeted maximum likelihood learning. The international journal of biostatistics. 28;2(1).
# [4] Gruber, S., Van der Laan, M. J. 2011. tmle: An R package for targeted maximum likelihood estimation.
# [5] van der Laan, M.J. and Luedtke, A.R., 2014. Targeted learning of an optimal dynamic treatment, and statistical inference for its mean outcome.



#-----------------------------------------------------------------------------------------------------
# rm(list=ls())


############################################################################################################
#------------------------------------------------ Stage 1 --------------------------------------------------
############################################################################################################

# bound x with min and max values of bounds
.bound <- function(x, bounds){
  x[x<min(bounds)] <- min(bounds); x[x>max(bounds)] <- max(bounds)
  return(x)
}


#--------------------------- helper function to produce ps via GLM -----------------------------
# input: exposure-- A(binary); set of confounders--W
#        gform-- optional regression formula (default: A~W)
# return ps.pred-- p(A=1|W) for all subjects
#        ps.A1.pred-- p(A=1|W) for subgroups with A=1
#        ps.A0.pred-- p(A=1|W) for subgroups with A=0

.ps_GLM <- function(A, W, gform, gbound){

  ps <- glm(gform, data = data.frame(A, W), family = "binomial")
  designM <- model.matrix(as.formula(gform), data.frame(A, W))
  ps.pred <- plogis( designM[, !is.na(ps$coeff)] %*% ps$coeff[!is.na(ps$coeff)] )
  ps.pred <- .bound(ps.pred, c(gbound, 1-gbound))
  ps_glm <- data.frame(A, ps.pred)
  return(list(ps_glm = ps_glm, fit_glm = summary(ps)))
}


#--------------------------- helper function to produce ps via SL -----------------------------
# input: exposure-- A(binary); set of confounders--W
#        SL.library-- optional SL functions (no default)
# return ps.pred-- p(A=1|W) for all subjects
#        ps.A1.pred-- p(A=1|W) for subgroups with A=1
#        ps.A0.pred-- p(A=1|W) for subgroups with A=0

.ps_SL <- function(A, W, SL.library, gbound){
  # use SL to compute p(A=1|W)
  library(SuperLearner)
  ps <- SuperLearner(Y = A, X = W, family = binomial(), SL.library = SL.library)
  ps.pred <- .bound(ps$SL.predict, c(gbound, 1-gbound))
  ps_sl <- data.frame(A, ps.pred)
  return(list(ps_sl = ps_sl, fit_SL = ps))
}




#-------------- helper function to check all input variables --------------------
.check_var1 <- function(A, W, gform1, gform2, SL.library1, SL.library2, gbound){
  # check the dimension of W and A
  if(length(A)!= nrow(W) ){stop("The data dimensions do not correspond.")}

  # check missing values
  miss <- sapply(data.frame(A, W), function(x)sum(is.na(x)))
  if (sum(miss) > 0) {stop("There are missing values in the dataset. ")}

  # check gbound
  if (gbound>1 || gbound < 0){stop("gbound must be between 0 and 1.")}

  # check binary treatment A
  if (is.factor(A)){if (length(levels(A)) != 2L) stop("A must be a binary variable. ")}
  if (is.numeric(A)) {
    uniq_a <- unique(A)
    if ( length(unique(uniq_a)) > 2L || any(as.integer(uniq_a) != uniq_a) ) {stop("A must be a binary variable. ")}
    if (any(! uniq_a %in% c(0,1))) stop("A must be defined as 1 for treatment and 0 for control.")
  }
  if (!is.numeric(A) & !is.factor(A) & !is.logical(A)) {stop("A must be a binary variable. ")}

  # check gform1 and gform2
  if (!is.null(gform1)){
    noSpace.gform <- stringr::str_replace_all(string = gform1, pattern = " ", repl = "")
    if (substr(noSpace.gform, 1, 2)!= "A~"){stop("Error inputs for gform1 occur. ")}
  }
  if (!is.null(gform2)){
    noSpace.gform <- stringr::str_replace_all(string = gform2, pattern = " ", repl = "")
    if (substr(noSpace.gform, 1, 2)!= "A~"){stop("Error inputs for gform2 occur. ")}
  }

  if (is.null(gform1) & is.null(SL.library1) & is.null(gform2) & is.null(SL.library2))
    stop("Please specify the gforms and/or SL.librarys. ")
}



#--------------------------- main Summary function of ps for all subjects -----------------------------

# input: Y, A, W, gGLM, gorm(default is NULL), SL.library(no default value)
# return: summary of p(A=1|W) for all subjects with GLM and SL

summary_ps <- function(A, W, gform1 = NULL, gform2 = NULL, SL.library1 = NULL, SL.library2 = NULL, gbound = 0, verbose = FALSE){
  .check_var1(A, W, gform1, gform2, SL.library1, SL.library2, gbound)
  # check verbose
  if (!is.logical(verbose)) {stop("Verbose should be TRUE or FALSE. ")}

  # always have gform1 specified if gform2 is specified
  if (is.null(gform1) & !is.null(gform2)) {gform1 = gform2; gform2 = NULL}

  # always have SL.library1 specified if SL.library2 is specified
  if (is.null(SL.library1) & !is.null(SL.library2)) {SL.library1 = SL.library2; SL.library2 = NULL}

  # create a case vector
  case_vec <- sapply(1:4, function(x)as.numeric(!is.null(list(gform1, gform2, SL.library1, SL.library2)[[x]])))

  # create summary names
  rowName <- c("P(A=1|W) for all subjects by ", "P(A=0|W) for all subjects by ", "P(A=1|W) in subgroups A=1 by ", "P(A=0|W) in subgroups A=0 by ")
  colName <- c("Min.", "1st Qu.","Median", "Mean", "3rd Qu.", "Max.")

  # compute ps on case_vec
  ps_glm_list <- sapply(1:2, function(x)if (case_vec[x]==1) {
    ps.GLM <- .ps_GLM(A = A, W = W, gform = eval(as.name(paste0("gform",x))), gbound)$ps_glm
    data_glm <- data.frame(rbind(summary(ps.GLM$ps.pred), summary(1-ps.GLM$ps.pred),
                                 summary(ps.GLM[A ==1,]$ps.pred), summary(1-ps.GLM[A ==0,]$ps.pred)))
    dimnames(data_glm) <- list(c(paste0(rowName,"GLM",x)),colName)
    list(data_glm)})
  if (is.null(ps_glm_list[[2]]) & !is.null(ps_glm_list[[1]])){
    rownames(ps_glm_list[[1]][[1]]) <- gsub('.{1}$', '', rownames(ps_glm_list[[1]][[1]]))}

  ps_SL_list <- sapply(1:2, function(x)if (case_vec[x+2] ==1) {
    ps.SL <- .ps_SL(A = A, W = W, SL.library = eval(as.name(paste0("SL.library",x))), gbound)$ps_sl
    data_SL <- data.frame(rbind(summary(ps.SL$ps.pred), summary(1-ps.SL$ps.pred),
                                summary(ps.SL[A ==1,]$ps.pred), summary(1-ps.SL[A ==0,]$ps.pred)))
    dimnames(data_SL) <- list(c(paste0(rowName,"SL",x)),colName)
    list(data_SL)})
  if (is.null(ps_SL_list[[2]]) & !is.null(ps_SL_list[[1]])){
    rownames(ps_SL_list[[1]][[1]]) <- gsub('.{1}$', '', rownames(ps_SL_list[[1]][[1]]))}

  # create all possible case vectors
  case_all_vec <- list(c(1,0,0,0), c(1,1,0,0), c(1,0,1,0), c(1,1,1,0), c(1,0,1,1), c(1,1,1,1), c(0,0,1,0), c(0,0,1,1))

  data_name <- c("ps_glm_list[[1]][[1]]",
                 "rbind(ps_glm_list[[1]], ps_glm_list[[2]])",
                 "rbind(ps_glm_list[[1]][[1]], ps_SL_list[[1]][[1]])",
                 "rbind(ps_glm_list[[1]], ps_glm_list[[2]], ps_SL_list[[1]][[1]])",
                 "rbind(ps_glm_list[[1]][[1]], ps_SL_list[[1]], ps_SL_list[[2]])",
                 "rbind(ps_glm_list[[1]], ps_glm_list[[2]], ps_SL_list[[1]], ps_SL_list[[2]])",
                 "ps_SL_list[[1]][[1]]",
                 "rbind(ps_SL_list[[1]], ps_SL_list[[2]])")

  index <- as.numeric(unlist(sapply(1:8, function(x) if (all(case_all_vec[[x]] == case_vec)) return(x))))
  ps_all <- eval(parse(text = data_name[index]))


  if (verbose == TRUE) {
    try(fit_glm_list <- sapply(1:2, function(x) if (case_vec[x]==1) {
      fit_GLM <- .ps_GLM(A = A, W = W, gform = eval(as.name(paste0("gform",x))), gbound)$fit_glm
      list(fit_GLM)}))
    try(fit_SL_list <- sapply(1:2, function(x) if (case_vec[x+2] ==1) {
      fit_SL <- .ps_SL(A = A, W = W, SL.library = eval(as.name(paste0("SL.library",x))), gbound)$fit_SL
      list(fit_SL)}))

    summ_name <- c("fit_glm_list[[1]][[1]]",
                   "list(fit_glm_list[[1]], fit_SL_list[[2]])",
                   "list(fit_glm_list[[1]][[1]], fit_SL_list[[1]][[1]])",
                   "list(fit_glm_list[[1]], fit_glm_list[[2]], fit_SL_list[[1]][[1]])",
                   "list(fit_glm_list[[1]][[1]], fit_SL_list[[1]], fit_SL_list[[2]])",
                   "list(fit_glm_list[[1]], fit_glm_list[[2]], fit_SL_list[[1]], fit_SL_list[[2]])",
                   "fit_SL_list[[1]][[1]]",
                   "list(fit_SL_list[[1]], fit_SL_list[[2]])")
    fit_all <- eval(parse(text = summ_name[index]))
    return(list(probabilities = ps_all, summaries = fit_all))
  }
  else {return(ps_all)}
}



############################### main function: density plots of log of weights ###############################
# input: A, W, gorm(default is NULL), SL.library(default is NULL)
# return: 2 x 2 density plots (A, B, C, D) of the log of the estimated weights for
#         A--treatment A=1 in subset of subjects with A=1
#         B--treatment A=0 in subset of subjects with A=0
#         C--treatment A=1 for all subjects
#         D--treatment A=0 for all subjects
plot_ps <- function(A, W, gform1 = NULL, gform2 = NULL, SL.library1 = NULL, SL.library2 = NULL, gbound = 0,...){
  .check_var1(A, W, gform1, gform2, SL.library1, SL.library2, gbound)

  # always have gform1 specified if gform2 is specified
  if (is.null(gform1) & !is.null(gform2)) {gform1 = gform2; gform2 = NULL}

  # always have SL.library1 specified if SL.library2 is specified
  if (is.null(SL.library1) & !is.null(SL.library2)) {SL.library1 = SL.library2; SL.library2 = NULL}

  # create a case vector
  case_vec <- sapply(1:4, function(x)as.numeric(!is.null(list(gform1, gform2, SL.library1, SL.library2)[[x]])))

  # compute ps on case_vec
  ps_glm_list <- sapply(1:2, function(x)if (case_vec[x]==1) {
    ps.GLM <- .ps_GLM(A = A, W = W, gform = eval(as.name(paste0("gform",x))), gbound)$ps_glm
    n.GLM <- length(ps.GLM$ps.pred)
    ps1.GLM <- data.frame(as.factor(c(rep(paste0("GLM",x), n.GLM))),c(log(1/ps.GLM$ps.pred)))
    ps0.GLM <- data.frame(as.factor(c(rep(paste0("GLM",x), n.GLM))),c(log(1/(1-ps.GLM$ps.pred))))
    ps.A1.GLM <- data.frame(as.factor(c(rep(paste0("GLM",x), length(ps.GLM[A ==1,]$ps.pred)))),c(log(1/ps.GLM[A ==1,]$ps.pred)))
    ps.A0.GLM <- data.frame(as.factor(c(rep(paste0("GLM",x), length(ps.GLM[A ==0,]$ps.pred)))),c(log(1/(1-ps.GLM[A ==0,]$ps.pred))))
    colnames(ps1.GLM) <- colnames(ps0.GLM) <- colnames(ps.A1.GLM) <- colnames(ps.A0.GLM) <- c("Method", "Propensity")
    list(list(ps1.GLM, ps0.GLM, ps.A1.GLM, ps.A0.GLM))})
  if (is.null(ps_glm_list[[2]])) {for (i in 1:4)  {ps_glm_list[[1]][[1]][[i]]$Method = gsub('.{1}$', '', ps_glm_list[[1]][[1]][[i]]$Method) }}

  ps_SL_list <- sapply(1:2, function(x)if (case_vec[x+2] ==1) {
    ps.SL <- .ps_SL(A=A, W=W, SL.library = eval(as.name(paste0("SL.library",x))), gbound)$ps_sl
    n.SL <- length(ps.SL$ps.pred)
    ps1.SL <- data.frame(as.factor(c(rep(paste0("SL",x), n.SL))),c(log(1/ps.SL$ps.pred)))
    ps0.SL <- data.frame(as.factor(c(rep(paste0("SL",x), n.SL))),c(log(1/(1-ps.SL$ps.pred))))
    ps.A1.SL <- data.frame(as.factor(c(rep(paste0("SL",x), length(ps.SL[A ==1,]$ps.pred)))),c(log(1/ps.SL[A ==1,]$ps.pred)))
    ps.A0.SL <- data.frame(as.factor(c(rep(paste0("SL",x), length(ps.SL[A ==0,]$ps.pred)))),c(log(1/(1-ps.SL[A ==0,]$ps.pred))))
    colnames(ps1.SL) <- colnames(ps0.SL) <- colnames(ps.A1.SL) <- colnames(ps.A0.SL) <- c("Method", "Propensity")
    list(list(ps1.SL, ps0.SL, ps.A1.SL, ps.A0.SL))})
  if (is.null(ps_SL_list[[2]])) {for (i in 1:4)  {ps_SL_list[[1]][[1]][[i]]$Method = gsub('.{1}$', '', ps_SL_list[[1]][[1]][[i]]$Method) }}

  # adjust the case vect and produce index number among all case vector
  case_all_vec <- list(c(1,0,0,0), c(1,1,0,0), c(1,0,1,0), c(1,1,1,0), c(1,0,1,1), c(1,1,1,1), c(0,0,1,0), c(0,0,1,1))
  index <- as.numeric(unlist(sapply(1:8, function(x) if (all(case_all_vec[[x]] == case_vec)) return(x))))

  # create corresponding datas to case_all_vec
  data_name <- list(c("ps1 = ps_glm_list[[1]][[1]][[1]]", "ps0 = ps_glm_list[[1]][[1]][[2]]",
                      "psA1 = ps_glm_list[[1]][[1]][[3]]", "psA0 = ps_glm_list[[1]][[1]][[4]]"),

                    c("ps1 = rbind(ps_glm_list[[1]][[1]], ps_glm_list[[2]][[1]])",
                      "ps0 = rbind(ps_glm_list[[1]][[2]], ps_glm_list[[2]][[2]])",
                      "psA1 = rbind(ps_glm_list[[1]][[3]], ps_glm_list[[2]][[3]])",
                      "psA0 = rbind(ps_glm_list[[1]][[4]], ps_glm_list[[2]][[4]])"),

                    c("ps1 = rbind(ps_glm_list[[1]][[1]][[1]], ps_SL_list[[1]][[1]][[1]])",
                      "ps0 = rbind(ps_glm_list[[1]][[1]][[2]], ps_SL_list[[1]][[1]][[2]])",
                      "psA1 = rbind(ps_glm_list[[1]][[1]][[3]], ps_SL_list[[1]][[1]][[3]])",
                      "psA0 = rbind(ps_glm_list[[1]][[1]][[4]], ps_SL_list[[1]][[1]][[4]])"),

                    c("ps1 = rbind(ps_glm_list[[1]][[1]], ps_glm_list[[2]][[1]], ps_SL_list[[1]][[1]][[1]])",
                      "ps0 = rbind(ps_glm_list[[1]][[2]], ps_glm_list[[2]][[2]], ps_SL_list[[1]][[1]][[2]])",
                      "psA1 = rbind(ps_glm_list[[1]][[3]], ps_glm_list[[2]][[3]], ps_SL_list[[1]][[1]][[3]])",
                      "psA0 = rbind(ps_glm_list[[1]][[4]], ps_glm_list[[2]][[4]], ps_SL_list[[1]][[1]][[4]])"),

                    c("ps1 = rbind(ps_glm_list[[1]][[1]][[1]], ps_SL_list[[1]][[1]], ps_SL_list[[2]][[1]])",
                      "ps0 = rbind(ps_glm_list[[1]][[1]][[2]], ps_SL_list[[1]][[2]], ps_SL_list[[2]][[2]])",
                      "psA1 = rbind(ps_glm_list[[1]][[1]][[3]], ps_SL_list[[1]][[3]], ps_SL_list[[2]][[3]])",
                      "psA0 = rbind(ps_glm_list[[1]][[1]][[4]], ps_SL_list[[1]][[4]], ps_SL_list[[2]][[4]])"),

                    c("ps1 = rbind(ps_glm_list[[1]][[1]], ps_glm_list[[2]][[1]], ps_SL_list[[1]][[1]], ps_SL_list[[2]][[1]])",
                      "ps0 = rbind(ps_glm_list[[1]][[2]], ps_glm_list[[2]][[2]], ps_SL_list[[1]][[2]], ps_SL_list[[2]][[2]])",
                      "psA1 = rbind(ps_glm_list[[1]][[3]], ps_glm_list[[2]][[3]], ps_SL_list[[1]][[3]], ps_SL_list[[2]][[3]])",
                      "psA0 = rbind(ps_glm_list[[1]][[4]], ps_glm_list[[2]][[4]], ps_SL_list[[1]][[4]], ps_SL_list[[2]][[4]])"),

                    c("ps1 = ps_SL_list[[1]][[1]][[1]]", "ps0 = ps_SL_list[[1]][[1]][[2]]",
                      "psA1 = ps_SL_list[[1]][[1]][[3]]", "psA0 = ps_SL_list[[1]][[1]][[4]]"),

                    c("ps1 = rbind(ps_SL_list[[1]][[1]], ps_SL_list[[2]][[1]])",
                      "ps0 = rbind(ps_SL_list[[1]][[2]], ps_SL_list[[2]][[2]])",
                      "psA1 = rbind(ps_SL_list[[1]][[3]], ps_SL_list[[2]][[3]])",
                      "psA0 = rbind(ps_SL_list[[1]][[4]], ps_SL_list[[2]][[4]])"))

  # produce datas on index number
  ps_data <- sapply(1:4, function(x) list(eval(parse(text = data_name[[index]][[x]]))))
  xlab_name <- c("log(1/g_1)", "log(1/g_0)", "log(1/g_1)[A=1]", "log(1/g_1)[A=0]")

  # create 4 plots and combine to one plot
  library(ggplot2)
  cols <- c("red", "green", "blue", "orange")
  ps_plot <- sapply(1:4, function(y) list(ggplot(ps_data[[y]], aes(fill = ps_data[[y]]$Method, x = ps_data[[y]]$Propensity))+ geom_density(alpha=0.3)+
                                            theme(legend.title = element_blank(),legend.position="top",legend.key.size = unit(0.3, "cm"))+
                                            scale_fill_manual(values = cols, labels = c(as.character(unique(ps_data[[y]]$Method))))+xlab(xlab_name[y])))

  ps_plots <- ggpubr::ggarrange(ps_plot[[3]], ps_plot[[4]], ps_plot[[1]], ps_plot[[2]], labels = c("A","B","C","D"), ncol=2, nrow=2)
  par(ask = FALSE)
  return(ps_plots)
}




############################################################################################################
#------------------------------------------------ Stage 2 --------------------------------------------------
############################################################################################################

#---------------------- function to fit glm on Y in the observed data -------------------------
# input: Y, A, W, Yform, outcome type(continuous or binary)
# output: the coefficients of the glm and valide variable names for two treatment groups

.fit_glm_Y <- function(Y, A, W, Yform, outcome_type){
  ObsData <- data.frame(Y, A, W)

  # split the observed dataset to two treatment groups
  base1 <- ObsData[ObsData$A ==1,]
  base0 <- ObsData[ObsData$A ==0,]

  # regress the outcome separately in two groups
  if (identical(outcome_type,"continuous")){
    regY1 <- glm(Yform, data = base1, family = "gaussian")
    regY0 <- glm(Yform, data = base0, family = "gaussian")
  }
  if (identical(outcome_type,"binary")){
    regY0 <- glm(Yform, data = base0, family = "binomial")
    regY1 <- glm(Yform, data = base1, family = "binomial")
  }
  return(list(fit_y0=regY0, fit_y1=regY1))
}




#----------------------------- function to produce one bootstrap data -------------------------------

# This function is using the coeffs of the outcome regression in the original dataset to predict
#       the outcome in the bootstrap data then create the bootstrap dataset
# input: A, W, outcome type(continuous or binary), and the coeffs of outcome regression and valide
#        variable names for two treatment groups
# output: one bootstrap dataset

.sim_data <- function(A, W, Yform, outcome_type, fit_y1, fit_y0)
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
  noSpace.Yform <- substring(stringr::str_replace_all(string = Yform, pattern = " ", repl= ""), 2)
  designM1 <- model.matrix (as.formula(noSpace.Yform), BaseBData1)
  designM0 <- model.matrix (as.formula(noSpace.Yform), BaseBData0)

  # get the new outcome based on the coeffs of glm on original data
  if (identical(outcome_type,"continuous")){
    Ypred1 <- designM1[, !is.na(coef(fit_y1))] %*% fit_y1$coeff[!is.na(coef(fit_y1))]
    Ypred0 <- designM0[, !is.na(coef(fit_y0))] %*% fit_y0$coeff[!is.na(coef(fit_y0))]
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
#------------------------------------------------ Stage 3 --------------------------------------------------
############################################################################################################

#-------------------------- main function to compute the true effect of bootstrap data ---------------------

# This function is to compute the true effect of one bootstrap dataset.
# input: one bootstrap dataset, outcome type (binary or continuous) and the coeffs of outcome regression
#        and valid variable names for two treatment groups
# output: one numerical true value

.trueVal <- function (W, Yform, outcome_type, fit_y1, fit_y0){

  # create design matrix on Yform
  noSpace.Yform <- substring(stringr::str_replace_all(string = Yform, pattern = " ", repl = ""), 2)
  designM <- model.matrix (as.formula(noSpace.Yform), W)

  # compute two potential outcomes of original data
  if (identical(outcome_type,"continuous")){
    yy1 <- designM[, !is.na(coef(fit_y1))] %*% fit_y1$coeff[!is.na(coef(fit_y1))]
    yy0 <- designM[, !is.na(coef(fit_y0))] %*% fit_y0$coeff[!is.na(coef(fit_y0))]
  }
  if (identical(outcome_type,"binary")){
    yy1 <- rbinom(dim(W)[1],1, plogis( designM[, !is.na(coef(fit_y1))] %*% coef(fit_y1)[!is.na(coef(fit_y1))] ))
    yy0 <- rbinom(dim(W)[1],1, plogis( designM[, !is.na(coef(fit_y0))] %*% coef(fit_y0)[!is.na(coef(fit_y0))] ))
  }

  # compute the difference of the mean of two potential outcomes
  trueVal <- mean(yy1, na.rm = TRUE) - mean(yy0, na.rm = TRUE)
  return(trueVal)
}




############################################################################################################
#------------------------------------------------ Stage 4 --------------------------------------------------
############################################################################################################

#------------------------- helper function to estimate ate by TMLE ------------------------------------
# gform: user defined form
# SL.library: user defined library in SL
# output: the estimated ATE and confidence intervals

.tmle_mod <- function(Y, A, W, Qform, gform = NULL, SL.library, gbound = gbound, outcome_type){
  library(tmle)
  if (identical(outcome_type,"continuous")){
    tmlemod <- tmle(Y = Y, A = A, W = W, Qform = Qform, gform = gform, g.SL.library = SL.library, gbound = gbound, family = "gaussian")}

  if (identical(outcome_type,"binary")){
    tmlemod <- tmle(Y = Y, A = A, W = W, Qform = Qform, gform = gform, g.SL.library = SL.library, gbound = gbound, family = "binomial")}

  ate <- tmlemod$estimates$ATE$psi
  ci <- tmlemod$estimates$ATE$CI
  return(list(ATE = ate, CI = ci))
}




#------------------------ main function of TMLE to estimate ATE with flexible g (SL or GLM) ----------------------

# input: one observed dataset with outcome- Y(binary and/or continuous),binary treatment A and potential confounders W
# gform -- optional regression formula for estimating the  P(A=1|W)
# gbound -- bound on P(A=1|W), defaults to 0.025
# gGLM -- logical: if TRUE, use GLM; otherwise use SL
# User could define gform in GLM, default is A~W
# SL.library -- prediction algorithms for data adaptive estimation of g
# returns: the estimated ate and 95% confidence interval

TMLE_ate <- function (ObsData, Yform, outcome_type, gGLM, gbound = 0.025, gform = NULL, SL.library){

  library(tmle)
  Y <- ObsData$Y
  A <- ObsData$A
  W <- ObsData[, !names(ObsData) %in% c("Y", "A")]
  noSpace.Yform <- substring(stringr::str_replace_all(string = Yform, pattern = " ", repl = ""), 3)
  Qform <- paste("Y~A+", paste(noSpace.Yform, collapse = "+"))

  # according to the logical value of gGLM, use GLM or SL to estimate ate
  if (gGLM){ # if gGLM is true, use GLM
    if(!is.null(gform)){ # use user supplied gform
      eff_tmle <- .tmle_mod(Y, A, W, gbound = gbound, outcome_type = outcome_type, Qform = Qform, gform = gform) }
    else { # use default gform
      gform <- paste("A~", paste(colnames(W), collapse="+"))
      eff_tmle <- .tmle_mod(Y, A, W, gbound = gbound, outcome_type = outcome_type, Qform = Qform, gform = gform) }
  }
  else { # if gGLM is FALSE, use SL
    library(SuperLearner)
    eff_tmle <- .tmle_mod(Y, A, W, gbound = gbound, outcome_type = outcome_type, Qform = Qform, SL.library = SL.library)
  }
  return(eff_tmle) # return list type values: ate and CI
}




#--------------------------- helper functions to estimate ate by AIPTW --------------------------------

#----------- compute weight according to different gbound --------------

# use the bound function to compute the weight
.weight <- function (x, gbound){
  if (gbound ==0){weigh <- 1/x}
  if (gbound!=0){weigh <- 1/.bound(x, c(gbound, 1-gbound))}
  return(weigh)
}

#-------------- compute the weight via GLM ----------------
# given dataset, gform and bounds, use GLM to compute weights
.weight.glm <- function(ObsData, gform, gbound){

  ps <- glm(gform, data = ObsData, family = "binomial")
  designM <- model.matrix (as.formula(gform), ObsData)
  ps.pred <- plogis( designM[, !is.na(coef(ps))] %*% coef(ps)[!is.na(coef(ps))] )

  ps.obs <- ifelse(ObsData$A == 0, 1- ps.pred, ps.pred)
  w <- .weight(ps.obs, gbound)# use the weight function
  return(w)
}


#-------------------------- main function of AIPTW to estimate ATE with flexible g (SL or GLM) -----------------------

# input: one observed dataset with outcome- Y(binary and/or continuous),binary treatment A and potential confounders W
# gform -- optional regression formula for estimating the  P(A=1|W)
# gbound -- bound on P(A=1|W), defaults to 0.025
# gGLM -- logical: if TRUE, use GLM; otherwise use SL
# User could define gform in GLM, default is A~W
# SL.library -- prediction algorithms for data adaptive estimation of g
# returns: the estimated ate and 95% confidence interval
#
AIPTW_ate <- function (ObsData, Yform, outcome_type, gGLM, gbound = 0.025, gform = NULL, SL.library){

  Y <- ObsData$Y
  A <- ObsData$A
  W <- ObsData[, !names(ObsData) %in% c("Y", "A")]
  n <- length(Y)

  # create design matrix on simulated data
  noSpace.Yform <- substring(stringr::str_replace_all(string = Yform, pattern = " ", repl = ""), 2)
  designM <- model.matrix (as.formula(noSpace.Yform), ObsData)

  # get the potential outcome based on the coeffs of glm on original data
  if (identical(outcome_type,"continuous")){

    reg1 <- glm(Yform, data = ObsData[A ==1,], family = "gaussian")
    reg0 <- glm(Yform, data = ObsData[A ==0,], family = "gaussian")

    Y1_hat <- designM[, !is.na(coef(reg1))] %*% reg1$coeff[!is.na(coef(reg1))]
    Y0_hat <- designM[, !is.na(coef(reg0))] %*% reg0$coeff[!is.na(coef(reg0))]
  }

  if (identical(outcome_type,"binary")){

    reg1 <- glm(Yform, data = ObsData[A ==1,], family = "binomial")
    reg0 <- glm(Yform, data = ObsData[A ==0,], family = "binomial")

    Y1_hat <- rbinom(n, 1, plogis(designM[, !is.na(coef(reg1))] %*% reg1$coeff[!is.na(coef(reg1))]))
    Y0_hat <- rbinom(n, 1, plogis(designM[, !is.na(coef(reg0))] %*% reg0$coeff[!is.na(coef(reg0))]))
  }

  # estimate propensity score
  if (gGLM){
    if(!is.null(gform)){
      w <- .weight.glm(ObsData, gform, gbound) }
    else {
      gform=paste("A~", paste(colnames(W), collapse = "+"))
      w <- .weight.glm(ObsData, gform, gbound)}
  }
  else {
    library(SuperLearner)
    ps.SL <- suppressWarnings(SuperLearner(Y = A, X = W, family = binomial(), SL.library = SL.library))
    ps.obs <- ifelse(ObsData$A == 0, 1- ps.SL$SL.predict, ps.SL$SL.predict)
    w <- .weight(ps.obs, gbound)
  }

  # estimate ate
  mu1 <- sum((Y - Y1_hat)*A*w)/n + mean(Y1_hat); mu0 <- sum((Y - Y0_hat)*(1-A)*w)/n + mean(Y0_hat)
  res_ate <- mu1 - mu0
  IF <- ((Y-Y1_hat)*A*w + Y1_hat)-((Y - Y0_hat)*(1-A)*w + Y0_hat) - res_ate # compute the influence function
  ci <- res_ate + c(-1,1)*1.96*sqrt( sum(IF^2)/(n^2)) # compute the standard error then the 95% CI

  eff_aiptw <- list(ATE = res_ate, CI = ci)
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
# returns: the estimated ate and 95% confidence interval

IPTW_ate <- function (ObsData, outcome_type, gGLM, gbound = 0.025, gform = NULL, SL.library){
  Y <- ObsData$Y
  A <- ObsData$A
  W <- ObsData[, !names(ObsData) %in% c("Y", "A")]

  # estimate propensity score
  if (gGLM){
    if(!is.null(gform)){
      w <- .weight.glm(ObsData, gform, gbound) }
    else {
      gform <- paste("A~", paste(colnames(W), collapse = "+"))
      w <- .weight.glm(ObsData, gform, gbound) }
  }
  else {
    library(SuperLearner)
    ps.SL <- SuperLearner(Y = A, X = W, family = binomial(), SL.library = SL.library)
    ps.obs <- ifelse(ObsData$A == 0, 1- ps.SL$SL.predict, ps.SL$SL.predict)
    w <- .weight(ps.obs, gbound)
  }
  # estimate ate
  mod <- lm(Y~ A,data = ObsData,weights = w)
  res_ate <-summary(mod)$coef[2,1]
  # compute 95% CI
  CI <- res_ate + 1.96*c(-1,1)*sqrt(sandwich::vcovHC(mod)[2,2])
  eff_iptw <- list(ATE = res_ate, CI = CI)
  return(eff_iptw)# return list type values: ate, CI
}





############################################################################################################
#------------------------------------------------ Stage 5 --------------------------------------------------
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
  # if (class(verbose)!="logical"){stop("verbose must be TRUE or FALSE")}

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
# gform -- optional regression formula for estimating the  P(A=1|W)
# gbound -- bound on P(A=1|W), defaults to 0.025
# gGLM -- logical: if TRUE, use GLM; otherwise use SL
# User could define gform in GLM, default is A~W
# SL.library -- prediction algorithms for data adaptive estimation of g
# returns: the estimated ate and 95% confidence intervals for three methods: IPTW, AIPTW, TMLE

bdt <- function (Y, A, W, outcome_type, M, gbound = 0.025, gGLM, gform = NULL, SL.library = NULL){

  bias_TMLE <- NULL; bias_IPTW <- NULL; bias_AIPTW <- NULL
  cov_TMLE <- cov_IPTW <- cov_AIPTW <- 0
  ps.all.GLM <- NULL; ps.all.SL <- NULL

  # check all inputs
  .check_var2(Y, A, W, outcome_type, M, gbound, gGLM, gform, SL.library)

  # create the outcome regression formula
  Yform <- paste ("Y~1+", paste(colnames(W), collapse = "+"))

  # compute the coeffs derived from the glm of outcome for original dataset
  fit_y0 <- .fit_glm_Y(Y, A, W, Yform, outcome_type)$fit_y0
  fit_y1 <- .fit_glm_Y(Y, A, W, Yform, outcome_type)$fit_y1

  true_eff <- .trueVal(W, Yform = Yform, outcome_type = outcome_type, fit_y1, fit_y0)

  # loops to create N bootstrap datasets
  for (i in 1:M){
    bootdata <- .sim_data( A, W, Yform = Yform, outcome_type = outcome_type, fit_y1, fit_y0)

    if (!is.null(gform) & gGLM){
      ps.GLM <- data.frame(.ps_GLM(A = bootdata$A, W = bootdata[, !names(bootdata) %in% c("Y", "A")], gform = gform, gbound = gbound)$ps_glm)
      names(ps.GLM) <- c("A", "ps")
      ps.all.GLM <- rbind(ps.all.GLM, ps.GLM)
    }
    if (!is.null(SL.library) & !gGLM){
      ps.SL <- data.frame(.ps_SL(A= bootdata$A, W = bootdata[, !names(bootdata) %in% c("Y", "A")], SL.library = SL.library, gbound = gbound)$ps_sl)
      names(ps.SL) <- c("A", "ps")
      ps.all.SL <- rbind(ps.all.SL, ps.SL)
    }

    # compute estimates of tmle, iptw, and aiptw of the bootstrap dataset
    tmle_res <- TMLE_ate(bootdata, Yform = Yform, outcome_type = outcome_type, gGLM, gbound = gbound, gform = gform, SL.library = SL.library)
    iptw_res <- IPTW_ate(bootdata, outcome_type = outcome_type, gGLM, gbound = gbound, gform = gform, SL.library = SL.library)
    aiptw_res <- AIPTW_ate(bootdata, Yform = Yform, outcome_type = outcome_type, gGLM, gbound = gbound, gform = gform, SL.library = SL.library)

    # compute the bias
    bias_TMLE[i] <- tmle_res$ATE-true_eff
    bias_IPTW[i] <- iptw_res$ATE-true_eff
    bias_AIPTW[i] <- aiptw_res$ATE-true_eff

    # compute the times that true value is within the intervals of three methods
    cov_TMLE <- cov_TMLE + ifelse(true_eff >= tmle_res$CI[1] & true_eff <= tmle_res$CI[2], 1, 0 )
    cov_IPTW <- cov_IPTW + ifelse(true_eff >= iptw_res$CI[1] & true_eff <= iptw_res$CI[2], 1, 0 )
    cov_AIPTW <- cov_AIPTW + ifelse(true_eff >= aiptw_res$CI[1] & true_eff <= aiptw_res$CI[2], 1, 0 )
  }
  # print status message
  if (gGLM){
    if (!is.null(gform)){
      if (!is.null(SL.library)) {# if user set gGLM is TRUE and define SL function simultaneously
        cat("Estimating ate using user-supplied regression formula and ignoring supplied SuperLearner functions for g \n\n")}
      else (cat("Estimating ate using user-supplied regression formula for g \n\n"))
    }
    else (cat("Estimating ate using default regression formula for g: A ~ W \n\n"))
  }
  else {
    if (!is.null(SL.library)){
      if (!is.null(gform)){
        cat("Estimating ate using user-supplied SuperLearner functions and ignoring supplied regression formula for g \n\n")}
      else {cat("Estimating ate using user-supplied SuperLearner functions for g \n\n")}
    }
    else (stop("Please specify SuperLearner library for g. \n\n"))
  }

  results <- list(true_Effect = true_eff, gbound =  gbound, bias_IPTW = bias_IPTW, bias_AIPTW = bias_AIPTW, bias_TMLE = bias_TMLE,
                  cov_IPTW = cov_IPTW/M, cov_AIPTW = cov_AIPTW/M, cov_TMLE = cov_TMLE/M,
                  ps_GLM = ps.all.GLM, ps_SL = ps.all.SL)
  class(results) <- "bdt"
  return(results)
}




#----------------- print functions to output the estimated ate bias of AIPTW, IPTW, TMLE -------------------
# input the results of main bias_boot function
# return the bias of the three estimated ate
bias.bdt <- function(object,...){
  if(identical(class(object), "bdt")){
    cat("Eestimated bias of ATE with IPTW, AIPTW, TMLE:\n\n")

    bias.ate <- data.frame(cbind(object$bias_IPTW, object$bias_AIPTW, object$bias_TMLE))
    colnames(bias.ate) <- c("bias_IPTW", "bias_AIPTW", "bias_TMLE")
  } else {
    stop("Object must have class 'bdt'")}
  return(bias.ate)
}
bias <- function(x){ UseMethod("bias",x)}




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

    summary.bdt <- list(true.effect = object$true_Effect, gbound = object$gbound,
                        bias.iptw = summary(object$bias_IPTW), bias.aiptw = summary(object$bias_AIPTW), bias.tmle = summary(object$bias_TMLE),
                        cov.iptw = object$cov_IPTW, cov.aiptw = object$cov_AIPTW, cov.tmle = object$cov_TMLE,
                        ps1 = ps1.sum, ps0 = ps0.sum, ps11 = ps11.sum, ps00 = ps00.sum)
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

    i <- specify_decimal(x$bias.iptw, 5); a <- specify_decimal(x$bias.aiptw, 5); t <- specify_decimal(x$bias.tmle, 5)
    cat("\nSummary of ATE bias based on IPTW, AIPTW and TMLE:\n\n")
    bias_data <- data.frame(c(" IPTW:","AIPTW:"," TMLE:"), rbind(i, a, t))
    names(bias_data) <- c("      ","   Min.  ",  " 1st Qu. "," Median ", "  Mean  ", " 3rd Qu. ", "  Max.  ")
    gdata::write.fwf(bias_data)

    cat("\n\nCoverage rate of IPTW: ", x$cov.iptw, " \n")
    cat("Coverage rate of AIPTW: ", x$cov.aiptw, " \n")
    cat("Coverage rate of TMLE: ", x$cov.tmle, " \n")

    cat("\n\nSummary of propensity scores truncated at gbound", x$gbound, " \n")
    p1 <- specify_decimal(x$ps1,6); p0 <- specify_decimal(x$ps0,6)
    p11 <- specify_decimal(x$ps11,6); p00 <- specify_decimal(x$ps00,6)
    cat("\n                                  Min.       1st Qu.     Median      Mean      3rd Qu.   Max.")
    cat("\n     P(A=1|X) for all subjects: ", p1[1], " ", p1[2], " ",p1[3], " ",p1[4], " ",p1[5], " ",p1[6])
    cat("\n     P(A=0|X) for all subjects: ", p0[1], " ", p0[2], " ",p0[3], " ",p0[4], " ",p0[5], " ",p0[6])
    cat("\n    P(A=1|X) for subgroups A=1: ", p11[1], " ", p11[2], " ",p11[3], " ",p11[4], " ",p11[5], " ",p11[6])
    cat("\n    P(A=0|X) for subgroups A=0: ", p00[1], " ", p00[2], " ",p00[3], " ",p00[4], " ",p00[5], " ",p00[6])
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
plot.bdt <- function(object,...){
  if(identical(class(object), "bdt")){

    #--------------------- boxplot of bias ate ------------------
    N <- length(object$bias_TMLE)
    type_est <- c(rep("IPTW",N), rep("AIPTW",N),rep("TMLE",N) )
    est_data <- data.frame(type_est, c(object$bias_IPTW, object$bias_AIPTW, object$bias_TMLE))
    colnames(est_data) <- c("type_est", "Bias_ATE")
    title <- paste("Boxplots of the ATE bias with different estimators")
    anno.data <- data.frame(type = c("IPTW", "AIPTW", "TMLE"), high = min(est_data$Bias_ATE)-0.03,
                            coverage = round(c(object$cov_IPTW, object$cov_AIPTW, object$cov_TMLE), 2))
    library(ggplot2)
    est_plot <- ggplot(est_data, aes(x = as.factor(type_est), y = Bias_ATE))+
      geom_boxplot(outlier.size = 0.4, fill= "cornflowerblue", notch = FALSE)+
      geom_rug(color = "black") + ggtitle(title)+
      geom_point(position = "jitter", size=0.7,color="blue", alpha=0.5)+
      geom_label(data = anno.data, aes(x = type,  y = high, label = paste("Cov: ", coverage,sep="")),
                 label.padding = unit(0.35, "lines"), # Rectangle size around label
                 label.size = 0.25,colour = "red", size = 2.5,fill = "white") +
      geom_hline(yintercept = 0,col = "darkred",linetype = "dotted") +
      ylim(min(est_data$Bias_ATE)-0.03,max(est_data$Bias_ATE) + 0.03) +
      ylab("ATE bias") + theme(axis.title.x = element_blank(),axis.text.x = element_text(angle=0))


    #------------------- density plot of ps ------------------

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
  } else {
    stop("Object must have class 'bdt'")
    fig_box_density <- NULL
  }
  par(ask=TRUE)
  return(fig_box_density)
}





