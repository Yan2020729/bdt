

############################################################################################################
#------------------------------------------------ Stage 1 --------------------------------------------------
############################################################################################################

# bound x with min and max values of bounds
.bound <- function(x, bounds){
  x[x<min(bounds)] <- min(bounds); x[x>max(bounds)] <- max(bounds)
  return(x)
  # y <- x; q <- quantile(x, bounds)
  # y[x < q[1]] <- q[1]; y[x > q[2]] <- q[2]
  # return(y)
}


.convertCategorPS <- function(W, gform, remove_first_dummy = FALSE, remove_most_frequent_dummy = FALSE){
  if (any(lapply(W, is.factor) == TRUE)){
    cat_names <- names(W)[lapply(W, is.factor) == TRUE]
    # extract baseline variable names from gform
    noSpace.gform <- stringr::str_replace_all(string=gform, pattern=" ", repl="")
    g_val <- chartr(old = "+", new = " ", substring(noSpace.gform, 3))# get a string of variables in gform
    g_val <- unlist(strsplit( g_val, split = " "))# split the string

    W = fastDummies::dummy_cols(W, remove_first_dummy = remove_first_dummy,
                                remove_most_frequent_dummy=remove_most_frequent_dummy)
    W = W[, !names(W) %in% c(cat_names)]

    gform_cat <- stringr::str_c(c(unlist(sapply(1:length(g_val), function(x)
      names(W)[stringr::str_detect(names(W), g_val[x])]))), collapse = "+")
    gform <- paste0("A~", gform_cat)
  }
  return(list(W=W, gform=gform))
}


#--------------------------- helper function to produce ps via GLM -----------------------------
# input: exposure-- A(binary); set of confounders--W
#        gform-- optional regression formula (default: A~W)
# return ps.pred-- p(A=1|W) for all subjects
#        ps.A1.pred-- p(A=1|W) for subgroups with A=1
#        ps.A0.pred-- p(A=1|W) for subgroups with A=0

.ps_GLM <- function(A, W, gform, gbound, remove_first_dummy, remove_most_frequent_dummy){

  coverCate <- .convertCategorPS(W, gform, remove_first_dummy, remove_most_frequent_dummy)
  W = coverCate$W; gform = coverCate$gform

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
  # library(SuperLearner)
  ps <- SuperLearner(Y = A, X = W, family = binomial(), SL.library = SL.library)
  ps.pred <- .bound(ps$SL.predict, c(gbound, 1-gbound))
  ps_sl <- data.frame(A, ps.pred)
  return(list(ps_sl = ps_sl, fit_SL = ps))
}




#-------------- helper function to check all input variables --------------------
.check_var1 <- function(A, W, gform1, gform2, SL.library1, SL.library2, gbound, verbose){
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

  # check verbose
  if (class(verbose)!="logical"){stop("verbose must be TRUE or FALSE")}
}



#--------------------------- main function of ps for all subjects -----------------------------

# input: Y, A, W, gGLM, gorm(default is NULL), SL.library(no default value), verbose (default is FALSE)
# return: summary of p(A=1|W) for all subjects with GLM and SL and summary of fit coefficients
# if verbose is TRUE, it will show the fit summaries.

ps <- function(A, W, gform1 = NULL, gform2 = NULL, SL.library1 = NULL, SL.library2 = NULL, gbound = 0,
               verbose = TRUE, remove_first_dummy = FALSE, remove_most_frequent_dummy = FALSE){
  .check_var1(A, W, gform1, gform2, SL.library1, SL.library2, gbound, verbose)

  # always have gform1 specified if gform2 is specified
  if (is.null(gform1) & !is.null(gform2)) {gform1 = gform2; gform2 = NULL}

  # always have SL.library1 specified if SL.library2 is specified
  if (is.null(SL.library1) & !is.null(SL.library2)) {SL.library1 = SL.library2; SL.library2 = NULL}

  # create a case vector
  case_vec <- sapply(1:4, function(x)as.numeric(!is.null(list(gform1, gform2, SL.library1, SL.library2)[[x]])))

  # compute ps on case_vec
  ps_glm_list <- sapply(1:2, function(x)if (case_vec[x]==1) {
    ps.GLM <- .ps_GLM(A = A, W = W, gform = eval(as.name(paste0("gform",x))), gbound,
                      remove_first_dummy = remove_first_dummy,
                      remove_most_frequent_dummy = remove_most_frequent_dummy)$ps_glm
    n.GLM <- length(ps.GLM$ps.pred)
    ps1.GLM <- data.frame(as.factor(c(rep(paste0("GLM",x), n.GLM))),c(ps.GLM$ps.pred))
    ps0.GLM <- data.frame(as.factor(c(rep(paste0("GLM",x), n.GLM))),c((1-ps.GLM$ps.pred)))
    ps.A1.GLM <- data.frame(as.factor(c(rep(paste0("GLM",x), length(ps.GLM[A ==1,]$ps.pred)))),c(ps.GLM[A ==1,]$ps.pred))
    ps.A0.GLM <- data.frame(as.factor(c(rep(paste0("GLM",x), length(ps.GLM[A ==0,]$ps.pred)))),c((1-ps.GLM[A ==0,]$ps.pred)))
    colnames(ps1.GLM) <- colnames(ps0.GLM) <- colnames(ps.A1.GLM) <- colnames(ps.A0.GLM) <- c("Method", "Propensity")
    list(list(ps1.GLM, ps0.GLM, ps.A1.GLM, ps.A0.GLM))})
  if (is.null(ps_glm_list[[2]]) & !is.null(ps_glm_list[[1]])) {
    for (i in 1:4)  {ps_glm_list[[1]][[1]][[i]]$Method = gsub('.{1}$', '', ps_glm_list[[1]][[1]][[i]]$Method) }}

  ps_SL_list <- sapply(1:2, function(x)if (case_vec[x+2] ==1) {
    ps.SL <- .ps_SL(A=A, W=W, SL.library = eval(as.name(paste0("SL.library",x))), gbound)$ps_sl
    n.SL <- length(ps.SL$ps.pred)
    ps1.SL <- data.frame(as.factor(c(rep(paste0("SL",x), n.SL))),c(ps.SL$ps.pred))
    ps0.SL <- data.frame(as.factor(c(rep(paste0("SL",x), n.SL))),c((1-ps.SL$ps.pred)))
    ps.A1.SL <- data.frame(as.factor(c(rep(paste0("SL",x), length(ps.SL[A ==1,]$ps.pred)))),c(ps.SL[A ==1,]$ps.pred))
    ps.A0.SL <- data.frame(as.factor(c(rep(paste0("SL",x), length(ps.SL[A ==0,]$ps.pred)))),c((1-ps.SL[A ==0,]$ps.pred)))
    colnames(ps1.SL) <- colnames(ps0.SL) <- colnames(ps.A1.SL) <- colnames(ps.A0.SL) <- c("Method", "Propensity")
    list(list(ps1.SL, ps0.SL, ps.A1.SL, ps.A0.SL))})
  if (is.null(ps_SL_list[[2]])  & !is.null(ps_SL_list[[1]])) {
    for (i in 1:4)  {ps_SL_list[[1]][[1]][[i]]$Method = gsub('.{1}$', '', ps_SL_list[[1]][[1]][[i]]$Method) }}

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

  if (verbose == TRUE) { # if verbose is TRUE, return the fit summaries
    try(fit_glm_list <- sapply(1:2, function(x) if (case_vec[x]==1) {
      fit_GLM <- .ps_GLM(A = A, W = W, gform = eval(as.name(paste0("gform",x))), gbound,
                         remove_first_dummy = remove_first_dummy,
                         remove_most_frequent_dummy = remove_most_frequent_dummy)$fit_glm
      list(fit_GLM)}))
    try(fit_SL_list <- sapply(1:2, function(x) if (case_vec[x+2] ==1) {
      fit_SL <- .ps_SL(A = A, W = W, SL.library = eval(as.name(paste0("SL.library",x))), gbound)$fit_SL
      list(fit_SL)}))

    summ_name <- c("fit_glm_list[[1]][[1]]",
                   "list(fit_glm_list[[1]], fit_glm_list[[2]])",
                   "list(fit_glm_list[[1]][[1]], fit_SL_list[[1]][[1]])",
                   "list(fit_glm_list[[1]], fit_glm_list[[2]], fit_SL_list[[1]][[1]])",
                   "list(fit_glm_list[[1]][[1]], fit_SL_list[[1]], fit_SL_list[[2]])",
                   "list(fit_glm_list[[1]], fit_glm_list[[2]], fit_SL_list[[1]], fit_SL_list[[2]])",
                   "fit_SL_list[[1]][[1]]",
                   "list(fit_SL_list[[1]], fit_SL_list[[2]])")
    fit_all <- eval(parse(text = summ_name[index]))
    ps_all <- list(probabilities = ps_data, fit_summaries = fit_all)
    class(ps_all) = "ps"
    return(ps_all)
  }
  else {
    ps_all <- list(probabilities = ps_data)
    class(ps_all) = "ps"
    return(ps_all)}
}


#------------------------summary of ps ------------------------

summary.ps <- function(object,...){
  if(identical(class(object), "ps")){
    rowName <- c("P(A=1|W) for all subjects by ", "P(A=0|W) for all subjects by ", "P(A=1|W) in subgroups A=1 by ", "P(A=0|W) in subgroups A=0 by ")
    colName <- c("Min.", "1st Qu.","Median", "Mean", "3rd Qu.", "Max.")
    methods <- unique(object$probabilities[[1]]$Method)
    ALL_methods = summ_list = GLM1 = GLM2 = SL1 = SL2 = GLM = SL = NULL
    for (i in 1:length(methods)){
      summ <- t(data.frame(sapply(1:4, function(x)summary(object$probabilities[[x]][which(object$probabilities[[x]]$Method == methods[[i]]), ]$Propensity))))
      dimnames(summ) <- list(c(paste0(rowName,as.character(methods[[i]]))),colName)
      assign(as.character(methods[[i]]), summ)
      summ_list <- rbind(summ_list, summ)
    }
    if (!is.null(GLM1)&!is.null(GLM2)&!is.null(SL1)&!is.null(SL2)) {
      summary.ps = list(ALL_methods = summ_list, GLM1 = GLM1, GLM2 = GLM2, SL1 = SL1, SL2 = SL2)
    } else if (!is.null(GLM1)&!is.null(GLM2)&!is.null(SL)) {
      summary.ps = list(ALL_methods = summ_list, GLM1 = GLM1, GLM2 = GLM2, SL= SL, SL1 = SL, SL2 = SL)
    } else if (!is.null(GLM1)&!is.null(GLM2)&is.null(SL)){
      summary.ps = list(ALL_methods = summ_list, GLM1 = GLM1, GLM2 = GLM2)
    } else if (!is.null(GLM)&!is.null(SL1)&!is.null(SL2)) {
      summary.ps = list(ALL_methods = summ_list, GLM = GLM, GLM1 = GLM, GLM2 = GLM, SL1= SL1, SL2 = SL2)
    } else if(is.null(GLM) &!is.null(SL1)&!is.null(SL2)) {
      summary.ps = list(ALL_methods = summ_list, SL1= SL1, SL2 = SL2)
    } else if (!is.null(GLM) & !is.null(SL)){
      summary.ps = list(ALL_methods = summ_list, GLM = GLM, GLM1 = GLM, GLM2 = GLM, SL = SL, SL1= SL, SL2 = SL)
    } else if (!is.null(GLM) & is.null(SL)){
      summary.ps = list(ALL_methods = summ_list, GLM = GLM, GLM1 = GLM, GLM2 = GLM)
    } else {summary.ps = list(ALL_methods = summ_list, SL = SL, SL1= SL, SL2 = SL)}
    class(summary.ps) <- "summary.ps"
  } else {
    stop("Object must have class 'ps'")
    summary.ps <- NULL}
    return(summary.ps)
}


print.summary.ps <- function(x,...){
  if(identical(class(x), "summary.ps")){
    cat("\nSummary of estimated probabilities for given method(s):  \n\n")
    print(x$ALL_methods)
  }
}



############################### main function: density plots of log of weights ###############################
# input: A, W, gorm(default is NULL), SL.library(default is NULL)
# return: 2 x 2 density plots (A, B, C, D) of the log of the estimated weights for
#         A--treatment A=1 in subset of subjects with A=1
#         B--treatment A=0 in subset of subjects with A=0
#         C--treatment A=1 for all subjects
#         D--treatment A=0 for all subjects
plot.ps <- function(x, xlab_name = c("log(1/g_1)", "log(1/g_0)", "log(1/g_1)[A=1]", "log(1/g_0)[A=0]"),
                    cols = c("red", "green", "blue", "orange"), ...){
  if(identical(class(x), "ps")){
    # produce datas on index number
    ps_data <- x$probabilities

    # create 4 plots and combine to one plot
    # library(ggplot2)
    ps_plot <- sapply(1:4, function(y) list(ggplot(ps_data[[y]], aes(fill = ps_data[[y]]$Method, x = log(1/ps_data[[y]]$Propensity)))+ geom_density(alpha=0.3)+
                                              theme(legend.title = element_blank(),legend.position="top",legend.key.size = unit(0.3, "cm"))+
                                              scale_fill_manual(values = cols, labels = c(as.character(unique(ps_data[[y]]$Method))))+xlab(xlab_name[y])))

    ps_plots <- ggpubr::ggarrange(ps_plot[[3]], ps_plot[[4]], ps_plot[[1]], ps_plot[[2]], labels = c("A","B","C","D"), ncol=2, nrow=2)
    par(ask = FALSE)
  }else {
    stop("Object must have class 'ps'")
    ps_plots <- NULL}
  return(ps_plots)
}
# plot(dd)



