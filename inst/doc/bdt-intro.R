
library(bdt)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

#---------------------------

set.seed(1250)
n = 1000
sigma = matrix(c(2, 1, 1, 1), ncol = 2)
W = matrix(rnorm(n*2), ncol = nrow(sigma)) %*% chol(sigma)
W = W + matrix(rep(c(.5, 1),each = n), byrow = FALSE, ncol = 2)
bound = function(x, bounds){
  x[x < min(bounds)] = min(bounds); x[x > max(bounds)] = max(bounds)
  return(x)}
I1 = bound(rnorm(n,mean = 1, sd = 2), c(-3,3))
I2 = bound(rnorm(n,mean = 1, sd = 1.9), c(-3,3))
P = bound(rnorm(n,mean = 1, sd  =1.5), c(-3,3))
X = data.frame(W, I1, I2, P)
colnames(X) = c("W1", "W2", "I1", "I2",  "P")
X$W1 <- bound(X$W1,c(-3,4)); X$W2 <- bound(X$W2,c(-3,4))
A = rbinom(n, 1, plogis(0.2 + X[,"W1"] + 0.3*X[,"I1"] + X[,"W1"]*X[,"I1"]
                        - 0.2*(X[,"W2"] + X[,"I2"])^2 ))
Y = 1 + A + X[,"W1"] + 2*X[,"W2"] + 0.5*(X[,"W1"] + X[,"P"])^2 + rnorm(n)

#--------------------------

ps = ps(A, X, gform1 = "A~ W1+W2+I1+I2", gform2 = "A~ W1*I1+W2*I2",
        SL.library1 = NULL,
        SL.library2 = c("SL.glm", "SL.glmnet", "SL.gam", "SL.glm.interaction"),
        verbose = TRUE)
summary(ps)
plot(ps)
ps$fit_summaries

#--------------------------

qform = "Y~A+W1+W2+P"
gform = "A~W1+W2+I1+I2"
bdt_glm = bdt(Y, A, X, M = 10, outcome_type = "continuous",
              Qform = qform, gbound = 0, gGLM = TRUE, gform = gform)
summary(bdt_glm)
plot(bdt_glm)

lib = c("SL.gam", "SL.glm.interaction")
bdt_sl = bdt(Y, A, X, M = 10, outcome_type = "continuous", Qform = qform,
             gbound = 0, gGLM = FALSE, SL.library = lib)
summary(bdt_sl)
plot(bdt_sl)

