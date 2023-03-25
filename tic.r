n <- c(102,  99,   108,	 76,   81,   90)
p <- c(0.539,0.525,0.528,0.724,0.617,0.555)
x <- round(n*p)
## x <- n*p
y <- cbind(x,n-x)
f <- rep(c(40,150,350),2)
(g <- gl(2,3))
m1 <- glm(y ~ g + f, family=binomial(link="logit"))
n <- nrow(y)
library(sandwich) # for OPG and Robust covariance matrices 
library(pracma) # for pseudo-inverse function 
library(psych) # for trace funcion 



GAIC <- function(model){
    opg_vcov <- vcovOPG(model)
    robust_vcov <- vcovHC(model,type="HC0")
    n <- nobs(model)
    gaic = -logLik(model)[1]/n + tr(pinv(opg_vcov) %*% robust_vcov)/n
    return(gaic)
}

GAIC(m1)

 opg_vcov <- vcovOPG(m1)
    robust_vcov <- vcovHC(m1,type="HC0")

opg_eig <- eigen(opg_vcov)$values
robust_eig <- eigen(robust_vcov)$values

min(opg_eig)
min(opg_eig) / max(opg_eig)

min(robust_eig)
min(robust_eig) / max(robust_eig)


XBIC <- function(model){
    opg_vcov <- vcovOPG(model)
    robust_vcov <- vcovHC(model,type="HC0")
    
    npar <- length(coef(model))
    n <- nobs(model)

    A <- pinv(robust_vcov) * n
    B = pinv(opg_vcov) * n 

    XBIC <- BIC(model) + tr(A %*% pinv(B)) - npar*log(2*pi) + log(det(pinv(A*n)))
   return(XBIC)
}

XBIC(m1)
