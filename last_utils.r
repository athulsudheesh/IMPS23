library(pROC)
library(sandwich) # for OPG and Robust covariance matrices 
library(pracma) # for pseudo-inverse function 
library(psych) # for trace funcion 
getReport <- function(data){
    colnames(data) <- c("true","pred")
    trues<- as.data.frame(factor(data$true, levels = c('1','0')))
    preds<- as.data.frame(factor(data$pred, levels = c('1','0')))
    cfd <- data.frame(trues,preds)
    colnames(cfd) <- c("trues","preds")

    cf <- confusionMatrix(reference=cfd$trues, data=cfd$preds)

    acc <- cf$overall[1]
    spe <- cf$byClass[2]
    sen<- cf$byClass[1]
    roc <- roc(as.numeric(cfd$trues), as.numeric(cfd$preds))
    auc <- auc(roc)
    out = list(acc=acc,
            spe = spe,
            sen = sen,
            auc = auc)
    return(out)
}

GAIC <- function(model){
    opg_vcov <- vcovOPG(model)
    robust_vcov <- vcovHC(model,type="HC0")

    opg_eig <- eigen(opg_vcov)$values
    robust_eig <- eigen(robust_vcov)$values
    
    if (min(opg_eig) <= 0){
        print("Eigenvalues of OPGvcov less than or equal to 0")
    }
    if (min(robust_eig) <= 0){
        print("Eigenvalues of Hessianvcov less than or equal to 0")
    }
    if (min(opg_eig) / max(opg_eig) < 0.00001){
        print("OPG Condition No very small")
    }
    if (min(robust_eig) / max(robust_eig) < 0.00001){
        print("Hessian Condition No very small")
    }

    n <- nobs(model)
    gaic = -logLik(model)[1]/n + tr(robust_vcov %*% pinv(opg_vcov))/n
    return(gaic)
}

XBIC <- function(model){
    opg_vcov <- vcovOPG(model)
    robust_vcov <- vcovHC(model,type="HC0")
    
    npar <- length(coef(model))
    n <- nobs(model)

    A <- pinv(robust_vcov) * n
    B = pinv(opg_vcov) * n 

    XBIC <- BIC(model) + tr(A %*% pinv(B)) - npar*log(2*pi) + log(det(pinv(A*n)))

}