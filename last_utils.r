library(pROC)
library(sandwich) # for OPG and Robust covariance matrices 
library(pracma) # for pseudo-inverse function 
library(psych) # for trace funcion 
library(matrixcalc) # for checking singularity of matrices 
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

    if(is.singular.matrix(bread(model))){
        print("Bread is a singular matrix")
    }
    if(is.singular.matrix(meat(model))){
        print("Meat is a singular matrix")
    }
    if(is.singular.matrix(robust_vcov)){
        print("RobustVCOV is a singular matrix")
    }
   # if (min(opg_eig) <= 0){
   ##     print("Eigenvalues of OPGvcov less than or equal to 0")
    #}
    #if (min(robust_eig) <= 0){
    #    print("Eigenvalues of Hessianvcov less than or equal to 0")
    #}
    #if (min(opg_eig) / max(opg_eig) < 0.00001){
    #    print("OPG Condition No very small")
    #}
    #if (min(robust_eig) / max(robust_eig) < 0.00001){
    #    print("Hessian Condition No very small")
    #}

    n <- nobs(model)
    gaic = -logLik(model)[1]/n + tr(bread(model)*meat(model))/n
    return(n*gaic)
}

XBIC <- function(model){
    npar <- length(model$coefficients)
    xbic <- BIC(model) + tr(bread(model)*meat(model)) - npar*log(2*pi) - log(det(meat(model)))
}

DIMT <- function(model){
    npar <- length(model$coefficients)
    dimt <- log(det(bread(model)))/npar + log(det(meat(model)))/npar
}


check.non_singular <- function(A){
    eigvals <- eigen(A)$values

    if((min(eigvals) > 0 & max(eigvals) > 0.0001 & min(eigvals)/max(eigvals) < 0.001)){
        return(TRUE)
    }
    else {
       return(FALSE)
    }
}