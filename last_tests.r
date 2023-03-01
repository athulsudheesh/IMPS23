source("last_utils.r")
#=========================================================================
#                  Run all tests for independence
#=========================================================================
# Data in each field of out looks like below
#          AUC Acc Spe Sen
# thresh1   0   0   0   0
# thresh2   0   0   0   0
run.tests <- function(simdata, pvals,thresholds){
    x2res <- mires <- aicres <- bicres <- vector()

    for (i in 1:length(pvals)){

        nptres <- nptests(simdata, pvals[i])
        x2 <- nptres$x2
        mi <- nptres$mi
 
        x2res <- rbind(x2res, c(x2$acc, x2$spe, x2$sen, x2$auc))
        mires <- rbind(mires, c(mi$acc, mi$spe, mi$sen, mi$auc))
    }
    for (i in 1:length(thresholds)){

        ptres <- ptests(simdata, thresholds[i])
        aic <- ptres$aic
        bic <- ptres$bic

        aicres <- rbind(aicres, c(aic$acc, aic$spe, aic$sen, aic$auc))
        bicres <- rbind(bicres, c(bic$acc, bic$spe, bic$sen, bic$auc))
    }

    # Formatting the results data structure 
    res <- lapply(list(x2res= x2res,
                      mires = mires,
                      aicres = aicres,
                      bicres = bicres),function(x) as.data.frame(x))

    colnames(res$x2res) <- 
        colnames(res$mires) <- 
        colnames(res$aicres) <- 
        colnames(res$bicres) <- c("Acc", "Spe","Sen","AUC")
     
    out = list(
        x2 = res$x2res,
        mi = res$mires,
        aic = res$aicres,
        bic = res$bicres
    )
    return(out)
}


#=========================================================================
#                  Run Chi-Sq and MI test for the given pval
#=========================================================================
nptests <- function(simdata,pval){

    dep_res <- npdep(simdata,pval)
    ind_res <- npind2(simdata,pval)
  
    x2pred <- as.data.frame(rbind(dep_res$x2,ind_res$x2))
    mipred <- as.data.frame(rbind(dep_res$mi,ind_res$mi))

    x2pred <- as.data.frame(lapply(x2pred, as.factor))
    mipred <- as.data.frame(lapply(mipred, as.factor))
    

    x2 = getReport(x2pred)
    mi = getReport(mipred)
    out = list(x2 = x2, mi = mi)
    return(out) 

}

#=========================================================================
#                  Run AIC and BIC test for given threshold
#=========================================================================
ptests <- function(simdata,threshold){
    
    dep_res <- pdep(simdata,threshold)
    ind_res <- pind2(simdata,threshold)

    aicpred <- as.data.frame(rbind(dep_res$aic,ind_res$aic))
    bicpred <- as.data.frame(rbind(dep_res$bic,ind_res$bic))

    aicpred <- as.data.frame(lapply(aicpred, as.factor))
    bicpred <- as.data.frame(lapply(bicpred, as.factor))
    

    aic = getReport(aicpred)
    bic = getReport(bicpred)
    out = list(aic = aic, bic = bic)
    return(out) 

}
#=========================== DEPENDENCE TESTING ==========================
#=========================================================================
#                  Run Chi-Sq and MI test to check for dependence
#=========================================================================
npdep <- function(simdata, pval){
     x2 <- vector()
     mi <- vector()
     
     QMat <- simdata$QMat
     XMat <- simdata$X
     items <- simdata$itemnames
     itemcomb <- t(combn(items,2))

    for (i in 1:ncol(QMat)){
        corr_items <-items[as.logical(QMat[,i])]
        test_set <- as.data.frame(t(combn(corr_items,2)))

        for (j in 1:nrow(test_set)){
            x2_pred <- cd.x2(test_set$V1[j], test_set$V2[j],XMat,pval)
            mi_pred <- cd.mi(test_set$V1[j], test_set$V2[j],XMat,pval)
            
            x2 <- rbind(x2,c("1", x2_pred))
            mi <- rbind(mi, c("1", mi_pred))
        }
    }
     

     out = list(x2=x2,mi=mi)
     return(out)
}

#=========================================================================
#                  Run AIC and BIC test to check for dependence
#=========================================================================
pdep <- function(simdata,threshold) {
     aic <- vector()
     bic <- vector()
     
     QMat <- simdata$QMat
     XMat <- simdata$X
     items <- simdata$itemnames
     itemcomb <- t(combn(items,2))

    for (i in 1:ncol(QMat)){
        corr_items <-items[as.logical(QMat[,i])]
        test_set <- as.data.frame(t(combn(corr_items,2)))

        for (j in 1:nrow(test_set)){
            nested_models <- competing_dep_models(test_set$V1[i], test_set$V2[i], XMat)
                
            aic_pred <- ci.aic(nested_models$full, nested_models$reduced, threshold)
            bic_pred <- ci.bic(nested_models$full, nested_models$reduced, threshold)    
            
            aic <- rbind(aic, c("1", aic_pred))
            bic <- rbind(bic,c("1", bic_pred))
        }
    }
     

     out = list(aic=aic,bic=bic)
     return(out)
}

#=========================== INDEPENDENCE TESTING ========================
#=========================================================================
#                  Run Chi-Sq and MI test to check for independence
#=========================================================================
npind2 <- function(simdata, pval){
    x2 <- vector()
    mi <- vector()
    
    QMat <- simdata$QMat
    XMat <- simdata$X
    itemnames <- simdata$itemnames
    skills <- simdata$skills
    rownames(QMat) <- itemnames
    itemcomb <- t(combn(itemnames,2))
    for (i in 1:nrow(itemcomb)){
        indc <- as.logical(QMat[itemcomb[i,][1],] * QMat[itemcomb[i,][2],])
        
        if (sum(indc) == 0)
        {
            x2_pred <- cd.x2(itemcomb[i,][1], itemcomb[i,][2],XMat,pval)
            mi_pred <- cd.mi(itemcomb[i,][1], itemcomb[i,][2],XMat,pval)
      

        #x2_pred <- ci.x2(itemcomb[i,][1], itemcomb[i,][2],skills[indc], dataset = data, pval)
        #mi_pred <- ci.mi(itemcomb[i,][1], itemcomb[i,][2],skills[indc], dataset = data, pval)

        x2 <- rbind(x2, c("0", x2_pred))
        mi <- rbind(mi,c("0", mi_pred))
        }
    }
    out = list(x2=x2,mi=mi)
} 

npind <- function(simdata, pval){
    x2 <- vector()
    mi <- vector()
    
    QMat <- simdata$QMat
    data <- simdata$XplusA
    itemnames <- simdata$itemnames
    skills <- simdata$skills
    rownames(QMat) <- itemnames
    itemcomb <- t(combn(itemnames,2))
    for (i in 1:nrow(itemcomb)){
        indc <- as.logical(QMat[itemcomb[i,][1],] * QMat[itemcomb[i,][2],])
        
        x2_pred <- ci.x2(itemcomb[i,][1], itemcomb[i,][2],skills[indc], dataset = data, pval)
        mi_pred <- ci.mi(itemcomb[i,][1], itemcomb[i,][2],skills[indc], dataset = data, pval)

        x2 <- rbind(x2, c("0", x2_pred))
        mi <- rbind(mi,c("0", mi_pred))
    }
    out = list(x2=x2,mi=mi)
    return(out)
}

#=========================================================================
#                  Run AIC and BIC test to check for independence
#=========================================================================
pind <- function(simdata, threshold){
    aic <- vector()
    bic <- vector()
    
    QMat <- simdata$QMat
    data <- simdata$XplusA
    itemnames <- simdata$itemnames
    rownames(QMat) <- itemnames
    skills <- simdata$skills
    itemcomb <- t(combn(itemnames,2))

    for (i in 1:nrow(itemcomb)){
        indc <- as.logical(QMat[itemcomb[i,][1],] * QMat[itemcomb[i,][2],])
        
        nested_models <- competing_indep_models(itemcomb[i,][1], itemcomb[i,][2],skills[indc], dataset = data)

        aic_pred <- ci.aic(nested_models$full, nested_models$reduced, threshold)
        bic_pred <- ci.bic(nested_models$full, nested_models$reduced, threshold)

        aic <- rbind(aic,c("0", aic_pred))
        bic <- rbind(bic,c("0", bic_pred))
    }
    
    out = list(aic=aic,bic=bic)
    return(out)
}

pind2 <- function(simdata, threshold){
    aic <- vector()
    bic <- vector()
    
    QMat <- simdata$QMat
    XMat <- simdata$X
    itemnames <- simdata$itemnames
    rownames(QMat) <- itemnames
    skills <- simdata$skills
    itemcomb <- t(combn(itemnames,2))

    for (i in 1:nrow(itemcomb)){
        indc <- as.logical(QMat[itemcomb[i,][1],] * QMat[itemcomb[i,][2],])
        
        nested_models <- competing_dep_models(itemcomb[i,][1], itemcomb[i,][2], XMat)
        #nested_models <- competing_indep_models(itemcomb[i,][1], itemcomb[i,][2],skills[indc], dataset = data)
        if(sum(indc)==0){
        
        aic_pred <- ci.aic(nested_models$full, nested_models$reduced, threshold)
        bic_pred <- ci.bic(nested_models$full, nested_models$reduced, threshold)
        
        aic <- rbind(aic,c("0", aic_pred))
        bic <- rbind(bic,c("0", bic_pred))
        }
    }
    
    out = list(aic=aic,bic=bic)
    return(out)
}



cd.x2 <- function(var1, var2, dataset, pval){
    test_result <- ci.test(var1, var2, data =dataset, test = "x2")
            if (test_result$p.val <= pval){
            # The null hypothesis that the given set of nodes are independent
            # When we have p <0.05, we reject the null
                dep <- "1"
            }
            else {
                dep <- "0"
            }
    return(dep)
}
# #Mutual Information based test for dependency 
cd.mi <- function(var1, var2, dataset, pval){
    test_result  <- ci.test(var1, var2, data =dataset, test = "mi")
             if (test_result$p.val <= pval){
            # The null hypothesis that the given set of nodes are independent
            # When we have p <0.05, we reject the null
                dep <- "1"
            }
            else {
                dep <- "0"
            }
    return(dep)
}

cd.aic <- function(fullmodel, reducedmodel, threshold){
    
    delta <- AIC(fullmodel) - AIC(reducedmodel)
    decision <- exp(0.5*delta)
    if (decision < threshold){
        dep <- "1"
    }

    else{
        dep <- "0"
    }
}

cd.bic <- function(fullmodel, reducedmodel, threshold){
    delta <- BIC(fullmodel) - BIC(reducedmodel)
    decision <- exp(0.5*delta)
    if (decision < threshold){
        dep <- "1"
    }

    else{
        dep <- "0"
    }
}


competing_dep_models <- function(var1, var2, dataset){
    fullmodel <- glm(get(var1) ~ 1 + get(var2), data = dataset, family = "binomial")
    reducedmodel <- glm(get(var1) ~ 1, data = dataset, family = "binomial")
    models <- list(
        full = fullmodel,
        reduced = reducedmodel
    )
}

ci.x2 <- function(var1, var2, cond_set,dataset, pval){
    test_result <- ci.test(var1, var2, cond_set, data =dataset, test = "x2")
            if (test_result$p.val <= pval){
            # The null hypothesis that the given set of nodes are independent
            # When we have p <0.05, we reject the null
                dep <- "1"
            }
            else {
                dep <- "0"
            }
    return(dep)
}

## Mutual Information test for independence 
ci.mi <- function(var1, var2, cond_set,dataset, pval){
    test_result <- ci.test(var1, var2, cond_set, data =dataset, test = "mi")
            if (test_result$p.val <= pval){
            # The null hypothesis that the given set of nodes are independent
            # When we have p <0.05, we reject the null
                dep <- "1"
            }
            else {
                dep <- "0"
            }
    return(dep)
}

ci.aic <- function(fullmodel, reducedmodel, threshold){
    
    delta <- AIC(fullmodel) - AIC(reducedmodel)
    decision <- exp(0.5*delta)
    if (decision < threshold){
        dep <- "1"
    }

    else{
        dep <- "0"
    }
}

ci.bic <- function(fullmodel, reducedmodel, threshold){
    
    delta <- BIC(fullmodel) - BIC(reducedmodel)
    decision <- exp(0.5*delta)
    if (decision < threshold){
        dep <- "1"
    }

    else{
        dep <- "0"
    }
}

competing_indep_models <- function(var1, var2, cond_set, dataset){
        full_features <- append(append(var2, cond_set),1)
        target <- var1

        rhs<- paste(full_features,collapse = "+")
        fullmodel_formula <- as.formula(paste(target, rhs, sep ="~"))
        fullmodel <- glm(fullmodel_formula,data = dataset, family = "binomial")

        reduced_features <- append(cond_set,1)
        rhs <- paste(reduced_features, collapse = "+")
        reduced_formula <- as.formula(paste(target, rhs, sep ="~"))
        reducedmodel <- glm(reduced_formula,data = dataset, family = "binomial")

        out = list(full = fullmodel,
                   reduced = reducedmodel)
        return(out)
}