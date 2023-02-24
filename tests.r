# Unconditional Tests for Dependency 
# #Chi-Sq test for dependency 

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

check.dependencies <- function(dat, pval,threshold){
    
        x2predicted_dep <- vector()
        mipredicted_dep <- vector()
        aicpredicted_dep <- vector()
        bicpredicted_dep <- vector()

        QMat <- dat$QMat
        XMat <- as.data.frame(dat$XMatrix)
        items <- dat$itemnames
        itemcomb <- t(combn(items,2))

        for (i in 1:ncol(QMat)){
            corr_items <-items[as.logical(QMat[,i])]
            test_set <- as.data.frame(t(combn(corr_items,2)))

            for (i in 1:nrow(test_set)){
        
                # Testing for dependencies 
                x2_pred <- cd.x2(test_set$V1[i], test_set$V2[i],XMat,pval)
                mi_pred <- cd.mi(test_set$V1[i], test_set$V2[i],XMat,pval)

                nested_models <- competing_dep_models(test_set$V1[i], test_set$V2[i], XMat)
                
                aic_pred <- cd.aic(nested_models$full, nested_models$reduced, threshold)
                bic_pred <- cd.bic(nested_models$full, nested_models$reduced, threshold)
                
                
                x2predicted_dep <- append(x2predicted_dep, x2_pred)
                mipredicted_dep <- append(mipredicted_dep, mi_pred)
                aicpredicted_dep <- append(aicpredicted_dep, aic_pred)
                bicpredicted_dep <- append(bicpredicted_dep, bic_pred)
            } # end of nrow(test_set)
        }
    truedep <- rep("1", length(x2predicted_dep))
    
    
    out <- list(
        truedep = truedep, 
        x2pred = x2predicted_dep,
        mipred = mipredicted_dep,
        aic_pred = aicpredicted_dep,
        bic_pred = bicpredicted_dep)
    return(out)
}

competing_dep_models <- function(var1, var2, dataset){
    fullmodel <- glm(get(var1) ~ 1 + get(var2), data = dataset, family = "binomial")
    reducedmodel <- glm(get(var1) ~ 1, data = dataset, family = "binomial")
    models <- list(
        full = fullmodel,
        reduced = reducedmodel
    )
}

# Independence testing 
## Chi-Sq test for independence 
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

check.independencies <- function(dataset, pval,threshold){

    x2predicted_dep <- vector()
    mipredicted_dep <- vector()
    aicpredicted_dep <- vector()
    bicpredicted_dep <- vector()
    
    QMat <- dataset$QMat
    data <- dataset$XplusA
    itemnames <- colnames(data)[1:nrow(QMat)]
    rownames(QMat) <- itemnames
    skills <- colnames(QMat)
    itemcomb <- t(combn(itemnames,2))
    
    for (i in 1:nrow(itemcomb)){
        indc <- as.logical(QMat[itemcomb[i,][1],] * QMat[itemcomb[i,][2],])

        ####### Chi-Sq Test
        x2_pred <- ci.x2(itemcomb[i,][1], itemcomb[i,][2],skills[indc], dataset = data, pval)
        mi_pred <- ci.mi(itemcomb[i,][1], itemcomb[i,][2],skills[indc], dataset = data, pval)
        
        nested_models <- competing_indep_models(itemcomb[i,][1], itemcomb[i,][2],skills[indc], dataset = data)
        
        aic_pred <- ci.aic(nested_models$full, nested_models$reduced, threshold)
        bic_pred <- ci.bic(nested_models$full, nested_models$reduced, threshold)

        x2predicted_dep <- append(x2predicted_dep, x2_pred)
        mipredicted_dep <- append(mipredicted_dep, mi_pred)
        aicpredicted_dep <- append(aicpredicted_dep, aic_pred)
        bicpredicted_dep <- append(bicpredicted_dep, bic_pred)    
    }
    truedep <- rep("0", length(x2predicted_dep))
    
    out <- list(
        truedep = truedep, 
        x2pred = x2predicted_dep,
        mipred = mipredicted_dep,
        aic_pred = aicpredicted_dep,
        bic_pred = bicpredicted_dep)
    return(out)
}