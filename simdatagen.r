library(GDINA)

sigma <- matrix(rep(0.3, times=25), ncol=5)
diag(sigma) <- 1
QMat <- read.csv("q.csv")
guess <- rep(0.3, nrow(QMat));
slip <- guess
n_students <- 50
sim_data <- function(sigma=sigma, QMat=QMat, guess=guess, slip=slip, model="DINA", n_students){
    itemnames <- rownames(QMat)
    skills <- colnames(QMat)
    n_items <- nrow(QMat)
    n_skills <- ncol(QMat)
    
    simdat <- CDM::sim.din(N=n_students,
                        QMat, guess = guess, slip = guess,
                        Sigma=sigma, rule=model)

    XMatrix <- simdat$dat
    colnames(XMatrix) <- itemnames
    Xdf<- as.data.frame(lapply(as.data.frame(XMatrix),factor))
    colnames(Xdf) <- itemnames

    # X Matrix plus latent skill profile
    XplusA <- as.data.frame(cbind(simdat$dat, simdat$alpha))
    colnames(XplusA) <- c(itemnames,skills)
    XplusAdf <- data.frame(lapply(XplusA,factor))
    itnames <- colnames(XplusAdf)[1:nrow(QMat)]
    data = list(itemnames = itemnames,
                skills = skills,
                n_items = n_items,
                n_skills = n_skills,
                XMatrix = Xdf, 
                QMat = QMat,
                XplusA = XplusAdf, 
                itnames = itnames
    )
    return(data)
}

run_simulation <- function(data,QMat,pval,threshold){

    dep_res = as.data.frame(check.dependencies(data,pval,threshold))
    ind_res = as.data.frame(check.independencies(data,pval,threshold))
    results = rbind(dep_res, ind_res)
    return(results)
}

sim <- function(N, pval, thres){
    QMat <- read.csv("q.csv")
    sigma <- matrix(rep(0.3, times=25), ncol=5)
    diag(sigma) <- 1
    guess <- rep(0.3, nrow(QMat));
    slip <- guess

    x2acc <- vector()
    miacc <- vector()
    aicacc <- vector()
    bicacc <- vector()

    x2spe <- vector()
    mispe <- vector()
    aicspe <- vector()
    bicspe <- vector()

    x2sen <- vector()
    misen <- vector()
    aicsen <- vector()
    bicsen <- vector()

    x2auc <- vector()
    miauc <- vector()
    aicauc <- vector()
    bicauc <- vector()

    for (i in 1:1){
    dat <- sim_data(sigma = sigma,QMat = QMat, guess = guess, slip = slip, n_students = N)
    res = run_simulation(dat,QMat, pval,thres)
    res = as.data.frame(lapply(res, as.factor))

    x2 <- confusionMatrix(na.omit(res$truedep), na.omit(res$x2pred))
    mi <- confusionMatrix(na.omit(res$truedep), na.omit(res$mipred))
    aic <- confusionMatrix(na.omit(res$truedep), na.omit(res$aic_pred))
    bic <-confusionMatrix(na.omit(res$truedep), na.omit(res$bic_pred))
    x2roc <- roc(as.numeric(res$truedep), as.numeric(res$x2pred))
    miroc <- roc(as.numeric(res$truedep), as.numeric(res$mipred))
    aicroc <- roc(as.numeric(res$truedep), as.numeric(res$aic_pred))
    bicroc<- roc(as.numeric(res$truedep), as.numeric(res$bic_pred))

    x2auc <- append(x2auc,auc(x2roc))
    miauc <- append(miauc,auc(miroc))
    aicauc <- append(aicauc,auc(aicroc))
    bicauc <- append(bicauc,auc(bicroc))


    x2acc <- append(x2acc,x2$overall[1])
    miacc <- append(miacc,mi$overall[1])
    aicacc <- append(aicacc,aic$overall[1])
    bicacc <- append(bicacc,bic$overall[1])

    x2spe <- append(x2acc,x2$byClass[2])
    mispe <- append(miacc,mi$byClass[2])
    aicspe <- append(aicacc,aic$byClass[2])
    bicspe <- append(bicacc,bic$byClass[2])

    x2sen <- append(x2acc,x2$byClass[1])
    misen <- append(miacc,mi$byClass[1])
    aicsen <- append(aicacc,aic$byClass[1])
    bicsen <- append(bicacc,bic$byClass[1])

    }
    x2meanacc <- mean(as.numeric(x2acc), na.rm=TRUE)
    mimeanacc <- mean(as.numeric(miacc), na.rm=TRUE)
    aicmeanacc <- mean(as.numeric(aicacc), na.rm=TRUE)
    bicmeanacc <- mean(as.numeric(bicacc), na.rm=TRUE)

    x2meanspe <- mean(as.numeric(x2spe), na.rm=TRUE)
    mimeanspe <- mean(as.numeric(mispe), na.rm=TRUE)
    aicmeanspe <- mean(as.numeric(aicspe), na.rm=TRUE)
    bicmeanspe <- mean(as.numeric(bicspe), na.rm=TRUE)

    x2meansen <- mean(as.numeric(x2sen), na.rm=TRUE)
    mimeansen <- mean(as.numeric(misen), na.rm=TRUE)
    aicmeansen <- mean(as.numeric(aicsen), na.rm=TRUE)
    bicmeansen <- mean(as.numeric(bicsen), na.rm=TRUE)

    x2meanauc <- mean(as.numeric(x2auc), na.rm=TRUE)
    mimeanauc <- mean(as.numeric(miauc), na.rm=TRUE)
    aicmeanauc <- mean(as.numeric(aicauc), na.rm=TRUE)
    bicmeanauc <- mean(as.numeric(bicauc), na.rm=TRUE)

out = list(accuracy = c(x2meanacc, mimeanacc, aicmeanacc, bicmeanacc),
           sensitivity= c(x2meansen, mimeansen, aicmeansen, bicmeansen),
            specificity = c(x2meanspe, mimeanspe, aicmeanspe,bicmeanspe),
            auc=c(x2meanauc, mimeanauc, aicmeanauc, bicmeanauc))
            return(out)
}
