rm(list = ls())
library(bnlearn)
source("tests.r")
source("simdatagen.r")
library(caret)

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

    for (i in 1:10){
    dat <- sim_data(sigma = sigma,QMat = QMat, guess = guess, slip = slip, n_students = N)
    res = run_simulation(dat,QMat, pval,thres)
    res = as.data.frame(lapply(res, as.factor))

    x2 <- confusionMatrix(res$truedep, res$x2pred)
    mi <- confusionMatrix(res$truedep, res$mipred)
    aic <- confusionMatrix(res$truedep, res$aic_pred)
    bic <-confusionMatrix(res$truedep, res$bic_pred)
    
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

out = list(accuracy = c(x2meanacc, mimeanacc, aicmeanacc, bicmeanacc),
           sensitivity= c(x2meansen, mimeansen, aicmeansen, bicmeansen),
            specificity = c(x2meanspe, mimeanspe, aicmeanspe,bicmeanspe))
            return(out)
}

res1000 <- sim(1000,0.05,1)
res500 <- sim(500,0.05,1)
res200 <- sim(200,0.05,1)
res100 <- sim(200,0.05,1)
res50 <- sim(50,0.05,1)
res20 <- sim(20,0.05,1)
siz <- c(1000,500,200,100,50,20)

acc<- rbind(res1000$accuracy, res500$accuracy, res200$accuracy, res100$accuracy, res50$accuracy, res20$accuracy)
write.csv(acc, "acc_at0.05nd1.csv")
plot(siz, acc[,1], type="l")

#res = as.data.frame(check.dependencies(dat,0.05,0.05))

#res2 = as.data.frame(check.independencies(dat,0.05,1))

#rbind(res,res2)

#run_simulation <- function(dat, pval, threshold){
#    dep_res = as.data.frame(check.dependencies(dat,pval,threshold))
#    ind_res = as.data.frame(check.independencies(dat,pval,threshold))
#    results = rbind(dep_res, ind_res)
#    return(results)
#}
