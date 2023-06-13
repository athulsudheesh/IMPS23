# Loading all libraries 
library(bnlearn)
library(caret)
library(pROC)
rm(list = ls())
library(tictoc)
defaultW <- getOption("warn") 
options(warn = -1) 

source("las_simdat.r")
counter <- 0
pvals <- c(0.001, 0.005,0.01,.05,0.1,0.2,0.3, 0.4, 0.5, 0.7, 0.8, 0.9)
thresholds <- c(0.0001, 0.001, 0.005,0.01,.05,0.1,0.2,0.3, 0.4, 0.5, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 2, 2.5, 3, 5, 10, 20, 150) 


siz <- c(200, 500,1000,2000)
tic("Run Time")
res = run.analysis2(pvals, thresholds, siz)
toc()




get_vals <- function(i,j){
x2 <- rbind(res$mean[[1]]$x2[i,],
                res$mean[[2]]$x2[i,],
                res$mean[[3]]$x2[i,])
                #res$mean[[4]]$x2[i,],
                #res$mean[[5]]$x2[i,])
                #res$mean[[6]]$x2[i,])
    mi <- rbind(res$mean[[1]]$mi[i,],
                res$mean[[2]]$mi[i,],
                res$mean[[3]]$mi[i,])
                #res$mean[[4]]$mi[i,],
                #res$mean[[5]]$mi[i,])
                #res$mean[[6]]$mi[i,])

    aic <- rbind(res$mean[[1]]$aic[j,],
                res$mean[[2]]$aic[j,],
                res$mean[[3]]$aic[j,])
                #res$mean[[4]]$aic[j,],
                #res$mean[[5]]$aic[j,])
                #res$mean[[6]]$aic[j,])

    bic <- rbind(res$mean[[1]]$bic[j,],
                res$mean[[2]]$bic[j,],
                res$mean[[3]]$bic[j,])
                #res$mean[[4]]$bic[j,],
                #res$mean[[5]]$bic[j,])
                #res$mean[[6]]$bic[j,])

    write.csv(x2, "x2.csv")
    write.csv(mi, "mi.csv")
    write.csv(aic, "aic.csv")
    write.csv(bic, "bic.csv")

    x2sd <- rbind(res$sd[[1]]$x2[i,],
                res$sd[[2]]$x2[i,],
                res$sd[[3]]$x2[i,])
                #res$sd[[4]]$x2[i,],
                #res$sd[[5]]$x2[i,])
                #res$sd[[6]]$x2[i,])
    misd <- rbind(res$sd[[1]]$mi[i,],
                res$sd[[2]]$mi[i,],
                res$sd[[3]]$mi[i,])
                #res$sd[[4]]$mi[i,],
                #res$sd[[5]]$mi[i,])
                #res$sd[[6]]$mi[i,])

    aicsd <- rbind(res$sd[[1]]$aic[j,],
                res$sd[[2]]$aic[j,],
                res$sd[[3]]$aic[j,])
                #res$sd[[4]]$aic[j,],
                #res$sd[[5]]$aic[j,])
                #res$sd[[6]]$aic[j,])

    bicsd <- rbind(res$sd[[1]]$bic[j,],
                res$sd[[2]]$bic[j,],
                res$sd[[3]]$bic[j,])
                #res$sd[[4]]$bic[j,],
                #res$sd[[5]]$bic[j,])
                #res$sd[[6]]$bic[j,])

    write.csv(x2sd, "x2sd.csv")
    write.csv(misd, "misd.csv")
    write.csv(aicsd, "aicsd.csv")
    write.csv(bicsd, "bicsd.csv")

}

get_vals(4,13)
library(pracma)
get_auc <- function(){
    auc <- vector()
    aucsd <- vector()
    for (i in 1:4){

    x2Spe <- res$mean[[i]]$x2[,2]
    miSpe <- res$mean[[i]]$mi[,2]
    aicSpe <- res$mean[[i]]$aic[,2]
    bicSpe <- res$mean[[i]]$bic[,2]
    gaicSpe <- res$mean[[i]]$gaic[,2]
    xbicSpe <- res$mean[[i]]$xbic[,2]
    dimtSpe <- res$mean[[i]]$dimt[,2]

    x2Spesd <- res$sd[[i]]$x2[,2]
    miSpesd <- res$sd[[i]]$mi[,2]
    aicSpesd <- res$sd[[i]]$aic[,2]
    bicSpesd <- res$sd[[i]]$bic[,2]
    gaicSpesd <- res$sd[[i]]$gaic[,2]
    xbicSpesd <- res$sd[[i]]$xbic[,2]
    dimtSpesd <- res$sd[[i]]$dimt[,2]

    x2Sen <- res$mean[[i]]$x2[,3]
    miSen <- res$mean[[i]]$mi[,3]
    aicSen <- res$mean[[i]]$aic[,3]
    bicSen <- res$mean[[i]]$bic[,3]
    gaicSen <- res$mean[[i]]$gaic[,3]
    xbicSen <- res$mean[[i]]$xbic[,3]
    dimtSen <- res$mean[[i]]$dimt[,3]

    x2Sensd <- res$sd[[i]]$x2[,3]
    miSensd <- res$sd[[i]]$mi[,3]
    aicSensd <- res$sd[[i]]$aic[,3]
    bicSensd <- res$sd[[i]]$bic[,3]
    gaicSensd <- res$sd[[i]]$gaic[,3]
    xbicSensd <- res$sd[[i]]$xbic[,3]
    dimtSensd <- res$sd[[i]]$dimt[,3]

    x2 <- trapz(1 - x2Spe, x2Sen)
    mi <- trapz(1 - miSpe, miSen)
    aic <- trapz(1 - aicSpe, aicSen)
    bic <- trapz(1 - bicSpe, bicSen)
    gaic <- trapz(1 - gaicSpe, gaicSen)
    xbic <- trapz(1 - xbicSpe, xbicSen)
    dimt <- trapz(1 - dimtSpe, dimtSen)

    x2sd <- trapz(x2Spesd, x2Sensd)
    misd <- trapz(miSpesd, miSensd)
    aicsd <- trapz(aicSpesd, aicSensd)
    bicsd <- trapz(bicSpesd, bicSensd)
    gaicsd <- trapz(gaicSpesd, gaicSensd)
    xbicsd <- trapz(xbicSpesd, xbicSensd)
    dimtsd <- trapz(dimtSpesd, dimtSensd)
    auc <- rbind(auc, c(x2,mi,aic,bic, gaic, xbic, dimt))
    aucsd <- rbind(aucsd, c(x2sd,misd,aicsd,bicsd, gaicsd, xbicsd, dimtsd))
    }
    colnames(auc) = c("X2","MI","AIC","BIC", "GAIC", "XBIC", "DIMT")
    colnames(aucsd) = c("X2","MI","AIC","BIC","GIAC", "XBIC", "DIMT")
    
    write.csv(auc, "auc.csv")
    write.csv(aucsd, "aucsd.csv")
    
    #out = list(auc=auc, sd = aucsd)
    #return(out)
    #write.csv(x2, "x2roc.csv")
    #write.csv(mi, "miroc.csv")
    #write.csv(aic, "aicroc.csv")
    #write.csv(bic, "bicroc.csv")

    #write.csv(x2sd, "x2rocsd.csv")
    #write.csv(misd, "mirocsd.csv")
    #write.csv(aicsd, "aicrocsd.csv")
    #write.csv(bicsd, "bicrocsd.csv")
}
get_auc()
get_roc()
a <- get_roc()

get_roc(3)

siz <- c(200, 500,1000,2000, 5000)
auc <-function(tpr,fpr){
    dFPR <- c(diff(fpr),0)
    dTPR <- c(diff(tpr),0)
    AUC <- sum(fpr * tpr) + sum(dFPR * dTPR)/2
}

get_auc <- function(i){
    x2_fpr <- 1 - res$mean[[i]]$x2[,2]
    mi_fpr <- 1 - res$mean[[i]]$mi[,2]
    aic_fpr <- 1 - res$mean[[i]]$aic[,2]
    bic_fpr <- 1- res$mean[[i]]$bic[,2]
    print(x2_fpr)
    print(" ")
    x2_tpr <- res$mean[[i]]$x2[,3]
    mi_tpr <- res$mean[[i]]$mi[,3]
    aic_tpr <- res$mean[[i]]$aic[,3]
    bic_tpr <- res$mean[[i]]$bic[,3]

    x2_roc <- data.frame(TPR=x2_tpr, FPR = x2_fpr)
    x2_roc <- transform(x2_roc, 
                    dFPR = c(diff(FPR), 0),
                    dTPR = c(diff(TPR), 0))

    print(x2_tpr)
        print(" ")
    x2_fpr_sd <- res$sd[[i]]$x2[,2]
    mi_fpr_sd <- res$sd[[i]]$mi[,2]
    aic_fpr_sd <- res$sd[[i]]$aic[,2]
    bic_fpr_sd <- res$sd[[i]]$bic[,2]

    x2_tpr_sd <- res$sd[[i]]$x2[,3]
    mi_tpr_sd <- res$sd[[i]]$mi[,3]
    aic_tpr_sd <- res$sd[[i]]$aic[,3]
    bic_tpr_sd <- res$sd[[i]]$bic[,3]
    
    x2AUC <- with(x2_roc,auc(x2_tpr,x2_fpr))
    return(x2AUC)
}

save(res,file="500res.Rdata")
save(res,file="abstractResults.Rdata")
get_roc <- function(i){


    x2 <- res$mean[[i]]$x2[,c(2,3)]
    mi <- res$mean[[i]]$mi[,c(2,3)]
    aic <- res$mean[[i]]$aic[,c(2,3)]
    bic <- res$mean[[i]]$bic[,c(2,3)]

    x2sd <- res$sd[[i]]$x2[,c(2,3)]
    misd <- res$sd[[i]]$mi[,c(2,3)]
    aicsd <- res$sd[[i]]$aic[,c(2,3)]
    bicsd <- res$sd[[i]]$bic[,c(2,3)]

    
    #out = list(auc=auc, sd = aucsd)
    #return(out)
    write.csv(x2, "x2roc.csv")
    write.csv(mi, "miroc.csv")
    write.csv(aic, "aicroc.csv")
    write.csv(bic, "bicroc.csv")

    write.csv(x2sd, "x2rocsd.csv")
    write.csv(misd, "mirocsd.csv")
    write.csv(aicsd, "aicrocsd.csv")
    write.csv(bicsd, "bicrocsd.csv")
}

get_roc(4)

plot([rand(10), rand(20)], color = [:black :orange], line = (:dot, 4), marker = ([:hex :d], 12, 0.8, Plots.stroke(3, :gray)))
