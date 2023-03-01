# Loading all libraries 
library(bnlearn)
library(caret)
library(pROC)
rm(list = ls())
library(tictoc)
defaultW <- getOption("warn") 
options(warn = -1) 

source("las_simdat.r")

pvals <- c(0.001, 0.005,0.01,.05,0.1,0.2,0.3, 0.4, 0.5, 0.7, 0.8, 0.9)
thresholds <- c(0.001, 0.005,0.01,.05,0.1,0.2,0.3, 0.4, 0.5, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 2, 2.5, 3, 5, 10, 20) 


siz <- c(20,50,100,200,500,1000,2000, 5000)
tic("Run Time")
res = run.analysis2(pvals, thresholds, siz)
toc()
get_vals <- function(i,j){
x2 <- rbind(res$mean[[1]]$x2[i,],
                res$mean[[2]]$x2[i,],
                res$mean[[3]]$x2[i,],
                res$mean[[4]]$x2[i,],
                res$mean[[5]]$x2[i,])
                #res$mean[[6]]$x2[i,])
    mi <- rbind(res$mean[[1]]$mi[i,],
                res$mean[[2]]$mi[i,],
                res$mean[[3]]$mi[i,],
                res$mean[[4]]$mi[i,],
                res$mean[[5]]$mi[i,])
                #res$mean[[6]]$mi[i,])

    aic <- rbind(res$mean[[1]]$aic[j,],
                res$mean[[2]]$aic[j,],
                res$mean[[3]]$aic[j,],
                res$mean[[4]]$aic[j,],
                res$mean[[5]]$aic[j,])
                #res$mean[[6]]$aic[j,])

    bic <- rbind(res$mean[[1]]$bic[j,],
                res$mean[[2]]$bic[j,],
                res$mean[[3]]$bic[j,],
                res$mean[[4]]$bic[j,],
                res$mean[[5]]$bic[j,])
                #res$mean[[6]]$bic[j,])

    write.csv(x2, "x2.csv")
    write.csv(mi, "mi.csv")
    write.csv(aic, "aic.csv")
    write.csv(bic, "bic.csv")

    x2sd <- rbind(res$sd[[1]]$x2[i,],
                res$sd[[2]]$x2[i,],
                res$sd[[3]]$x2[i,],
                res$sd[[4]]$x2[i,],
                res$sd[[5]]$x2[i,])
                #res$sd[[6]]$x2[i,])
    misd <- rbind(res$sd[[1]]$mi[i,],
                res$sd[[2]]$mi[i,],
                res$sd[[3]]$mi[i,],
                res$sd[[4]]$mi[i,],
                res$sd[[5]]$mi[i,])
                #res$sd[[6]]$mi[i,])

    aicsd <- rbind(res$sd[[1]]$aic[j,],
                res$sd[[2]]$aic[j,],
                res$sd[[3]]$aic[j,],
                res$sd[[4]]$aic[j,],
                res$sd[[5]]$aic[j,])
                #res$sd[[6]]$aic[j,])

    bicsd <- rbind(res$sd[[1]]$bic[j,],
                res$sd[[2]]$bic[j,],
                res$sd[[3]]$bic[j,],
                res$sd[[4]]$bic[j,],
                res$sd[[5]]$bic[j,])
                #res$sd[[6]]$bic[j,])

    write.csv(x2sd, "x2sd.csv")
    write.csv(misd, "misd.csv")
    write.csv(aicsd, "aicsd.csv")
    write.csv(bicsd, "bicsd.csv")

}

get_vals(4,13)

get_roc <- function(i){

    x2 <- res$mean[[i]]$x2[,c(2,3)]
    mi <- res$mean[[i]]$mi[,c(2,3)]
    aic <- res$mean[[i]]$aic[,c(2,3)]
    bic <- res$mean[[i]]$bic[,c(2,3)]

    x2sd <- res$sd[[i]]$x2[,c(2,3)]
    misd <- res$sd[[i]]$mi[,c(2,3)]
    aicsd <- res$sd[[i]]$aic[,c(2,3)]
    bicsd <- res$sd[[i]]$bic[,c(2,3)]

    write.csv(x2, "x2roc.csv")
    write.csv(mi, "miroc.csv")
    write.csv(aic, "aicroc.csv")
    write.csv(bic, "bicroc.csv")

    write.csv(x2sd, "x2rocsd.csv")
    write.csv(misd, "mirocsd.csv")
    write.csv(aicsd, "aicrocsd.csv")
    write.csv(bicsd, "bicrocsd.csv")
}

get_roc(8)

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
    roc_df <- transform(x2_roc, 
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

    x2_dFPR <- c(diff(x2_fpr),0)
    x2_dTPR <- c(diff(x2_tpr),0)
    print(x2_dFPR)
    print(x2_dTPR)
    x2AUC <- (sum(x2_fpr * x2_tpr) + sum(x2_dFPR * x2_dTPR))/2
    return(x2AUC)
}
