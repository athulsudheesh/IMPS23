library(pROC)
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