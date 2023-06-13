library(GDINA)
library(progress)
source("last_tests.r")
pb <- progress_bar$new(
  format = "  downloading [:bar] :percent eta: :eta",
  total = 100, clear = FALSE, width= 60)
#=========================================================================
#                  API for users to run the simulation 
#=========================================================================

# run_sim returns results for each sample size.
# 10 simulations are run for each sample
run.analysis <- function(pvals, thresholds, siz,
                    corr =0.3, g=0.3, s=0.3, 
                    QMat = read.csv("q.csv"),
                    model ="DINA"){

    # results hold all the mean performance measures for each simulated sample
    ## results[[i]] = results from sim data sample of siz[i]
    results <- vector("list",length(siz))

    for (i in 1:length(siz)){   
          pb$tick()
        x2perf <- miperf <- as.data.frame(matrix(0,length(pvals),4))
        aicperf <- bicperf <- as.data.frame(matrix(0, length(thresholds),4))

        colnames(x2perf) <- 
                colnames(miperf) <- 
                colnames(aicperf) <- 
                colnames(bicperf) <- c("Acc", "Spe","Sen","AUC")

        for (n in 1:3){
            simdata <- gen.data(QMat = QMat, g=g, s=s, corr=corr,
                            model = model,N=siz[i])

            sim_i <- run.tests(simdata,pvals, thresholds)

            x2perf <- x2perf + sim_i$x2
            miperf <- miperf + sim_i$mi
            aicperf <- aicperf + sim_i$aic
            bicperf <- bicperf + sim_i$bic
        }# end of loop for simulation runs 


        results[[i]] <- list(x2=x2perf/10, mi=miperf/10, 
                          aic=aicperf/10, bic=bicperf/10)
    } # end of loop through siz vector
return(results)
} # end of run_sim function 

run.analysis2 <- function(pvals, thresholds, siz,
                    corr =0.3, g=0.3, s=0.3, 
                    QMat = read.csv("q.csv"),
                    model ="DINA"){

    # results hold all the mean performance measures for each simulated sample
    ## results[[i]] = results from sim data sample of siz[i]
    mean_results <- sd_results <- vector("list",length(siz))
    print("500 runs!!!")
    for (i in 1:length(siz)){   
          pb$tick()
        #x2perf <- miperf <- as.data.frame(matrix(0,length(pvals),4))
        #aicperf <- bicperf <- as.data.frame(matrix(0, length(thresholds),4))
        N = 10 
        x2perf <- miperf <- aicperf <- bicperf <- gaicperf <- xbicperf<- dimtperf <- vector("list",N)
        #colnames(x2perf) <- 
        #        colnames(miperf) <- 
        #        colnames(aicperf) <- 
        #        colnames(bicperf) <- c("Acc", "Spe","Sen","AUC")

        for (n in 1:N){
            simdata <- gen.data(QMat = QMat, g=g, s=s, corr=corr,
                            model = model,N=siz[i])

            sim_i <- run.tests(simdata,pvals, thresholds)

            x2perf[[n]] <- sim_i$x2
            miperf[[n]] <- sim_i$mi
            aicperf[[n]] <- sim_i$aic
            bicperf[[n]] <- sim_i$bic
            gaicperf[[n]] <- sim_i$gaic
            xbicperf[[n]] <- sim_i$xbic
            dimtperf[[n]] <- sim_i$dimt
        }# end of loop for simulation runs 
        
        x2mean <- apply( abind::abind(x2perf, along=3),  1:2, mean)
        mimean <- apply( abind::abind(miperf, along=3),  1:2, mean)
        aicmean <- apply( abind::abind(aicperf, along=3),  1:2, mean)
        bicmean <- apply( abind::abind(bicperf, along=3),  1:2, mean)
        gaicmean <- apply( abind::abind(gaicperf, along=3),  1:2, mean)
        xbicmean <- apply( abind::abind(xbicperf, along=3),  1:2, mean)
        dimtmean <- apply( abind::abind(dimtperf, along=3),  1:2, mean)
        
        x2sd <- apply( abind::abind(x2perf, along=3),  1:2, sd)
        misd <- apply( abind::abind(miperf, along=3),  1:2, sd)
        aicsd <- apply( abind::abind(aicperf, along=3),  1:2, sd)
        bicsd <- apply( abind::abind(bicperf, along=3),  1:2, sd)
        gaicsd <- apply( abind::abind(gaicperf, along=3),  1:2, sd)
        xbicsd <- apply( abind::abind(xbicperf, along=3),  1:2, sd)
        dimtsd <- apply( abind::abind(dimtperf, along=3),  1:2, sd)

        mean_results[[i]] <- list(x2=x2mean, mi=mimean, 
                          aic=aicmean, bic=bicmean, gaic=gaicmean, xbic=xbicmean, dimt=dimtmean)

        sd_results[[i]] <- list(x2=x2sd, mi=misd, 
                          aic=aicsd, bic=bicsd, gaic=gaicsd, xbic = xbicsd, dimt = dimtsd)
    } # end of loop through siz vector
    results = list(mean=mean_results, sd=sd_results)
return(results)
} # end of run_sim function 

#=========================================================================
#                   Function to generate the simulated data 
#=========================================================================
gen.data <- function(QMat = QMat, g=g, s=s, corr=corr,
                            model = model,N=siz){

    # Parameters for generating the simulated data 
    sigma <- matrix(rep(corr, times=25), ncol=5)
    diag(sigma) <- 1
    guess <- rep(g, nrow(QMat));
    slip <- rep(s, nrow(QMat));

    simdata <- CDM::sim.din(N=N,
                        QMat, guess = guess, slip = slip,
                        Sigma=sigma, rule=model)
    
    n_items <- nrow(QMat)
    n_skills <- ncol(QMat)
    itemnames <- paste("X", 1:n_items, sep = "")
    skills <- colnames(QMat)

    XMatrix <- as.data.frame(lapply(as.data.frame(simdata$dat),factor))
    colnames(XMatrix) <- itemnames
    
    # X Matrix + Skill profiles 
    XplusA <- as.data.frame(cbind(simdata$dat, simdata$alpha))
    XplusAdf <- data.frame(lapply(XplusA,factor))
    colnames(XplusAdf) <- c(itemnames,skills)

    data = list(
        n_items = n_items,
        n_skills = n_skills,
        itemnames = itemnames,
        skills = skills,
        QMat = QMat, 
        X = XMatrix,
        XplusA = XplusAdf
    )
    return(data)
}