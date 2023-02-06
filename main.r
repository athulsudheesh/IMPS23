# Simulation Studies for testing Conditional Independence in CDMs ===================================
#  Ref.   Lim, Y. S., & Drasgow, F. (2019). Conditional independence and 
#           dimensionality of cognitive diagnostic models: a test for model 
#           fit. Journal of Classification, 36, 295-305.


sigma <- matrix(rep(0.3, times=25), ncol=5) # Covariance matrix with unit variance and rho = 0.3
diag(sigma) <- 1
QMat <- read.csv("q.csv")
guess <- rep(0.3, nrow(QMat)); # guess and slip ~ Uniform(0,0.3)

library(GDINA);

simdat <- CDM::sim.din(N=500, QMat,
                guess = guess,slip = guess,
                Sigma=sigma, rule="DINA") 
                              
# CDM as a BayesNet =========================================
library(bnlearn)
a = empty.graph(c("x1", "x2", "x3","x4", "x5",
                    "x6", "x7", "x8", "x9", "x10",
                    "x11", "x12", "x13", "x14", "x15",
                    "x16", "x17", "x18", "x19", "x20",
                    "a1", "a2", "a3", "a4", "a5"))

adjMatrix <- as.matrix(read.csv("newadj.csv"))
row.names(adjMatrix) <- c("a1", "a2", "a3", "a4", "a5",
                    "x1", "x2", "x3","x4", "x5",
                    "x6", "x7", "x8", "x9", "x10",
                    "x11", "x12", "x13", "x14", "x15",
                    "x16", "x17", "x18", "x19", "x20")


colnames(adjMatrix)<- c("a1", "a2", "a3", "a4", "a5",
                    "x1", "x2", "x3","x4", "x5",
                    "x6", "x7", "x8", "x9", "x10",
                    "x11", "x12", "x13", "x14", "x15",
                    "x16", "x17", "x18", "x19", "x20")

amat(a) <- adjMatrix

# Plotting the bn to visually check for correctness
library(Rgraphviz)
graphviz.plot(a)


## Conditional Independence Testing using Chi-Sq
chi_test <- function(n){
cntr = 0 
for (i in 1:100){

    # generate simulated data with given QMatrix, guess, and slip 
    simulate_data <- CDM::sim.din(N=n,
                                    QMat,
                                    guess = guess, slip = guess,
                                    Sigma=sigma, rule="DINA");


    # generating an empty bayesian network 
    a = empty.graph(c("x1", "x2", "x3","x4", "x5",
                    "x6", "x7", "x8", "x9", "x10",
                    "x11", "x12", "x13", "x14", "x15",
                    "x16", "x17", "x18", "x19", "x20",
                    "a1", "a2", "a3", "a4", "a5"))


    # Specifying the conditional independence relations based on given Q Matrix 
    amat(a) <- adjMatrix
    

    # preparing the data for training bayesian network 
    fulldat <- cbind(simulate_data$dat, simulate_data$alpha)
    colnames(fulldat)<- c("x1", "x2", "x3","x4", "x5",
                    "x6", "x7", "x8", "x9", "x10",
                    "x11", "x12", "x13", "x14", "x15",
                    "x16", "x17", "x18", "x19", "x20",
                    "a1", "a2","a3","a4","a5")
    fdf <- as.data.frame(fulldat)
    factoreddata <- data.frame(lapply(fdf, factor))

    # Checking conditional independence using Chi-Squared 
    res <- arc.strength(a, data = factoreddata, criterion = "x2" )                
    if (any(res$strength > 0.05))
    {
        cntr = cntr + 1 
    }
};
return(cntr)
}

# Computing Type I error for different sample size 
chi_val <-vector("list", 5)
j <-1 
for (i in seq(100,1000, by = 250)){
    chi_val[j] <- chi_test(i) 
    j = j + 1
}

mi_test <- function(n){
cntr = 0 
for (i in 1:100){

    # generate simulated data with given QMatrix, guess, and slip 
    simulate_data <- CDM::sim.din(N=n,
                                    QMat,
                                    guess = guess, slip = guess,
                                    Sigma=sigma, rule="DINA");


    # generating an empty bayesian network 
    a = empty.graph(c("x1", "x2", "x3","x4", "x5",
                    "x6", "x7", "x8", "x9", "x10",
                    "x11", "x12", "x13", "x14", "x15",
                    "x16", "x17", "x18", "x19", "x20",
                    "a1", "a2", "a3", "a4", "a5"))


    # Specifying the conditional independence relations based on given Q Matrix 
    amat(a) <- adjMatrix
    

    # preparing the data for training bayesian network 
    fulldat <- cbind(simulate_data$dat, simulate_data$alpha)
    colnames(fulldat)<- c("x1", "x2", "x3","x4", "x5",
                    "x6", "x7", "x8", "x9", "x10",
                    "x11", "x12", "x13", "x14", "x15",
                    "x16", "x17", "x18", "x19", "x20",
                    "a1", "a2","a3","a4","a5")
    fdf <- as.data.frame(fulldat)
    factoreddata <- data.frame(lapply(fdf, factor))

    # Checking conditional independence using Chi-Squared 
    res <- arc.strength(a, data = factoreddata, criterion = "mi" )                
    if (any(res$strength > 0.05))
    {
        cntr = cntr + 1 
    }
};
return(cntr)
}

mi_val <-vector("list", 5)
j <-1 
for (i in seq(100,1000, by = 250)){
    mi_val[j] = mi_test(i) 
    j = j + 1
}

# Type I error plot for Chi-Sq and MI 
plot(c(100,350,600,850,1100), c(100,96,56,18,0),type='l', xlab ="Sample Size", ylab="Type 1 error", 
             sub="No. of simulations = 100", col="red")
lines(c(100,350,600,850,1100), c(100,92,44,14,0), col ="darkgreen")
legend("topright", legend = c("Chi-sq", "MI"), lwd = 3,col = c("red", "darkgreen"))
