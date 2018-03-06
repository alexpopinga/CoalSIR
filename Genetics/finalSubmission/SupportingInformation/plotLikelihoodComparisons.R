source('likelihoodEstimator.R')

regenerate <- FALSE

###### MAIN ######

if (regenerate) {
    gamma <- 0.3
    beta <- 0.00075
    S0 <- 999
    origin <- 12.7808530307

    tree <- read.tree('VolzSIRgamma_truth.tree')

    gammaVec <- seq(.1,.7,by=.05)

    # Estimate STOCHASTIC coalescent likelihoods for different gammas

    Ntraj <- 10000
    Nensemb <- 10

    llensemb <- list()
    for (e in 1:Nensemb)
        llensemb[[e]] <- rep(0,length(gammaVec))
    Nensemb <- 10

    llensemb <- list()
    for (e in 1:Nensemb)
        llensemb[[e]] <- rep(0,length(gammaVec))

    for (i in 1:length(gammaVec)) {
        for (e in 1:Nensemb) {
            llensemb[[e]][i] <- getCoalescentTreeDensity(tree, beta, gammaVec[i], S0, origin, Ntraj)
        }
    }

    llmean <- rep(0,length(gammaVec))
    llsd <- rep(0, length(gammaVec))
    for (i in 1:length(gammaVec)) {
        thisEnsemble <- rep(0, Nensemb)
        for (e in 1:Nensemb) {
            thisEnsemble[e] <- llensemb[[e]][i]
        }
        llmean[i] <- mean(thisEnsemble)
        llsd[i] <- sd(thisEnsemble)
    }


    # Estimate DETERMINISTIC coalescent likelihoods for different gammas
    
    gammaVecDet <- seq(0.1,0.35,by=0.01)
    lldet <- rep(0, length(gammaVecDet))
    for (i in 1:length(gammaVecDet)) {
        lldet[i] <- getDeterministicCoalescentTreeDensity(tree, beta, gammaVecDet[i], S0, origin)
    }

} else {
    load(file='likelihoodResultsFromR10000_noCorrection.RData')
}

# Load in Java code results for same tree:
df <- read.table('likelihoodResultsFromJava_noCorrection.txt', header=T)
javaGamma <- df$gamma
javaLogP <- df$logP
javaSD <- apply(df[,3:12], 1, sd)


# Create figure
pdf('gammaLikelihoodComparison_noCorrection.pdf', width=7, height=5)

plot(gammaVec, llmean, 'o', #ylim=c(-440,-400),
     xlab=expression(gamma),
     ylab='Log likelihood',
     main='Log likelihoods from simulated tree',
     col='blue')
lines(gammaVec, llmean+2*llsd, lty=2, col='blue')
lines(gammaVec, llmean-2*llsd, lty=2, col='blue')

lines(javaGamma, javaLogP, 'o', col='red')
lines(javaGamma, javaLogP+2*javaSD, lty=2, col='red')
lines(javaGamma, javaLogP-2*javaSD, lty=2, col='red')

#lines(gammaVecDet, lldet, 'o', col='purple')
             
lines(c(0.3,0.3), c(-1e10,1e10), lty=2, col='grey', lwd=2)
#legend('bottomright', inset=.05, c('R','Java (10000)','R (det.)', 'Truth'), lty=c(1,1,1,2), pch=c(1,1,1,NA), lwd=c(1,1,1,2), col=c('blue','red','purple','grey'))
legend('bottomright', inset=.05, c('R 10*(10^4+)','Java 10*(10^4+)', '+/- 2*SD', 'Truth'), lty=c(1,1,2,2), pch=c(1,1,NA,NA), lwd=c(1,1,1,2), col=c('blue','red','black','grey'))
#legend('bottomright', inset=.05, c('R','Java (10000)', '+/- 2*SD', 'Truth'), lty=c(1,1,2,2), pch=c(1,1,NA,NA), lwd=c(1,1,1,2), col=c('blue','red','black','grey'))

dev.off()
