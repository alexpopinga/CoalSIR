simTraj <- function(R0, gamma, S0, T, Nsteps) {

    dt <- T/(Nsteps-1)
    
    beta <- R0*gamma/S0
    
    t <- 0
    I <- 1
    S <- S0

    res <- list()
    res$S <- rep(0,Nsteps)
    res$I <- rep(0,Nsteps)
    res$t <- rep(0,Nsteps)

    i <- 1
    
    while (i<=Nsteps) {
        
        ainfect <- beta*S*I
        aremove <- gamma*I
        a0 <- ainfect + aremove

        if (a0==0) {
            t <- Inf
        }
        else {
            t <- t + rexp(1, a0)
        }

        while (t>dt*(i-1)) {
            res$S[i] <- S
            res$I[i] <- I
            res$t[i] <- dt*(i-1)
            i <- i + 1

            if (i>Nsteps) {
                return(res)
            }
        }

        if (runif(1)*a0 < ainfect) {
            I <- I + 1
            S <- S - 1
        } else {
            I <- I - 1
        }
    }

    return(res)
}

getDeterministic <- function(R0, gamma, S0, T, Nsteps) {

    dt <- T/(Nsteps-1)
    
    beta <- R0*gamma/S0

    t <- 0
    I <- 1
    S <- S0
    i <- 1

    for (i in 1:(Nsteps-1)) {

        Sprime <- S[i]
        Iprime <- I[i]

        for (iter in 1:3) {
            Sprime <- S[i] + 0.5*dSdt(Sprime, Iprime, beta)*dt
            Iprime <- I[i] + 0.5*dIdt(Sprime, Iprime, beta, gamma)*dt
        }

        S[i+1] <- 2*Sprime - S[i]
        I[i+1] <- 2*Iprime - I[i]

        t[i+1] <- t[i] + dt
    }

    res <- list()
    res$t <- t
    res$I <- I

    return(res)
}

dSdt <- function(S,I,beta) {
    return(-beta*S*I)
}

dIdt <- function(S,I,beta,gamma) {
    return(beta*S*I - gamma*I)
}

getResidue <- function(stochTraj, detTraj, offset) {

    r <- 0
    for (i in min(1, 1+offset):max(Nsteps, Nsteps+offset)) {
        j <- i + offset
        
        detI <- 0
        if (j>0 && j<=Nsteps) {
            detI <- detTraj$I[j]
        } else {
            detI <- 0
        }

        if (i>0 && i<=Nsteps) {
            stochI <- stochTraj$I[i]
        } else {
            stochI <- 0
        }
        
        r <- r + abs(sqrt(stochI) - sqrt(detI))
    }
        
    return(r)
}

findOptimalOffset <- function(stochTraj, detTraj, coarseGrain=1) {

    residues <- NULL
    
    Nsteps <- length(detTraj$t)
    offsets <- seq(-Nsteps, Nsteps, by=coarseGrain)
    for (offset in offsets) {
        residues <- append(residues, getResidue(stochTraj, detTraj, offset))
    }

    res <- list()
    res$r <- residues
    res$o <- offsets

    optoIdx <- which.min(residues)
    res$opto <- offsets[optoIdx]
    res$minr <- residues[optoIdx]
    
    return(res)
}

## MAIN ##

R0 <- 1.5
gamma <- 0.3
S0 <- 999
T <- 150
Nsteps <- 2000
dt <- T/(Nsteps-1)
coarseGrain <- 20

infectMin <- 300

det <- getDeterministic(R0, gamma, S0, T, Nsteps)

traj <- list()
residual <- NULL
residualAdjusted <- NULL
offset <- NULL
for (i in 1:500) {
    cat(paste(i,"\n"))
    while (TRUE) {
        traj[[i]] <- simTraj(R0, gamma, S0, T, Nsteps)
        if (S0-traj[[i]]$S[Nsteps] >= infectMin)
            break;
    }
    residual[i] <- getResidue(traj[[i]], det, 0)
    opt <- findOptimalOffset(traj[[i]], det, coarseGrain=coarseGrain)
    residualAdjusted[i] <- opt$minr
    offset[i] <- dt*opt$opto
}

# It's useful to peek at these when setting infectMin
epiSizes <- NULL
for (i in 1:length(traj)) {epiSizes[i] <- 999 - traj[[i]]$S[length(traj[[i]]$S)]}

## Figure

pdf('originFitFigure.pdf', width=5, height=6)

par(mfrow=c(2,1))
par(mar=c(3,3,2,2))
par(mgp=c(1.5,0.5,0))

plot(1, type="n", xlim=c(-20,150), ylim=c(0,150),
     xlab="Time (after true origin)", ylab="Prevalence", main="(a)")

for (i in 1:length(traj)) {
    lines(traj[[i]]$t, traj[[i]]$I, col=rgb(0,0,0,0.1))
}

for (i in 1:length(traj)) {
    off <- offset[i] + (runif(1)+0.5)*dt*coarseGrain
    lines(det$t-off, det$I, col=rgb(0,0,1,0.1))
}

lines(det$t, det$I, col='red', lwd=2)

legend('topright', inset=0.03,
       c('Stochastic epidemic', 'Deterministic', 'Fitted deterministic'),
       lty=1, lwd=2, col=c('black','red','blue'), cex=0.8)

hres <- hist(residual, breaks=seq(0,10000, length.out=50), plot=F)
hresAdj <- hist(residualAdjusted, breaks=seq(0,10000, length.out=50), plot=F)

plot(hres$mids, hres$density, 'l', lwd=2,
     col='red', ylim=c(0, 0.0015), xlim=c(0,5e3),
     xlab="Residual", ylab="Sample density", main="(b)")
lines(hresAdj$mids, hresAdj$density, lwd=2, col='blue')

legend('topright', inset=0.03,
       c('No origin fitting', 'Origin fitted'),
       lwd=2, lty=1, col=c('red','blue'), cex=0.8)

dev.off()
