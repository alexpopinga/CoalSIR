# Monte Carlo estimate of SIR parameter likelihoods given tree.
# This script requires PEA (github.com/tgvaughan/PEA).

require(pea)

# Integrate the SIR ODEs over the specified domain
deterministicSIR <- function(beta, gamma, S0, T, steps, maxIter=3) {

    S <- S0
    I <- 1
    t <- 0

    dt <- T/(steps-1)

    for (tidx in 2:steps) {

        Sp <- S[tidx-1]
        Ip <- I[tidx-1]

        for (iter in 1:maxIter) {
            dSdt <- -beta*Sp*Ip
            dIdt <- beta*Sp*Ip - gamma*Ip

            Sp <- S[tidx-1] + 0.5*dt*dSdt
            Ip <- I[tidx-1] + 0.5*dt*dIdt
        }

        S[tidx] <- 2*Sp - S[tidx-1]
        I[tidx] <- 2*Ip - I[tidx-1]
        t[tidx] <- tidx*dt

    }

    res <- list()
    res$S <- S
    res$I <- I
    res$t <- t

    return(res)
}



# Stochastic simulation of an SIR trajectory
simSIRTraj <- function(beta, gamma, S0, T) {

    S <- S0
    I <- 1
    t <- 0

    idx <- 1
    while (TRUE) {

        ainfect <- beta*S[idx]*I[idx]
        aremove <- gamma*I[idx]
        atot <- ainfect + aremove

        if (atot==0) {
            tnext <- Inf
        } else {
            tnext <- t[idx] + rexp(1, atot)
        }

        if (tnext>T)
            break

        t[idx+1] <- tnext

        if (runif(1,min=0,max=atot)<ainfect) {
            S[idx+1] <- S[idx] - 1
            I[idx+1] <- I[idx] + 1
        } else {
            S[idx+1] <- S[idx]
            I[idx+1] <- I[idx] - 1
        }

        idx <- idx + 1
    }

    res <- list()
    res$S <- S
    res$I <- I
    res$t <- t

    return(res)
}

# Population size function used in coalescent calculation
getEffectivePopSize <- function(beta,S,I,corrected=TRUE) {
    if (correct)
        return ((I-1)/(2*beta*S))
    else
        return (I/(2*beta*S))
}

# Calculate probability density of tree event sequence given trajectory
getCoalescentTreeDensityForTraj <- function(treeEvents, traj, beta, origin, corrected=TRUE) {

    # Check for trajectories shorter than tree
    if (corrected)
        bottleNecks <- which(traj$I<=1)
    else
        bottleNecks <- which(traj$I<1)
    if (length(bottleNecks)>0 && (traj$t[bottleNecks[1]] < origin)) {
        return (-Inf)
    }

    logDensity <- 0
    
    Svec <- rev(traj$S)
    Ivec <- rev(traj$I)
    tvec <- origin - rev(traj$t)
    
    tidx <- 1
    t <- 0
    
    for (idx in 2:length(treeEvents$heights)) {
        while (treeEvents$heights[idx]>tvec[tidx]) {
            
            rate <- 1.0/getEffectivePopSize(beta, Svec[tidx], Ivec[tidx], corrected)
            
            # Waiting time contribution
            logDensity <- logDensity + -(tvec[tidx]-t)*choose(treeEvents$lineages[idx-1],2)*rate
            
            t <- tvec[tidx]
            tidx <- tidx + 1
        }

        rate <- 1.0/getEffectivePopSize(beta,Svec[tidx],Ivec[tidx], corrected)
        
        # Waiting time contribution
        logDensity <- logDensity +
            -(treeEvents$heights[idx]-t)*choose(treeEvents$lineages[idx-1],2)*rate
        
        # Event time contribution (only if coalescence)
        if (treeEvents$lineages[idx]<treeEvents$lineages[idx-1])
            logDensity <- logDensity + log(rate)
        
        t <- treeEvents$heights[idx]
    }

    return (logDensity)
}


# Likelihood estimate from the deterministic solution
getDeterministicCoalescentTreeDensity <- function(tree, beta, gamma, S0, origin, steps=1000, corrected=TRUE) {

    treeEvents <- getNodeHeights(tree)
    traj <- deterministicSIR(beta, gamma, S0, origin, steps)
    logDensity <- getCoalescentTreeDensityForTraj(tree, traj, beta, gamma, origin, corrected)

    return (logDensity)
}


# Stochastic likelihood estimate using a number of simulated SIR epidemics
getStochasticCoalescentTreeDensity <- function(tree, beta, gamma, S0, origin, Ntraj=1000, corrected=TRUE) {

    logDensity <- rep(0,Ntraj)

    treeEvents <- getNodeHeights(tree)

    trajIdx <- 0
    goodTrajIdx <- 0
    while (goodTrajIdx<Ntraj) {
        trajIdx <- trajIdx + 1
        if (trajIdx%%100 == 0)
            cat(paste("beta:",beta,"gamma:",gamma,"S0:",S0,"origin:",origin,"Trajectory",goodTrajIdx,"of",Ntraj,"\n"))

        traj <- simSIRTraj(beta, gamma, S0, origin)

        thisLogDensity <- getCoalescentTreeDensityForTraj(treeEvents, traj, beta, gamma, origin, corrected)
        
        # Increment goodTrajIdx when trajectory encompasses tree
        if (thisLogDensity>-Inf) {
            goodTrajIdx <- goodTrajIdx + 1
        }
        
        logDensity[trajIdx] <- thisLogDensity
    }

    # Calculate log of mean of densities
    
    maxLogDensity <- max(logDensity)

    logDensityShifted <- logDensity - maxLogDensity
    scaledDensities <- exp(logDensityShifted)
    meanScaledDensity <- mean(scaledDensities)

    return(log(meanScaledDensity) + maxLogDensity)
}
