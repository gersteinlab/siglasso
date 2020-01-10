randomMutGen <-
function (mutationNum, simNum, noise) 
{
    min_weight <- 0.05
    min_weight <- 0.02
    simSig <- sample(seq(2, 31), simNum)
    simPercentPool <- sort(runif(simNum - 1, max = 1 - noise - 
        simNum * min_weight), decreasing = T)
    simSigProfile <- rep(0, 30)
    simSpec <- sig[, simSig[1]] * (1 - noise - simNum * min_weight - 
        simPercentPool[1] + min_weight)
    simSigProfile[simSig[1] - 1] <- 1 - noise - simNum * min_weight - 
        simPercentPool[1] + min_weight
    simSpec <- simSpec + sig[, simSig[simNum]] * (simPercentPool[simNum - 
        1] + min_weight)
    simSigProfile[simSig[simNum] - 1] <- simPercentPool[simNum - 
        1] + min_weight
    if (simNum > 2) {
        for (i in 2:(simNum - 1)) {
            simSpec <- simSpec + sig[, simSig[i]] * (simPercentPool[i - 
                1] - simPercentPool[i] + min_weight)
            simSigProfile[simSig[i] - 1] <- simPercentPool[i - 
                1] - simPercentPool[i] + min_weight
        }
    }
    return(list(mutDraw(mutationNum, simSpec, noise), simSigProfile, 
        simSpec))
}
