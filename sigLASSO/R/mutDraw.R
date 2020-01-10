mutDraw <-
function (mutationNum, simSpec, noise) 
{
    noise_permute <- sample(seq(from = 1, to = 96), floor(runif(1, 
        min = 1, max = 25)))
    noise_random <- abs(rnorm(length(noise_permute)))
    noise_weights <- noise_random/sum(noise_random) * noise
    for (i in 1:length(noise_permute)) {
        simSpec[noise_permute[i]] <- simSpec[noise_permute[i]] + 
            noise_weights[i]
    }
    simAccm <- rep(NA, 96)
    simAccm[1] <- simSpec[1]
    for (i in 2:96) {
        simAccm[i] <- simAccm[i - 1] + simSpec[i]
    }
    simRandom <- runif(mutationNum)
    simMut <- rep(0, 96)
    simMut[1] <- length(which(simRandom <= simAccm[1]))
    for (i in 2:96) {
        simMut[i] = length(which(simRandom <= simAccm[i] & simRandom > 
            simAccm[i - 1]))
    }
    return(simMut)
}
