loglikelihood_cal <-
function (weights, observation, signatures) 
{
    weights <- apply(weights, 2, function(x) x/sum(x))
    P <- log10(signatures %*% weights + 10^(-10))
    if (is.vector(P)) {
        L = 1
        for (i in length(P)) {
            L = L * P[i]^observation[i]
        }
    }
    else {
        L = rep(0, ncol(P))
        for (i in c(1:nrow(P))) {
            L = L + P[i, ] * observation[i, ]
        }
    }
    return(L)
}
