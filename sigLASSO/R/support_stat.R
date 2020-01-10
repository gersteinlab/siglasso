support_stat <-
function (weight_hat, weight_true) 
{
    N = ncol(weight_hat)
    F_neg = rep(NA, N)
    F_poz = rep(NA, N)
    T_poz = rep(NA, N)
    for (i in 1:N) {
        inferred_set <- as.numeric(which(weight_hat[, i] > 0))
        true_set <- as.numeric(which(weight_true[, i] > 0))
        F_neg[i] <- length(which(!(true_set %in% inferred_set)))
        F_poz[i] <- length(which(!(inferred_set %in% true_set)))
        T_poz[i] <- length(which(inferred_set %in% true_set))
    }
    return(cbind(T_poz, F_poz, F_neg))
}
