siglasso <-
function (sample_spectrum, signature, conf = 1, prior, adaptive = T, 
    gamma = 1, alpha_min = 400, iter_max = Inf, sd_multiplier = 0.5, 
    elastic_net = F) 
{
    if (missing(sample_spectrum)) {
        stop("siglasso(spectrum, signature, prior, adaptive=T, gamma=1, \n                    alpha_min=400, iter_max=20, sd_multiplier=0.5")
    }
    if (missing(signature)) {
        print("No signature supplied, will use the COSMIC signatures")
        signature = load_sig()
    }
    if (nrow(signature) != length(sample_spectrum)) {
        stop("The number of rows of signatures do not equal the number of \n            rows of the sample spectrum")
    }
    if (length(unique(sample_spectrum)) == 1) {
        if (unique(sample_spectrum) == 0) {
            stop("The spectrum is an empty vector")
        }
        random_idx <- ceiling(length(sample_spectrum) * runif(2))
        sample_spectrum[random_idx[1]] <- sample_spectrum[random_idx[1]] - 
            1
        sample_spectrum[random_idx[2]] <- sample_spectrum[random_idx[2]] + 
            1
    }
    if (missing(prior)) {
        prior = rep(1, ncol(signature))
    }
    else if (length(prior) != ncol(signature)) {
        stop("The length of prior does not equal the number of signatures!")
    }
    return(siglasso_internal(sample_spectrum, signature, prior, 
        adaptive, elastic_net, gamma, alpha_min, iter_max, sd_multiplier, 
        conf))
}
