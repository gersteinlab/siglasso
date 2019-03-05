#' Attributing signatures using sigLASSO
#' 
#' Given a mutation spectrum of one or more sample, 
#' assign signatures with weights.
#' 
#' @param sample_spectrum A data.frame/matrix of mutation counts 
#' catelogued by context
#' @param signature A signature matrix, in the form of possibilities 
#' of mutations. If not supplied, it loads 30 COSMIC signatures
#' @param conf Confidence factor of the mutational fitting (default: 0.1)
#' @param prior A vector of vector of weights on different signatures. 
#' (default: a vector of 1s (flat prior)).
#' @param adaptive A boolean variable indicating whether adaptive lasso should 
#' be used
#' @param gamma The hyperparameter gamma in adapative lasso. We find it is not 
#' sensitive to the analysis (default: 1)
#' @param alpha_min The  minimal alpha value (default:400)
#' @param iter_max Maximum iterations (default: Inf)
#' @param sd_multiplier cut-off for lambda tuning; how many sd it tolerates 
#' from min MSE.
#' @param elastic_net oolean variable indicating whether elastic net should 
#' be used
#' @param makeplot Boolean variable; if it should barplot the results
#' @return A data.frame of weights of all signatures
#' @export

siglasso <- function(sample_spectrum, signature, conf = 0.1, prior, 
						adaptive = T, gamma = 1, alpha_min = 400, 
						iter_max = Inf, sd_multiplier = 0.5, elastic_net = F, 
						makeplot = T) {
    if (missing(sample_spectrum)) {
        stop("siglasso(spectrum, signature, prior, adaptive=T, gamma=1, 
                    alpha_min=400, iter_max=20, sd_multiplier=0.5")
    }
    if (missing(signature)) {
        print("No signature supplied, will use the COSMIC signatures")
        signature <- devtools::use_data(cosmic30sig)
    }
    if (nrow(signature) != length(sample_spectrum)) {
        stop("The number of rows of signatures do not equal the number of 
            rows of the sample spectrum")
    }
    
    if (missing(prior)) {
        prior <- rep(1, ncol(signature))
    } else if (length(prior) != ncol(signature)) {
        stop("The length of prior does not equal the number of signatures!")
    }

	sample_spectrum<-data.frame(sample_spectrum)
   
	if (ncol(sample_spectrum < 2)) {
		if (length(unique(sample_spectrum)) == 1) {
			if (unique(sample_spectrum) == 0) {
				stop("The spectrum is an empty vector")
		    }
			random_idx <- ceiling(length(sample_spectrum) * runif(2))
			sample_spectrum[random_idx[1]] <- sample_spectrum[random_idx[1]] - 1
			sample_spectrum[random_idx[2]] <- sample_spectrum[random_idx[2]] + 1
		}
		return_weights<-data.frame(siglasso_internal(sample_spectrum, 
														signature, prior, 
														adaptive, elastic_net, 
														gamma, alpha_min, 
														iter_max, sd_multiplier,
														conf))
	}
	else {
		sprintf("There are %d samples", ncol(sample_spectrum))
		return_weights<-NULL
		for (k in ncol(sample_spectrum)) {
			if (length(unique(sample_spectrum[,k])) == 1) {
				if (unique(sample_spectrum[,k]) == 0) {
					stop("The spectrum is an empty vector")
				}
				random_idx <- ceiling(length(sample_spectrum[,k]) * runif(2))
				sample_spectrum[random_idx[1]] <- sample_spectrum[
													random_idx[1],k] - 1
				sample_spectrum[random_idx[2]] <- sample_spectrum[
													random_idx[2],k] + 1
			}
			return_weights<-cbind(return_weights, siglasso_internal(
									sample_spectrum[,k], signature, prior, 
									adaptive, elastic_net, gamma, alpha_min, 
									iter_max, sd_multiplier, conf)) 
		}
	}
	if (makeplot){
		plot_sigs(return_weights)
	}
	return(return_weights)
}
