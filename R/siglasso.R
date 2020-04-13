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
#' from min MSE (default: 1.0).
#' @param elastic_net Boolean variable indicating whether elastic net should 
#' be used
#' @param plot Boolean variable; if it should barplot the results
#' @param plot_colors You can supply your own colors here (need to equal the
#' number of signatures), or by default the program will assign colors.
#' @param plot_orders How to reorder when plotting. It can be 
#' supplied as NULL, that is no reorder, or reordering by the sum of 
#' signature assignment (see reorder) or a vector, indicating the orders 
#' (default: NULL)  
#' @param reorder Boolean variable indicating whether reorder when plotting
#' or not. (default: False)
#' @param normalize Whethere to normalize the signatures before fitting. 
#' "none", "genome" or "exome". default: "none" (also recommended for COSMIC 
#' signatures). Note here because normalized mutation counts will be 
#' not discrete anymore. We applied the normalization on the signature.
#' They are mathmatically equivalent. 
#' @param default_sig Default signature type if no signatures were supplied;
#' three choices available so far, "cosmic_v2" (default) or "cosmic_v3_exo", 
#' "cosmic_v3_wholegenome" 
#' @return A data.frame of weights of all signatures
#' @export

siglasso <- function(sample_spectrum, signature, conf = 0.1, prior, 
                        adaptive = T, gamma = 1, alpha_min = 400, 
                        iter_max = Inf, sd_multiplier = 1.0, elastic_net = F, 
                        plot = T, plot_colors = NA, plot_orders = NULL,
                        reorder = F, normalize = "none", 
                        default_sig = "cosmic_v2") {
    if (missing(sample_spectrum)) {
        stop("siglasso(spectrum, signature, prior, adaptive=T, gamma=1, 
                    alpha_min=400, iter_max=20, sd_multiplier=0.5")
    }
    if (missing(signature)) {
        if (default_sig == "cosmic_v2") {
            print("No signature supplied, will use COSMIC v2 signatures")
            data(cosmic30sig)
            signature <- data.matrix(cosmic30sig)
        }
        else if (default_sig == "cosmic_v3_exo") {
            print("No signature supplied, will use COSMIC v3 (exome)")
            data(cosmic_v3_exo)
            signature <- data.matrix(cosmic_v3_exo)
        }
        else if (default_sig == "cosmic_v3_wholegenome") {
            print("No signature supplied, will use COSMIC v3 (whole genome)")
            data(cosmic_v3_wholegenome)
            signature <- data.matrix(cosmic_v3_wgs)
        }
        else {
            stop("Unknown default_sig!")
        }
    }
    if (nrow(signature) != nrow(sample_spectrum)) {
        stop("The number of rows of signatures do not equal the number of 
            rows of the sample spectrum")
    }
    
    if (missing(prior)) {
        prior <- rep(1, ncol(signature))
    } else if (length(prior) != ncol(signature)) {
        stop("The length of prior does not equal the number of signatures!")
    }

    if (normalize == "genome") {
        print("Will normalize using genome trinucleotide frequency")
        data(background_context)
        signature<-signature*background_context[,1]
        signature<-t(t(signature)/colSums(signature)) 
    }
    else if (normalize == "exome") {
        print("Will normalize using genome trinucleotide frequency")
        data(background_context)
        signature<-signature*background_context[,2]
        signature<-t(t(signature)/colSums(signature)) 
    }
    else if (normalize != "none") {
        stop("Unexpected normalization! Use 'genome', 'exome' or 'none'")
    }

    sample_spectrum<-data.frame(sample_spectrum)
   
    print(sprintf("There are %d samples", ncol(sample_spectrum)))
    return_weights<-NULL
    for (k in seq(ncol(sample_spectrum))) {
        if (length(unique(sample_spectrum[,k])) == 1) {
            if ((unique(sample_spectrum[,k]))[1] == 0) {
                stop("The spectrum is an empty vector")
            }
            random_idx <- ceiling(length(sample_spectrum[,k]) * runif(2))
            sample_spectrum[random_idx[1]] <- sample_spectrum[
                                                random_idx[1],k] - 1
            sample_spectrum[random_idx[2]] <- sample_spectrum[
                                                random_idx[2],k] + 1
        }
        return_weights <- cbind(return_weights, siglasso_internal(
                                sample_spectrum[,k], signature, prior, 
                                adaptive, elastic_net, gamma, alpha_min, 
                                iter_max, sd_multiplier, conf)) 
    }
    colnames(return_weights) <- colnames(sample_spectrum)
    if (plot){
        if (!is.na(plot_colors)) {
            if (length(plot_colors) != ncols(signature)){
                print(paste("The number of provided plotting colors does not",
                "equal the number of signatures"))
                print("There might be errors in plotting")
            }
        }
        if (length(plot_orders)>0) {
            if (length(plot_orders) != ncol(sample_spectrum)) {
                stop(paste("The number of provided plot_orders does not",
                "equal the number of samples!"))
            }
            plot_sigs(return_weights, plot_colors, sample_order = plot_orders)
        }
        else {
            plot_sigs(return_weights, plot_colors, re_order = )
        }
    }
    return(return_weights)
}
