#' Visualize the output of siglasso as barplots
#'
#' Show a barplot of the signature fractions in multiple samples
#'
#' @param sig_weights A matrix of weigths of signatures in multiple 
#' samples
#' @param sample_order A vector of how samples should be ordered. If not 
#' supplied, it will order by the sum of weights 
#' 
#' @export
#'
plot_sigs <- function(sig_weights, sample_order){
	par(mar = c(5.1,4.1,4.1,6.1))
	col_palette <- colorRampPalette(brewer.pal(11, "Spectral"))
	total_sigs <- length(rowSums(sig_weights)>0)
	
	if (missing(sample_order)){
		sample_order<-order(colSums(sig_weights))
	}
	if (length(sample_order) != ncol(sig_weights)){
		stop("The length of sample_order needs to match the number of columns 
                of sig_weights")
	}

	barplot(data.matrix(sig_weights[,sample_order]), xaxt = "n", col = 
            col_palette(total_sigs), xlab = "Samples", 
            ylab = "Mutation frac.", legend=NULL)
	legend(par("usr")[2], par("usr")[4], which(rowSums(sig_weights)>0), 
            col = col_palette(total_sigs), pch = 15, xpd = NA)
}

#' Visualize the output of siglasso, grouped and averaged, as a dotchart
#'
#' Show a dotchart of the signature fractions in samples, summarized
#' by groups. If groups are not supplied, it plots every sample 
#' 
#' @param sig_weights A matrix of weigths of signatures in multiple 
#' samples
#' @param groups A vector of factors indicating the groups of each
#' samples. If not supplied, it will plot all samples individually
#' @param stat How to aggregate the weights in a group. Choose "mean" or 
#' "median". Default is mean. 
#'
#' @export
plot_sigs_grouped <- function(sig_weights, groups, stat = "mean"){
	par(mar = c(5.1,4.1,6.1,3.1))
	
	if (missing(groups)){
		groups <- seq(1, ncol(sig_weights))
	}
	else if (length(groups) != ncol(sig_weights)){
		stop("The length of groups needs to match the number of columns 
                of sig_weights")
	}
	
	sig_weights <- t(sig_weights[rowSums(sig_weights)>0, ])
	
	col_palette <- colorRampPalette(brewer.pal(11, "Spectral"))
	if (stat == "mean") {
		plot_data <- aggregate(sig_weights, list(groups), mean)
	}
	else if (stat == "median") {
		plot_data <- aggregate(sig_weights, list(groups), median)
	}
	rownames(plot_data) <- plot_data[, 1]	
	plot_data <- data.matrix(plot_data[, -1])
	colnames(plot_data) <- seq(1:ncol(sig_weights))
	
	dotchart(plot_data, pch = 16, cex = 1.2, col = col_palette(nrow(plot_data)))
	legend(par("usr")[1], par("usr")[4] * 1.1, rownames(plot_data), 
            xpd = NA, col = col_palette(nrow(plot_data)), pch = 16, 
            ncol = nrow(plot_data))
}

#' Visualize mutational spectrums of a sample, or all samples summarized
#'
#' Show a plot of mutational spectrums of all samples. That is a plot of 
#' mutation counts catalogued by mutational contexts.
#'
#' @param sig_weights A matrix of mutational counts of different contexts
#  in multiple samples
#' @param stat How to aggregate the spectrums. Choose "mean" or 
#' "median". Default is mean. 
#' @param std96sigs Boolean variable, indicating if using 96
#' contexts and sorted by the contexts (from AC>AA to TT>GT)
#'
#' @export
plot_spectrum <- function(sig_weights, stat="mean", std96sigs = T){
	par(mar = c(5.1,4.1,6.1,3.1))	
	col_palette <- unlist(lapply(brewer.pal(6, "Spectral"), 
                                            function(x) rep(x, 16)))
	
	if (std96sigs && nrow(sig_weights) != 96){
		sprintf("The number of rows is %d, not 96!", nrow(sig_weights))
		std96sigs <- F
	}
	
	if (stat == "mean") {
		plot_data <- apply(sig_weights, 1, mean)
	}
	else if (stat == "median") {
		plot_data <- apply(sig_weights, 1, median)
	}

	if (std96sigs){
		mut <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
		plot(plot_data, type = "h", bty = 'n', xaxt = "n", 
            col = col_palette, lwd = 5, xlab = "", ylab = "Mutation counts")
		text(seq(8, 88, by = 16), par("usr")[4] * 1.1, xpd = NA, mut, 
                col = unique(col_palette), font = 2)
	}
	
	else{
		col_palette <- colorRampPalette(brewer.pal(11, "Spectral"))
		plot(plot_data, type = "h", bty = 'n', xaxt = "n", 
                col = col_palette(length(plot_data)), lwd = 5, xlab = "", 
                ylab = "Mutation counts")
	}
}


