#' Get the reverse complementary sequence
#'
#' Given a mutational context, get the context on its complementary
#' strand
#' 
#' @param mut_context A string of ref nucleotide context
#' @param alt_nuc The alternative nucleotide
#' 
#' @return Returns a string of nucleotide context of the same length
reverse_context <- function(mut_context, alt_nuc) {
	flipnuc <- function(nuc) {
		if (nuc == "A"){ return("T") }
		if (nuc == "T"){ return("A") }
		if (nuc == "C"){ return("G") }
		if (nuc == "G"){ return("C") }
		print(nuc)
		stop("Unknown nucleotide!")
	}
	
	if (nchar(mut_context) %% 2 != 1) {
		stop("The length of context need to be odd!")
	}
	if (nchar(alt_nuc) != 1) {
		stop("The length of alternative nucleotide need to be one")
	}

	return_context <- NULL
	for (i in seq(nchar(mut_context), by = -1)) {
		return_context <- paste(return_context, flipnuc(unlist(strsplit(
                                    mut_context, split=""))[i]), sep = "")
		if (i == nchar(mut_context) %/% 2 + 1) {
			return_context <- paste(return_context, flipnuc(alt_nuc), sep = ">")
		}
	}
	return(return_context)
}

#' Construct the spectrum from a mutation file with context
#'
#' The function converts mutation with reference contexts 
#' into a spectrum, which has mutation counts catalogued by 
#' contexts. It uses the fact that DNA has two complementary strands
#' and folds the A/G mutations into to T/C. 
#' 
#' @param input_muts A matrix of three columns; the first one is the context, 
#' with the central nucleotide mutated. The second one has the alternative
#' nucleotide. The third column stores the names of the samples. If it only has
#' two columns, the function assumes all mutations are from one single sample
#' 
#' @return Returns a matrix of mutation counts in different context, 
#' in each of the samples 
#'
#' @export

context2spec <- function(input_muts) {
	input_muts <- apply(data.frame(input_muts), c(1, 2), as.character)
	
	if (ncol(input_muts) == 2){
		print("No sample information! Assume all mutations from one sample")
	}
	if (nrow(input_muts) == 0){
		stop("The input is empty")
	}
	ref_length <- nchar(input_muts[1,1])
	if (ref_length %% 2 != 1) {
		stop("The length of ref context need to be odd!")
	}
	
	nucleotides <- c("A", "C", "G", "T")
	all_context <-  c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G") # initialize with six mutations
	for (i in seq(ref_length %/% 2)){
		all_context<-as.vector(sapply(sapply(all_context, function(x) sapply(
                        nucleotides, function(y) paste(y, x, sep=""))),
                        function(x) sapply(nucleotides, 
                        function(y) paste(x, y, sep=""))))
	}
	
	return_mat <- matrix(0, length(all_context), length(unique(input_muts[,3])))
	colnames(return_mat) <- unique(input_muts[,3])
	rownames(return_mat) <- all_context
	
	for (k in seq(nrow(input_muts))){
		if (nchar(input_muts[k, 1]) != ref_length){
			print("Warning: the ref context length is inconsistent: %s! 
                    Skipped", input_muts[k, 1])
			next
		}
		
		ref_context <- unlist(strsplit(input_muts[k, 1], ""))
		 
		if (ref_context[(ref_length %/% 2 + 1)] %in% c("A", "G")){
			row_idx <- which(all_context == reverse_context(input_muts[k, 1], 
                                input_muts[k, 2]))
			if (length(row_idx) == 0){
				print(reverse_context(input_muts[k, 1], input_muts[k, 2]))
				stop("Unexpected mutation context!")
			}
		}
		else {
			mut_context <- paste(c(ref_context[1:(ref_length %/% 2 + 1)], ">", 
                                input_muts[k, 2]), collapse="")
			mut_context <- paste(mut_context, ref_context[(ref_length %/% 2 + 
                                    1): ref_length], sep="")
			row_idx <- which(all_context == mut_context)
		}
		col_idx <- which(unique(input_muts[,3]) == input_muts[k,3])
		return_mat[row_idx, col_idx] <- return_mat[row_idx, col_idx] + 1
	}
	return(return_mat)
}

#' Construct the spectrum from a vcf file
#'
#' The function converts a vcf file into a spectrum matrix which has
#' mutation counts catalogued by contexts. It uses the fact that DNA has 
#' two complementary strands and folds the A/G mutations into to T/C. 
#' 
#' @param vcf_file A vcf file
#' @param ref_genome An fa file of reference genome sequence
#' @param context_length The length of context used. By default it is 3.
#' 
#' @return Returns a matrix of mutation counts in different context, 
#' in each of the samples 
#'
#' @export

vcf2spec <- function(vcf_file, ref_genome, context_length = 3){
	
	
	
}
