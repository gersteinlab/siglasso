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
#' @param plot a binary indicator whethere a plot of the specturm will be made 
#' (default: TRUE)
#' @param normalize Normalization choice ("none", "genome", "exome"). In general,
#' when fitting COSMIC signatures, which are derived from a mixture of WGS and WES 
#' smaples, we do not think there is additional benefits with normalization (default: none)
#' Normalization uses the trinucleotides frequencies in the genome.
#' Please do not use normalized frequencies for futher fitting (siglasso()) as it 
#' the mutation counts are no longer discrete. siglass() has an option for such 
#' normalization, operating on the signatures.
#' 
#' @return Returns a matrix of mutation counts in different context, 
#' in each of the samples 
#'
#' @export

context2spec <- function(input_muts, plot = TRUE, normalize = "none") {
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
	# initialize with six mutations
	all_context <-  c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G") 

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
                                    2): ref_length], sep="")
			row_idx <- which(all_context == mut_context)
		}
		col_idx <- which(unique(input_muts[,3]) == input_muts[k,3])
		return_mat[row_idx, col_idx] <- return_mat[row_idx, col_idx] + 1	
	}
	if (plot) {
		plot_spectrum(return_mat)
	}
	if (normalize == "genome") {
		print("Will normalize using genome trinucleotide frequency")
		data(background_context)
		return_mat<-return_mat/background_context[,1]
	}
	else if (normalize == "exome") {
		print("Will normalize using genome trinucleotide frequency")
		data(background_context)
		return_mat<-return_mat/background_context[,2]
	}
	else if (normalize != "none") {
		stop("Unexpected normalization! Use 'genome', 'exome' or 'none'")
	}
	return(return_mat)
}

#' Construct the spectrum from a VCF file
#'
#' The function converts a vcf file into a spectrum matrix which has
#' mutation counts catalogued by contexts. It uses the fact that DNA has 
#' two complementary strands and folds the A/G mutations into to T/C. 
#'  
#' It is a wrapper on bedtools(https://github.com/arq5x/bedtools2/releases). 
#'
#' @param bedtools_path Path of bedtool executable, e.g. "~/bedtools/bin/bedtool"
#' @param vcf_meta A comma(CSV)/space/tab(TAB)-delimited meta vcf file in the format of "vcf_path[comma, space or tab]sample_name",
#' Alternatively, you can just supply a list of paths to vcfs without sample_name.
#' It will then use vcf file names as sample name.:
#' @param ref_genome An fa file path of reference genome sequence
#' @param context_length The length of extension of each side of the mutation. 
#' By default it is 3 (trinucleotide context).
#' @param output_file An output file that contains the context matrix
#' @param overwrite Should it overwrite the output_file if existed
#' 
#' @return Returns a matrix of mutation counts in different context, 
#' in each of the samples 
#'
#' @export


vcf2spec <- function(bedtools_path = "bedtools", vcf_meta, ref_genome, output_file, context_length = 1, overwrite = F){
	if (!file.exists(vcf_meta)){
		stop("You need to supply a meta file specify the vcf file locations!")
	}
	if (!file.exists(ref_genome)){
		stop("You need to supply a reference genome (FASTA)!")
	}
	if (missing(output_file)){
		stop("You need to supply an output file where it write the context!")
	}
	if (file.exists(output_file)){
		if (!overwrite){
			stop("Output file exist! If you want to overwrite it, set overwrite = T")
		}
		file.remove(output_file)
	}
	
	system(paste(c("bash", system.file("/scripts/get_context.sh", package="siglasso"), vcf_meta, output_file, context_length, bedtools_path, ref_genome), collapse=" "))
	context_file<-read.table(output_file)
	return(context2spec(context_file))
}

#' Construct the spectrum from a MAF file
#'
#' The function converts a MAF file into a spectrum matrix which has
#' mutation counts catalogued by contexts. It uses the fact that DNA has 
#' two complementary strands and folds the A/G mutations into to T/C. 
#'  
#' It is a wrapper on bedtools(https://github.com/arq5x/bedtools2/releases). 
#'
#' @param bedtools_path Path of bedtool executable, e.g. "~/bedtools/bin/bedtool"
#' @param maf A maf file
#' @param ref_genome An fa file path of reference genome sequence
#' @param context_length The length of extension of each side of the mutation. 
#' By default it is 3 (trinucleotide context).
#' @param output_file An output file that contains the context matrix
#' 
#' @return Returns a matrix of mutation counts in different context, 
#' in each of the samples 
#'
#' @export


maf2spec <- function(bedtools_path = "bedtools", maf, ref_genome, context_length = 1, output_file){
	#command = paste(bedtools_path," flank -i ~/Downloads/hg001_hg38_1k.vcf -b", sep = "")
	#system("~/Downloads/bedtools2.27/bin/bedtools flank -i ~/Downloads/hg001_hg38_1k.vcf -b")
	
	if (!file.exists(maf)){
		stop("You need to supply a meta file specify the vcf file locations!")
	}
	if (!file.exists(ref_genome)){
		stop("You need to supply a reference genome (FASTA)!")
	}
	if (!file.exists(output_file)){
		stop("You need to supply an output file where it write the context!")
	}
	
	system(paste(c("bash", system.file("/scripts/get_context.sh", package="siglasso"), vcf_meta, output_file, context_length, bedtools_path, ref_genome), collapse=" "))
	context_file<-read.table(output_file)
	return(context2spec(context_file))
}
