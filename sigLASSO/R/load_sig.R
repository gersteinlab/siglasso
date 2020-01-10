load_sig <-
function () 
{
    signature <- read.table("./data/signatures_probabilities.txt", 
        sep = "\t", header = T)
    sig <- signature[, 3:33]
    sig$Somatic.Mutation.Type <- gsub("\\[", "", sig$Somatic.Mutation.Type)
    sig$Somatic.Mutation.Type <- gsub("]", "", sig$Somatic.Mutation.Type)
    sig$Somatic.Mutation.Type <- gsub(">", "", sig$Somatic.Mutation.Type)
    sig <- sig[order(substring(sig$Somatic.Mutation.Type, 4, 
        4)), ]
    sig <- sig[order(substring(sig$Somatic.Mutation.Type, 1, 
        1)), ]
    sig <- sig[order(substring(sig$Somatic.Mutation.Type, 2, 
        3)), ]
    return(data.matrix(sig[, -1]))
}
