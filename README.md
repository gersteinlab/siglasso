# SigLASSO: a LASSO approach jointly optimizing sampling likelihood and cancer mutation signatures

Multiple mutational processes drive carcinogenesis, leaving characteristic signatures on tumor genomes. Determining the active signatures from the full repertoire of potential ones can help elucidate the mechanisms underlying cancer initiation and development. This involves decomposing the frequency of cancer mutations categorized according to their trinucleotide context into a linear combination of known mutational signatures. We formulate this task as an optimization problem with L1 regularization and develop a software tool, sigLASSO, to carry it out efficiently. First, by explicitly adding multinomial sampling into the overall objective function, we jointly optimize the likelihood of sampling and signature fitting. This is especially important when mutation counts are low and sampling variance, high, such as the case in whole exome sequencing. sigLASSO uses L1 regularization to parsimoniously assign signatures to mutation profiles, leading to sparse and more biologically interpretable solutions. Additionally, instead of hard thresholding and choosing a priori, a discrete subset of active signatures, sigLASSO fine-tunes model complexity parameters, informed by the scale of the data and prior knowledge. Finally, it is challenging to evaluate sigLASSO’s signature assignments. To do this, we construct a set of criteria, which we can apply consistently across assignments. 


## Prerequisite
R packages "nnls" and "glmnet"

```
install.packages("nnls")
install.packages("glmnet")
```

## Usage
```
source("siglasso_code.R")
siglasso(spectrum, signature, prior, adaptive=T, gamma=1, alpha_min=400, iter_max=20, sd_multiplier=0.5)


### OUTPUTS
siglasso returns a vector of weights assigned to signatures

### INPUTS/PARAMETERS (only spectrum and signature are required)
spectrum: a vector of mutation counts of different nucleotide contexts in a single tumor sample
signature: a matrix of probablibilies on different nucleotide contexts of different signatures. 
            By default it uses 30 COSMIC mutation signatures
prior: a vector of weights on different signatures. By default it is a vector of 1s (flat prior).
adaptive: boolean variable indicating whether adaptive lasso should be used
gamma: the gamma paramether in adapative lasso. We find it is not sensitive to the analysis
alpha_min: the minimal alpha value. 
iter_max: maximum iterations. We find even in very low mutation (~20) scenarios, 
			sigLASSO convergs in about ten iterations.
sd_multiplier: cut-off for lambda tuning; how many sd should it be away from min MSE.
```
