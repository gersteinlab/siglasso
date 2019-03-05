# SigLASSO: Optimizing Cancer Mutation Signatures Jointly with Sampling Likelihood
![Image of siglasso](https://raw.githubusercontent.com/gersteinlab/siglasso/master/images/siglasso_schematics.png)

## Why you should use sigLASSO?
SigLASSO considers both sampling error(especially significant when the mutation count is low) and signature fitting. It parameterizes the model empirically. Let the data speak for itself. Moreover, you will be able to feed prior knowledge of the signatures into the model in a soft thresholding way. No more picking up signature subsets by hand! SigLASSO achieves signature selection by using L1 regularization.

## Introduction
Multiple mutational processes drive carcinogenesis, leaving characteristic signatures on tumor genomes. Determining the active signatures from the full repertoire of potential ones can help elucidate mechanisms underlying cancer initiation and development. This task involves decomposing the counts of cancer mutations, tabulated according to their trinucleotide context, into a linear combination of known mutational signatures. We formulate it as an optimization problem and develop sigLASSO, a software tool, to carry it out efficiently. SigLASSO features four key aspects: (1) By explicitly adding multinomial sampling into the overall objective function, it jointly optimizes the likelihood of sampling and signature fitting. Considering multinomial sampling is particularly important when mutation counts are low and sampling variance is high, such as in exome sequencing. (2) sigLASSO uses L1 regularization to parsimoniously assign signatures to mutation profiles, leading to sparse and more biologically interpretable solutions resembling previously well-characterized results. (3) sigLASSO fine-tunes model complexity, informed by the scale of the data and biological-knowledge based priors. In particular, instead of hard thresholding and choosing a priori a discrete subset of active signatures, sigLASSO allows continuous priors, which can be effectively learned from auxiliary information. (4) Because of this, sigLASSO can assess model uncertainty and abstain from making certain assignments in low-confidence contexts. Finally, to evaluate SigLASSO signature assignments in comparison to other approaches, we develop a set of reasonable expectations (e.g. sparsity, the ability to abstain, and robustness to noise) that we apply consistently in a variety of contexts.

## Dependencies
To fetch the package from GitHub (we are working on the CRAN submission!), you will need "devtools"
```
install.packages("devtools")
library("devtools")
```

## Install
Just one line and vollà!
```
devtools::install_github("gersteinlab/siglasso")
```

## Usage
```
siglasso(spectrum, signature, conf = 0.1, prior, adaptive = T, gamma =  1, alpha_min = 400, 
		iter_max = 20, sd_multiplier = 0.5, elastic_net = F, makeplot = T)


### OUTPUTS
siglasso returns a data.frame of weights assigned to signatures

### INPUTS/PARAMETERS (only spectrum and signature are required)
spectrum: a data.frame of mutation counts of different nucleotide contexts in a single tumor sample
signature: a matrix of probabilities on different nucleotide contexts of different signatures. 
            By default it uses 30 COSMIC mutation signatures
conf: confidence factor of the signature fitting
prior: a vector of weights on different signatures. By default, it is a vector of 1s (flat prior).
adaptive: boolean variable indicating whether adaptive lasso should be used
gamma: the gamma parameter in the adapative lasso. We find it is not sensitive to the analysis
alpha_min: the minimal alpha value. 
iter_max: maximum iterations. We find even in very low mutation (~20) scenarios, 
			sigLASSO converges in about ten iterations.
sd_multiplier: cut-off for lambda tuning; how many sd should it be away from min MSE.
elastic_net: should it use elastic_net instead of LASSO?
makeplot: should it make a barplot of the results?
```
