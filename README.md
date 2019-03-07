![Image of siglasso](https://raw.githubusercontent.com/gersteinlab/siglasso/master/images/siglasso_icon.png)
## SigLASSO: Optimizing Cancer Mutation Signatures Jointly with Sampling Likelihood
![Image of siglasso](https://raw.githubusercontent.com/gersteinlab/siglasso/master/images/siglasso_schematics.png)

## Why you should use sigLASSO?
First of all, it is because mutation counts matter.
Most other methods, including deconstructSigs, is invariant to different mutation counts. Let's say we observed 10,000 mutations in a cancer sample, the number is large therefore we should have a pretty good estimation of the mutational spectrum of the tumor. Here is a synthetic example.

| Context  | Counts | Percentage |
| -------- | --- | --- |
|  AC>AA   |   500  |   5  |
|  AC>AC   |   300  |   3  |
|  AC>AG   |   100  |   1  |
|  AC>AT   |   1000 |   10 |
|   ...    |   ...  |  ... |
|  TT>GT   |   0    |   0  |

But life is not always so good. This time we sequence another tumor, but observed only 100 mutations. It could be because of exome sequencing, sequencing depth is shallow or the tumor is "silent". Now we are highly unsure about the spectrum.

| Context  | Counts | Percentage |
| -------- | --- | --- |
|  AC>AA   |   5  |   5  |
|  AC>AC   |   3  |   3  |
|  AC>AG   |   1  |   1  |
|  AC>AT   |   1  |   10 |
|   ...    |  ... |  ... |
|  TT>GT   |   0  |   0  |

Although the percentage is exactly the same as the previous sample, we are much less certain about out estimation here. For example, in the first sample, with 0/10,000 TT>GT, we are pretty certain that this mutation context is very rare. However, in the second sample, the last zero is very likely due to undersampling rather than a real, very low mutation frequency. Therefore, we should expect these two samples to have different signature solutions. A more confident one for the first one, a less, coarse one for the second sampel to reflect the uncertainty in sampling. Most methods, will give identical solutions in these two situations.

SigLASSO considers both sampling error(especially significant when the mutation count is low) and signature fitting. 

Moreover, it parameterizes the model empirically. Let the data speak for itself. Moreover, you will be able to feed prior knowledge of the signatures into the model in a soft thresholding way. No more picking up signature subsets by hand! SigLASSO achieves signature selection by using L1 regularization.


## Dependencies
To fetch the package from GitHub, you will need "devtools"
```
install.packages("devtools")
library("devtools")
```

If start with vcf files, you will additionally need bedtools and a reference genome FASTA file. See the below.

## Install
Just one line and vollà!
```
devtools::install_github("gersteinlab/siglasso")
```

## Usage
It is as simple as 2 (or 3) steps! 

### 0. Starting with VCF
If started with VCF, that is you do not have flanking nucleotide context...). To make the pacakge transparent and lightweighted, we decide to wrap around a bash script to parse mutational context. To do this, you will need bedtools(https://bedtools.readthedocs.io). For Mac OS, I find pacakge managers, like homebrew(https://brew.sh/), greatly simplify the installation.
```
brew install bedtools
```
Then you also need to have a uncompressed reference genome that is compatible to your vcf ready. 

Now, before running, one last thing is you need to prepare a meta file specify the paths of all your vcf files. An optional 2nd column can be used to specify sample names. Otherwise, it will grab the filenames as sample names.
```
/Users/data/breast/brca_cancer_1.vcf, brca_a_1
/Users/data/breast/brca_cancer_1.vcf, brca_b_1
/Users/data/prostate/prca_cancer_1.vcf, prca_c_1
```
Run

```
my_spectrum <- vcf2spec(bedtools_path, vcf_meta, ref_genome, output_file, ...)
```
This function will write "output_file" to disk. It is a context file.

### 1. Starting with a context file
Alternatively, you can start with a context file. It is easier when you already know the flanking region context. Many annotation tools now give you this. 

The context file has three columns. The 1st column is the reference nucleotide context. Typically, it is a trinucleotide. The 2nd column stores the alternative alleles. The last column has the sample names (can be unordered).

```
ATT    C    brca_a_1
CGT    C    brca_a_2
TCC    A    brca_a_1
```

Now convert this raw mutation context file into spectrum. 

```
my_context_file <- read.table(path_to_context_file)
my_spectrum <- context2spec(my_context_file)
```
It will make plots of the spectrum, to know more, see ```plot_spectrum()```
![Image of plot_spectrum](https://raw.githubusercontent.com/gersteinlab/siglasso/master/images/spec.jpeg)

### 2. Apply siglasso to context
This step is straightforward. You can supply your own signature file, or it will use the COSMIC signature. 

```
my_sigs <- siglasso(my_spectrum, ...)
```
A useful thing is prior, you can pass a vector of the length of the number of signatures, with numbers between 0 (strong preference) and 1 (no preference). We supply a prior file that we curated from COSMIC. 

```
load(cosmic_priors)
colnames(cosmic_prior)
my_prior = ifelse(cosmic_prior$BRCA==0, 0.1, 1) #adjust the strength to 0.1, BRCA
my_sigs <- siglasso(my_spectrum, prior = my_prior...)
```

### 3. Visualization
By default, siglasso() automatically generates a barplot of every samples
```
plot_sigs(my_sigs, ...)
```
There is another option to generate a dotcharts to compare between samples or samples groups

```
plot_sigs_grouped(my_sigs, [groups]...)
```
![Image of siglasso](https://raw.githubusercontent.com/gersteinlab/siglasso/master/images/plots.jpeg)

## Longer Introduction
Multiple mutational processes drive carcinogenesis, leaving characteristic signatures on tumor genomes. Determining the active signatures from the full repertoire of potential ones can help elucidate mechanisms underlying cancer initiation and development. This task involves decomposing the counts of cancer mutations, tabulated according to their trinucleotide context, into a linear combination of known mutational signatures. We formulate it as an optimization problem and develop sigLASSO, a software tool, to carry it out efficiently. SigLASSO features four key aspects: (1) By explicitly adding multinomial sampling into the overall objective function, it jointly optimizes the likelihood of sampling and signature fitting. Considering multinomial sampling is particularly important when mutation counts are low and sampling variance is high, such as in exome sequencing. (2) sigLASSO uses L1 regularization to parsimoniously assign signatures to mutation profiles, leading to sparse and more biologically interpretable solutions resembling previously well-characterized results. (3) sigLASSO fine-tunes model complexity, informed by the scale of the data and biological-knowledge based priors. In particular, instead of hard thresholding and choosing a priori a discrete subset of active signatures, sigLASSO allows continuous priors, which can be effectively learned from auxiliary information. (4) Because of this, sigLASSO can assess model uncertainty and abstain from making certain assignments in low-confidence contexts. Finally, to evaluate SigLASSO signature assignments in comparison to other approaches, we develop a set of reasonable expectations (e.g. sparsity, the ability to abstain, and robustness to noise) that we apply consistently in a variety of contexts.

## Citation 
(coming soon...)
