# Overview
We calculate the posterior probability $PRP_{kt}$ that the latent causal indicator $I_{kt}$ for the *t*th trait and *k*th gene takes value 1 using the gene-based test statistic $T_{kt}$, i.e. that the *k*th gene is causal for the t*th* trait under the assumed causal model:

<p>
$$
PRP_{kt}:=P(I_{kt}=1|T_{kt})=\int_{\delta_t\in \Delta_t}\int_{\tau_t\in T_t} P(I_{kt}|T_{kt},\delta_t,\tau_t)p_{\delta_t}(\delta_t|T_{kt})p_{\tau_t}(\tau_t|T_{kt})d\delta_t d\tau_t .
$$
</p>

The parameters $\delta_t$ and $\tau_t$ define the causal model, and we account for imprecise estimation of their values by integrating over their prior distributions $p_{\tau_t}$ and $p_{\delta_t}$. These prior distributions are parameterized empirically using the Metropolis-Hastings algorithm. The parameter $\delta_t$ is the proportion of non-causal genes for the *t*th trait, and $\tau_t$ is the variance- and LD score-scaled SNP heritability for the *t*th trait contributed by each SNP. Since genes spanning the genome are not generally independent of each other, we apply an iterative composite likelihood approach and randomly select genes within independent blocks. We fit the model in each block then pool results across iterations to make the inference.

Please see our [preprint](https://www.medrxiv.org/content/10.1101/2025.03.17.25324106v1) for more details about the BPACT method and the results of applying it to 32 complex traits.

Lorincz-Comi, N. & Cheng, F. Bayesian estimation of shared polygenicity identifies drug targets and repurposable medicines for human complex diseases, *medRxiv*, doi: https://doi.org/10.1101/2025.03.17.25324106.

# Implementation
Below is an example of how to apply **B**ayesian estimation of **p**olygenic **a**rchitectures for **c**omplex **t**raits (BPACT) for Alzheimer's disease and Lewy body dementia as an example.

Here is a flowchart showing where the data that BPACT requires comes from:
<div align="center">
  <img src="https://github.com/noahlorinczcomi/bpact/blob/main/flowchart.png" width="1000"/>
</div>

## Formatting data
```bpact``` requires genome-wide gene-based association test statistics to be in a data frame which has the following columns present:
  - ```gene```: Any gene-specific identifier (e.g., HGNC symbol)
  - ```chr```: Chromosome of the gene
  - ```mid```: Midpoint of gene body, but could also be TSS, start, or end of gene body; just needs to be a gene-specific chromosomal location
  - ```pval```: P-value used to test the null hypothesis that no SNPs in the gene-specific set are associated with the phenotype, e.g. from the test in [GenT](https://github.com/noahlorinczcomi/gent)
  - ```m```: The number of SNPs in the gene-specific set
  - ```gent_sigma2_h0```: The variance of the gene-based test statistic under the null hypothesis; i.e., $$2\text{trace}($$**RR**$$)$$, where **R** is the LD matrix for SNPs in the gene-specific set
  - ```gent_test_statistic```: The sum of SNP chi-squares for SNPs in the gene-specific set, i.e., $$\sum_{i=1}^mZ_i^2$$ where $Z_i$ is the Z-statistic for the *i*th SNP

These columns are automatically outputted by gene-based association test methods such as [GenT](https://github.com/noahlorinczcomi/gent).

Here is what the example data ```ad.gent``` and ```lbd.gent```, included with the ```bpact``` R package, looks like:
```R
remotes::install_github('noahlorinczcomi/bpact')
library(bpact)
library(dplyr)
library(Matrix)

# load necessary data
data(gent.Rho)
data(ld.df)

# load example data
data(ad.gent)
data(lbd.gent)

# View format:
str(ad.gent)
'data.frame':	17166 obs. of  7 variables:
 $ gene               : chr  "A1BG" "A1CF" "A3GALT2" "A4GALT" ...
 $ chr                : chr  "19" "10" "1" "22" ...
 $ mid                : num  5.89e+07 5.26e+07 3.38e+07 4.31e+07 1.38e+08 ...
 $ pval               : num  0.536 0.392 0.337 0.941 0.712 ...
 $ m                  : int  277 263 269 630 221 203 497 448 573 330 ...
 $ gent_sigma2_h0     : num  17064 16817 18546 65351 28789 ...
 $ gent_test_statistic: num  246 277 305 288 109 ...
 - attr(*, "na.action")= 'omit' Named int [1:934] 3 4 48 62 140 162 179 186 219 356 ...
  ..- attr(*, "names")= chr [1:934] "3" "4" "48" "62" ...

str(lbd.gent)
'data.frame':	16954 obs. of  7 variables:
 $ gene               : chr  "A1BG" "A1CF" "A3GALT2" "A4GALT" ...
 $ chr                : chr  "19" "10" "1" "22" ...
 $ mid                : num  5.89e+07 5.26e+07 3.38e+07 4.31e+07 1.38e+08 ...
 $ pval               : num  0.616 0.86 0.376 0.515 0.568 ...
 $ m                  : int  190 247 262 593 208 195 370 428 543 311 ...
 $ gent_sigma2_h0     : num  6434 16519 17816 60467 27842 ...
 $ gent_test_statistic: num  157 118 282 551 142 ...
 - attr(*, "na.action")= 'omit' Named int [1:1146] 3 4 48 62 140 162 179 186 219 356 ...
  ..- attr(*, "names")= chr [1:1146] "3" "4" "48" "62" ...
```

## Running the analysis
Below we estimate priors and fit posterior models to (i) calculate posterior risk probabilities (PRPs) that each gene is causal under a certain heritability model and (ii) estimate the number of genes which are shared between Alzheimer's disease (AD) and Lewy body dementia (LBD).
```R
# Alzheimer's disease (AD) posterior probabilities
ad_priors=compositemh(ad.gent, ld.df, gent.Rho, 480000, verbose=TRUE, chain_length=1000)
ad_posteriors=posterior_gene(ad.gent, ld.df, 480000, ad_priors)

# Lewy body dementia (LBD) posterior probabilities
lbd_priors=compositemh(lbd.gent, ld.df, gent.Rho, 16516, verbose=TRUE, chain_length=1000)
lbd_posteriors=posterior_gene(lbd.gent, ld.df, 16516, lbd_priors)

# Estimate number of shared associated genes between AD and LBD
shared_genes=propshared(ad_posteriors,lbd_posteriors,480000,16516,gent.Rho)

# Adjusting shared associated gene counts for GWAS sample overlap.
## NOTE: Only necessary if there is nonzero sample overlap and the traits are correlated
## Our experience is that sample overlap generally has a negligible impact on estimated shared counts,
## even for phenotypically correlated traits with large GWAS sample overlap.
adj_shared_genes=shared_count_simex(
  M=8148, # number of jointly tested independent genes
  mcausalgenes1=shared_genes$count_results$estimated_count[1], # estimated number of non-causal AD genes
  mcausalgenes2=shared_genes$count_results$estimated_count[2], # estimated number of non-causal LBD genes
  ngwas1=gwasn1, # sample size of AD GWAS
  ngwas2=gwasn2, # sample size of LBD GWAS
  estimated_overlap_count=shared_genes$count_results$estimated_count[3], # data-estimated number of overlapping causal genes
  upsilon_overlap=0.2, # ~n01/sqrt(n1*n2)*Corr(t1,t2)
  # assumed simulation parameters
  h21=0.1, # SNP heritability AD
  h22=0.1, # SNP heritability LBD
  niter=100, # number of iterations
  m=50, # number of tested SNPs per gene
  mcausal_snps_per_gene=5, # number of causal SNPs per causal gene
  LD_type='AR', # type of LD correlation matrix (AR for autoregressive, CS for compound symmetry, I for independence)
  LD_rho=0.5, # correlation parameter
  nlambdas=10, # number of lambdas to evaluate in SIMEX
  doplot=TRUE, # should a plot be created?
  verbose=FALSE
)
```

Below is the image produced by the above SIMEX adjustment, which shows that GWAS sample overlap has no impact on the estimated number of shared genes between Alzheimer's disease and Lewy body dementia.
<div align="center">
  <img src="https://github.com/noahlorinczcomi/bpact_analysis/blob/main/AD_LBD_sharing_from_sample_overlap%20copy.png" width="500"/>
</div>

## Schizophrenia and Bipolar I/II example
We also show the example below between schizophrenia and Bipolar I/II disorder to illustrate that GWAS sample overlap is not generally negligible.
<div align="center">
  <img src="https://github.com/noahlorinczcomi/bpact_analysis/blob/main/SCZ_BIP_sharing_from_sample_overlap_updatedJan24%20copy.png" width="550"/>
</div>

