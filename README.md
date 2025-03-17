## Overview
We calculate posterior probabilites *PRP<sub>kt</sub>* that the latent causal indicator *I<sub>t</sub>* for the *t*th trait and *k*th gene takes value 1 using the gene-based test statistic $T_{kt}$, i.e. that the *k*th gene is causal for the *th* trait under the assumed causal model:

<p>
$$
PRP_{kt}:=P(I_{kt}=1|T_{kt})=\int_{\delta_t\in \Delta_t}\int_{\tau_t\in T_t} P(I_{kt}|T_{kt},\delta_t,\tau_t)p_{\delta_t}(\delta_t|T_{kt})p_{\tau_t}(\tau_t|T_{kt})d\delta_t d\tau_t .
$$
</p>

The parameters $\delta_t$ and $\tau_t$ define the causal model, and we account for imprecise estimation of their values by integrating over their prior distributions $p_{\tau_t}$ and $p_{\delta_t}$. These prior distributions are parameterized empirically using the Metropolis-Hastings algorithm. The parameter $\delta_t$ is the proportion of non-causal genes for the *t*th trait, and $\tau_t$ is the variance- and LD score-scaled SNP heritability for the *t*th trait contributed by each SNP. Since genes spanning the genome are not generally independent of each other, we apply an iterative composite likelihood approach in which independent sets of genes are randomly selected, in which the model is fitted and across which inference is made.

## Preprint
Please see our preprint for more details about the BPACT method and the results of applying it to 32 complex traits.

## Implementation
Below is an example of how to apply **B**ayesian estimation of **p**olygenic **a**rchitectures for **c**omplex **t**raits for Alzheimer's disease and Lewy body dementia (as an example).
```R
library(bpact)
# load necessary data
data(gent.Rho)
data(ld.df)

# load example data
data(ad.gent)
data(lbd.gent)

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
  estimated_overlap_count=ps$count_results$estimated_count[3], # data-estimated number of overlapping causal genes
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
