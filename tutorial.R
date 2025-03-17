

### Demonstrating for two traits locally







# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
### Demonstrating locally
gent.Rho=readRDS('/Volumes/chengflab/Cheng-Noah/reference_data/gent_stat_ld/full_matrices/full_matricesEUR.Rds')
ld.df=readRDS('/Volumes/chengflab/Cheng-Noah/reference_data/gene_ldscores/EUR/allchrs_maf_weighted.Rds')
## trait 1
gent.data1=readRDS('/Volumes/chengflab/Cheng-Noah/GenT_results/AD_EUR.Rds')
gp=fread('/Volumes/chengflab/Cheng-Noah/manuscripts/druggable_genes/data/genepositions_hg19.txt')
gent.data1=gent.data1 %>%
  left_join(gp %>% select(gene=symbol,chr,mid),by='gene') %>%
  mutate(gent_test_statistic=gent_mu_h1-gent_mu_h0,
         m=gent_mu_h0,
         pval=gent_p)
gwasn1=480000
prior.fit1=compositemh(gent.data1,ld.df,gent.Rho,gwasn1,verbose=TRUE,chain_length=200)
posteriors1=posterior_gene(gent.data1,ld.df,gwasn1,prior.fit1)
## trait 2
gent.data2=readRDS('/Volumes/chengflab/Cheng-Noah/GenT_results/LBD_EUR.Rds')
gent.data2=gent.data2 %>%
  left_join(gp %>% select(gene=symbol,chr,mid),by='gene') %>%
  mutate(gent_test_statistic=gent_mu_h1-gent_mu_h0,
         m=gent_mu_h0,
         pval=gent_p)
gwasn2=16516
prior.fit2=compositemh(gent.data2,ld.df,gent.Rho,gwasn2,verbose=TRUE,chain_length=200)
posteriors2=posterior_gene(gent.data2,ld.df,gwasn2,prior.fit2)
## shared proportions
ps=propshared(posteriors1,posteriors2,gwasn1,gwasn2,gent.Rho,thr_r2=0.5,max_size=100)
ps$count_results

# adjusting shared proportions for sample overlap
simex_res=shared_count_simex(
  M=8148, # number of jointly tested independent genes
  mcausalgenes1=ps$count_results$estimated_count[1], # estimated number of non-causal trait 1 genes
  mcausalgenes2=ps$count_results$estimated_count[2], # estimated number of non-causal trait 2 genes
  ngwas1=gwasn1, # sample size of trait 1 GWAS
  ngwas2=gwasn2, # sample size of trait 2 GWAS
  estimated_overlap_count=ps$count_results$estimated_count[3], # data-estimated number of overlapping causal genes
  upsilon_overlap=0.2, # n01/sqrt(n1*n2)*Corr(t1,t2)
  # assumed simulation parameters
  h21=0.1, # SNP heritability trait 1
  h22=0.1, # SNP heritability trait 1
  niter=100, # number of iterations
  m=50, # number of tested SNPs per gene
  mcausal_snps_per_gene=5, # number of causal SNPs per causal gene
  LD_type='AR', # type of LD correlation matrix (AR for autoregressive, CS for compound symmetry, I for independence)
  LD_rho=0.5, # correlation parameter
  nlambdas=10, # number of lambdas to evaluate in SIMEX
  doplot=TRUE, # should a plot be created?
  verbose=FALSE
)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
### Demonstrating remotely
library(data.table);library(dplyr);library(bigsnpr);
library(mvnfast,lib='/home/lorincn/Rpkgs')
# library(polygene,lib='/home/lorincn/Rpkgs')
gent.Rho=readRDS('~/isilon/Cheng-Noah/reference_data/gent_stat_ld/full_matrices/full_matricesEUR.Rds')
ld.df=readRDS('~/isilon/Cheng-Noah/reference_data/gene_ldscores/EUR/allchrs_maf_weighted.Rds')
## trait 1
gent.data1=readRDS('~/isilon/Cheng-Noah/GenT_results/AD_EUR.Rds')
gp=fread('~/isilon/Cheng-Noah/manuscripts/druggable_genes/data/genepositions_hg19.txt')
gent.data1=gent.data1 %>%
  left_join(gp %>% select(gene=symbol,chr,mid),by='gene') %>%
  mutate(gent_test_statistic=gent_mu_h1-gent_mu_h0,
         m=gent_mu_h0,
         pval=gent_p)
gwasn1=480000
t0=Sys.time()
prior.fit1=compositemh(gent.data1,ld.df,gent.Rho,gwasn1,verbose=FALSE,chain_length=1000)
posteriors1=posterior_gene(gent.data1,ld.df,gwasn1,prior.fit1,verbose=FALSE)
t1=Sys.time()
d1=as.numeric(t1-t0)
cat(d1,'\n')
## trait 2
gent.data2=readRDS('~/isilon/Cheng-Noah/GenT_results/LBD_EUR.Rds')
gent.data2=gent.data2 %>%
  left_join(gp %>% select(gene=symbol,chr,mid),by='gene') %>%
  mutate(gent_test_statistic=gent_mu_h1-gent_mu_h0,
         m=gent_mu_h0,
         pval=gent_p)
gwasn2=16516
t0=Sys.time()
prior.fit2=compositemh(gent.data2,ld.df,gent.Rho,gwasn2,verbose=FALSE,chain_length=1000)
posteriors2=posterior_gene(gent.data2,ld.df,gwasn2,prior.fit2,verbose=FALSE)
t1=Sys.time()
d2=d1+as.numeric(t1-t0)
cat(d2,'\n')
## shared proportions
t0=Sys.time()
ps=propshared(posteriors1,posteriors2,gwasn1,gwasn2,gent.Rho,thr_r2=0.5,max_size=100,verbose=FALSE)
t1=Sys.time()
d3=d2+as.numeric(t1-t0)
cat(d3,'\n')
ps$count_results

# adjusting shared proportions for sample overlap
t0=Sys.time()
simex_res=shared_count_simex(
  M=8148, # number of jointly tested independent genes
  mcausalgenes1=ps$count_results$estimated_count[1], # estimated number of non-causal trait 1 genes
  mcausalgenes2=ps$count_results$estimated_count[2], # estimated number of non-causal trait 2 genes
  ngwas1=gwasn1, # sample size of trait 1 GWAS
  ngwas2=gwasn2, # sample size of trait 2 GWAS
  estimated_overlap_count=ps$count_results$estimated_count[3], # data-estimated number of overlapping causal genes
  upsilon_overlap=0.2, # n01/sqrt(n1*n2)*Corr(t1,t2)
  # assumed simulation parameters
  h21=0.1, # SNP heritability trait 1
  h22=0.1, # SNP heritability trait 1
  niter=100, # number of iterations
  m=50, # number of tested SNPs per gene
  mcausal_snps_per_gene=5, # number of causal SNPs per causal gene
  LD_type='AR', # type of LD correlation matrix (AR for autoregressive, CS for compound symmetry, I for independence)
  LD_rho=0.5, # correlation parameter
  nlambdas=10, # number of lambdas to evaluate in SIMEX
  doplot=TRUE, # should a plot be created?
  verbose=FALSE
)
t1=Sys.time()
d4=d3+as.numeric(t1-t0)
cat(d4,'\n')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
