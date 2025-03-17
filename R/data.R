#' Gene-based test statistic correlations
#'
#' Each entry in this list is a matrix of correlations between gene-based test statistics.
#'
#' @format This dataset is a list where each element is a matrix of correlations between gene-based test statistics due to LD between SNPs. Rows/columns of these matrices are named by the HGNC symbol of the gene to which the row/column corresponds. Rows/columns are ordered by base pair location of the transcription start site of each gene in hg19 coordinates.
#' @source Genes tested using the 1000 Genomes Phase 3 European reference panel
#' @examples
#' data(gent.Rho)
#' str(my_data)
"gent.Rho"

#' Gene LD scores
#'
#' MAF-weighted gene-specific LD scores
#'
#' @format This dataset contains allele frequency-weighted sums of gene-specific LD scores, and their squared versions, for 17479 genes spanning th 22 autosomes.
#' \describe{
#'   \item{chr}{Chromosome}
#'   \item{gene}{HGNC gene symbol}
#'   \item{m}{Number of SNPs whose LD scores were summed for the gene}
#'   \item{sumwldscore}{Sum of MAF-weighted LD scores for the gene}
#'   \item{sumwld2score}{Sum of MAF-weighted squared LD scores for the gene}
#' }
#' @source LD scores were calculated from the 1000 Genomes Phase 3 European reference panel in +-1Mb windows around gene bodies in hg19 coordinates
#' @examples
#' data(ld.df)
#' head(ld.df)
"ld.df"

#' Gene-based results for Alzheimer's disease
#'
#' Gene-based association test results for Alzheimer's disease
#'
#' @format This dataset contains results from genome-wide gene-based association testing applied to Alzheimer's disease
#' \describe{
#'   \item{gene}{HGNC gene symbol}
#'   \item{chr}{Chromosome}
#'   \item{mid}{Midpoint of gene body in hg19 coordinates}
#'   \item{pval}{Gene-based association test P-value from GenT}
#'   \item{m}{Number of SNPs used in gene-based association testing}
#'   \item{gent_sigma2_h0}{Variance of the gene-based association test statistic under the null hypothesis, which is 2 times the trace of the squared LD matrix of SNPs used}
#'   \item{gent_test_statistic}{Gent-based association test statistic from GenT, i.e. the sum of SNP chi-square statistics}
#' }
#' @source The Alzheimer's disease GWAS was from Bellenguez et al., 2021 of 480K subjects and the LD reference panel was the 1000 Genomes Phase 3 European panel
#' @examples
#' data(ad.gent)
#' head(ad.gent)
"ad.gent"

#' Gene-based results for Lewy body dementia
#'
#' Gene-based association test results for Lewy body dementia
#'
#' @format This dataset contains results from genome-wide gene-based association testing applied to Lewy body dementia
#' \describe{
#'   \item{gene}{HGNC gene symbol}
#'   \item{chr}{Chromosome}
#'   \item{mid}{Midpoint of gene body in hg19 coordinates}
#'   \item{pval}{Gene-based association test P-value from GenT}
#'   \item{m}{Number of SNPs used in gene-based association testing}
#'   \item{gent_sigma2_h0}{Variance of the gene-based association test statistic under the null hypothesis, which is 2 times the trace of the squared LD matrix of SNPs used}
#'   \item{gent_test_statistic}{Gent-based association test statistic from GenT, i.e. the sum of SNP chi-square statistics}
#' }
#' @source The Lewy body dementia GWAS was from Chia et al., 2021 of 16516 subjects and the LD reference panel was the 1000 Genomes Phase 3 European panel
#' @examples
#' data(lbd.gent)
#' head(lbd.gent)
"lbd.gent"
