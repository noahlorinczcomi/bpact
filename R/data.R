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
#' @format This dataset is contains allele frequency-weighted sums of gene-specific LD scores, and their squared versions, for 17479 genes spanning th 22 autosomes.
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
