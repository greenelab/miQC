#' Basic scRNA-seq QC metrics from an ovarian tumor
#'
#' Basic QC metrics from a high-grade serous ovarian cancer (HGSOC) sample.
#' The information included is the minimum needed to generate a miQC model and
#' make plots. Count data for the full tumor is available through GEO,
#' accession number GSM4816047.  
#'
#' @format A data frame with 1000 rows and 2 variables:
#' \describe{
#'   \item{detected}{number of unique genes detected}
#'   \item{subsets_mito_percent}{% of reads mapping to mitochondrial genes}
#'   ...
#' }
#' @source GEO Accession GSM4816047
"metrics"
