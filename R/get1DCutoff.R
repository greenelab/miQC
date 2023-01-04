#' get1DCutoff
#'
#' When the model is generated based only on mitochondrial percentage, e.g.
#' run_model(sce, model_type = "one_dimensional"), there will be a discrete 
#' cutoff point at which to remove cells with at least that mitochondrial
#' percentage. This function identifies this cutoff based on a given posterior
#' probability threshold. This number can then be passed as a cutoff to other
#' modalities.
#'
#' @param sce (SingleCellExperiment) Input data object.
#'
#' @param model (flexmix) Output of mixtureModel function, which should be
#'   explicitly called first to ensure stability of model parameters.
#'   Default = NULL.
#'
#' @param posterior_cutoff (numeric) The posterior probability of a cell being
#'   part of the compromised distribution, a number between 0 and 1. Any cells
#'   below the appointed cutoff will be marked to keep.
#'   Default = 0.75
#'
#' @param subsets_mito_percent (character) Column name in sce giving the
#'   percent of reads mapping to mitochondrial genes. This name is inherited
#'   from scater's addPerCellQC() function, provided the subset "mito" with
#'   names of all mitochondrial genes is passed in. See examples for details.
#'
#' @return Returns a single numeric value, the percent mitochondrial cutoff
#' at which to filter cells.
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom flexmix parameters posterior
#'
#' @export
#'
#' @examples
#' library(scRNAseq)
#' library(scater)
#' sce <- ZeiselBrainData()
#' mt_genes <- grepl("^mt-",  rownames(sce))
#' feature_ctrls <- list(mito = rownames(sce)[mt_genes])
#' sce <- addPerCellQC(sce, subsets = feature_ctrls)
#' model <- mixtureModel(sce, model_type = "one_dimensional")
#' get1DCutoff(sce, model)

get1DCutoff <- function(sce, model = NULL, posterior_cutoff = 0.75,
                        subsets_mito_percent = "subsets_mito_percent") {
    metrics <- as.data.frame(colData(sce))
    
    if (is.null(model)) {
        warning("call 'mixtureModel' explicitly to get stable model features")
        model <- mixtureModel(sce, model_type = "one_dimensional")
    }
    
    intercept1 <- parameters(model, component = 1)[1]
    intercept2 <- parameters(model, component = 2)[1]
    if (intercept1 > intercept2) {
        compromised_dist <- 1
        intact_dist <- 2
    } else {
        compromised_dist <- 2
        intact_dist <- 1
    }
    
    post <- posterior(model)
    prob_compromised <- post[, compromised_dist]
    keep <- prob_compromised <= posterior_cutoff
    
    metrics <- cbind(metrics, prob_compromised = prob_compromised, keep = keep)
    
    min_discard <- min(metrics[!metrics$keep, ]$subsets_mito_percent)
    max_keep <- max(metrics[metrics$keep, ]$subsets_mito_percent)
    
    if (min_discard <= max_keep) {
        stop("Filtering recommendations differ based on library size. This
             is expected when mixtureModel is run with default parameters. To
             get a single threshold based only on mitochondrial percentage,
             run mixtureModel(sce, model_type = \"one_dimensional\")")
    }
    
    min_discard
}
