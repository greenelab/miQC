#' plotFiltering
#'
#' Function to plot which cells will be kept and removed given their posterior
#' probability of belonging to the compromised distribution.
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
#' @param palette (character) Color palette. A vector of length two containing
#'   custom colors.
#'   Default = c("#999999", "#E69F00").
#'
#' @return Returns a ggplot object. Additional plot elements can be added as
#'   ggplot elements (e.g. title, customized formatting, etc).
#'
#' @importFrom BiocParallel MulticoreParam
#' @importFrom SingleCellExperiment colData
#' @importFrom flexmix parameters posterior
#' @importFrom ggplot2 ggplot aes labs geom_point scale_color_manual theme
#'
#' @export
#'
#' @examples
#' library(scRNAseq)
#' library(scater)
#' library(BiocParallel)
#' sce <- ZeiselBrainData()
#' mt_genes <- grepl("^mt-",  rownames(sce))
#' feature_ctrls <- list(mito = rownames(sce)[mt_genes])
#' sce <- addPerCellQC(sce, subsets = feature_ctrls, BPPARAM = MulticoreParam())
#' model <- mixtureModel(sce)
#' sce <- plotFiltering(sce, model)

plotFiltering <- function(sce, model = NULL, posterior_cutoff = 0.75,
                            palette = c("#999999", "#E69F00")) {
    metrics <- as.data.frame(colData(sce))

    if(is.null(model)) {
        warning("call 'mixtureModel' explicitly to get stable model features")
        model <- mixtureModel(sce)
    }

    intercept1 <- parameters(model, component = 1)[1]
    intercept2 <- parameters(model, component = 2)[1]
    if (intercept1 > intercept2) {
        compromised_dist <- 1
    } else {
        compromised_dist <- 2
    }

    post <- posterior(model)
    metrics$prob_compromised <- post[, compromised_dist]
    metrics$keep <- metrics$prob_compromised <= posterior_cutoff

    p <- ggplot(metrics, aes(x = detected, y = subsets_mito_percent,
                                colour = keep)) +
        labs(x = "Unique genes found", y = "Percent reads mitochondrial",
                color = "Keep") +
        scale_color_manual(values = palette) + geom_point()

    p
}