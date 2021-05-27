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
#' @param keep_all_below_boundary (boolean) Ensures that no cells below the
#'   intact cell distribution are removed. This should almost always be set
#'   to true.
#'   Default = TRUE
#'
#' @param enforce_left_cutoff (boolean) Prevents a U-shape in the filtering 
#'   plot. Identifies the cell with the lowest mitochondrial fraction that is
#'   set to be discarded, it ensures that no cells with lower library 
#'   complexity (further left) and higher mitochondrial percentage (further up)
#'   than it are kept.
#'   Default = TRUE
#'
#' @param palette (character) Color palette. A vector of length two containing
#'   custom colors.
#'   Default = c("#999999", "#E69F00").
#'
#' @param detected (character) Column name in sce giving the number of unique
#'   genes detected per cell. This name is inherited by default from scater's
#'   addPerCellQC() function.
#'
#' @param subsets_mito_percent (character) Column name in sce giving the
#'   percent of reads mapping to mitochondrial genes. This name is inherited
#'   from scater's addPerCellQC() function, provided the subset "mito" with
#'   names of all mitochondrial genes is passed in. See examples for details.
#'
#' @return Returns a ggplot object. Additional plot elements can be added as
#'   ggplot elements (e.g. title, customized formatting, etc).
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom flexmix parameters posterior
#' @importFrom ggplot2 ggplot aes labs geom_point scale_color_manual theme
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
#' model <- mixtureModel(sce)
#' plotFiltering(sce, model)

plotFiltering <- function(sce, model = NULL, posterior_cutoff = 0.75,
                            keep_all_below_boundary = TRUE,
                            enforce_left_cutoff = TRUE,
                            palette = c("#999999", "#E69F00"),
                            detected = "detected",
                            subsets_mito_percent = "subsets_mito_percent") {
    metrics <- as.data.frame(colData(sce))

    if (is.null(model)) {
        warning("call 'mixtureModel' explicitly to get stable model features")
        model <- mixtureModel(sce)
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

    if (keep_all_below_boundary == TRUE) {
        predictions <- fitted(model)[, intact_dist]
        metrics$intact_prediction <- predictions
        metrics[metrics$subsets_mito_percent <
                    metrics$intact_prediction, ]$keep <- TRUE
    }

    if (enforce_left_cutoff == TRUE) {
        min_discard <- min(metrics[!metrics$keep, ]$subsets_mito_percent)
        min_index <- which(metrics$subsets_mito_percent == min_discard)[1]
        lib_complexity <- metrics[min_index, ]$detected
        metrics[metrics$detected <= lib_complexity &
                    metrics$subsets_mito_percent >= min_discard, ]$keep <- FALSE
    }

    p <- ggplot(metrics, aes(x = detected, y = subsets_mito_percent,
                                colour = keep)) +
        labs(x = "Unique genes found", y = "Percent reads mitochondrial",
                color = "Keep") +
        scale_color_manual(values = palette) + geom_point()

    p
}
