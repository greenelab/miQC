#' filterCells
#'
#' Find those cells probabilistically determined to be compromised by the
#' mixture model and remove them from the dataset.
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
#' @param preserve_lower_bound (boolean) Ensures that no cells below the intact
#'   cell distribution are removed. This should almost always be left as true.
#'   Default = TRUE
#'
#' @param enforce_left_bound (boolean) Prevents a U-shape in the filtering plot.
#'   For the cell with the lowest mitochondrial fraction that is set to be
#'   discarded, it ensures that no cells with lower library complexity (further
#'   left) and higher mitochondrial percentage (further up) than it are kept.
#'   Default = TRUE
#'
#' @param verbose (boolean) Whether to report how many cells (columns) are being
#'   removed from the SingleCellExperiment object.
#'   Default = TRUE
#'
#' @return Returns a SingleCellExperiment object, the same as the input except
#'   with a new column in colData, prob_compromised, and all cells with greater
#'   than the set posterior probability removed from the dataset.
#'
#' @importFrom BiocParallel MulticoreParam
#' @importFrom SingleCellExperiment colData
#' @importFrom flexmix parameters posterior fitted
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
#' sce <- filterCells(sce, model)



filterCells <- function(sce, model = NULL, posterior_cutoff = 0.75,
                        preserve_lower_bound = TRUE,
                        enforce_left_bound = TRUE, verbose = TRUE) {
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
        intact_dist <- 1
        compromised_dist <- 2
    }

    post <- posterior(model)
    metrics$prob_compromised <- post[, compromised_dist]
    sce$prob_compromised <- metrics$prob_compromised
    metrics$keep <- metrics$prob_compromised <= posterior_cutoff

    if (preserve_lower_bound == TRUE) {
        predictions <- fitted(model)[, intact_dist]
        metrics$intact_prediction <- predictions
        metrics[metrics$subsets_mito_percent <
                    metrics$intact_prediction, ]$keep <- TRUE
    }

    if (enforce_left_bound == TRUE) {
        min_discard <- min(metrics[!metrics$keep, ]$subsets_mito_percent)
        min_index <- which(metrics$subsets_mito_percent == min_discard)[1]
        lib_complexity <- metrics[min_index, ]$detected
        metrics[metrics$detected <= lib_complexity &
                    metrics$subsets_mito_percent >= min_discard, ]$keep <- FALSE
    }

    if (verbose == TRUE) {
        to_remove <- length(which(metrics$keep == FALSE))
        total <- length(metrics$keep)
        cat("Removing", to_remove, "out of", total, "cells.")
    }
    sce <- sce[, metrics$keep]

    sce
}
