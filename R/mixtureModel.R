#' mixtureModel
#'
#' Function to fit a two-distribution mixture model on a SingleCellExperiment
#' object.
#'
#'
#' @param sce (SingleCellExperiment) Input data object.
#'
#' @param model_type (character) What type of model to generate. A linear
#'   mixture model ("linear") is recommended, but currently b-spline ("spline")
#'   and two-degree polynomial ("polynomial") are also supported
#'   Default = "linear".
#'
#' @return Returns a flexmix object with mixture model parameters, which is used
#'   to calculate posterior probability for each cell being compromised and make
#'   final filtering decisions.
#'
#' @importFrom BiocParallel MulticoreParam
#' @importFrom SingleCellExperiment colData
#' @importFrom flexmix flexmix
#' @importFrom splines bs
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



mixtureModel <- function(sce, model_type = "linear") {
    metrics <- as.data.frame(colData(sce))

    if (model_type == "linear") {
        model <- flexmix(subsets_mito_percent~detected,
                            data = metrics, k = 2)
    } else if (model_type == "spline") {
        model <- flexmix(subsets_mito_percent~bs(detected),
                            data = metrics, k = 2)
    } else if (model_type == "polynomial") {
        model <- flexmix(subsets_mito_percent~poly(detected, degree = 2),
                            data = metrics, k = 2)
    }

    model
}
