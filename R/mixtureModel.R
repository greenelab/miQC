#' mixtureModel
#'
#' Function to fit a two-distribution mixture model on a SingleCellExperiment
#' object.
#'
#'
#' @param sce (SingleCellExperiment) Input data object.
#'
#' @param model_type (character) What type of model to generate. A linear
#'   mixture model ("linear") based on mitochondrial percentage and library
#'   complexity is recommended. B-spline ("spline") and two-degree polynomial
#'   ("polynomial") models are also supported. For a simpler model, a
#'   one-dimensional gaussian mixture model ("one_dimensional") based on 
#'   mitochondrial percentage only is available.
#'   Default = "linear".
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
#' @param group (character) Optional column name in \code{colData(sce)} used to
#'   fit separate miQC models per group (for example, annotated cell type or
#'   sample). When \code{NULL}, one model is fit across all cells.
#'
#' @return Returns a flexmix object with mixture model parameters, which is used
#'   to calculate posterior probability for each cell being compromised and make
#'   final filtering decisions. When \code{group} is provided, returns a named
#'   list of flexmix models, one per group.
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom flexmix flexmix
#' @importFrom splines bs
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



mixtureModel <- function(sce, model_type = "linear", detected = "detected",
                            subsets_mito_percent = "subsets_mito_percent",
                            group = NULL) {
    metrics <- as.data.frame(colData(sce))

    .fit_model <- function(metrics, model_type) {
        if (model_type == "linear") {
            model <- flexmix(subsets_mito_percent~detected,
                                data = metrics, k = 2)
        } else if (model_type == "spline") {
            model <- flexmix(subsets_mito_percent~bs(detected),
                                data = metrics, k = 2)
        } else if (model_type == "polynomial") {
            model <- flexmix(subsets_mito_percent~poly(detected, degree = 2),
                                data = metrics, k = 2)
        } else if (model_type == "one_dimensional") {
            model <- flexmix(subsets_mito_percent~1, data = metrics, k = 2)
        } else {
            stop("Unsupported model_type: ", model_type)
        }

        if (length(model@components) < 2) {
            warning("Unable to identify two distributions. Use plotMetrics function
                    to confirm assumptions of miQC are met.")
        }

        model
    }

    if (!is.null(group)) {
        if (!group %in% colnames(metrics)) {
            stop("group must be a valid column name in colData(sce).")
        }
        if (any(is.na(metrics[[group]]))) {
            stop("group column must not contain NA values.")
        }

        indices <- split(seq_len(nrow(metrics)), as.character(metrics[[group]]))
        model <- lapply(indices, function(i) .fit_model(metrics[i, , drop = FALSE],
                                                        model_type = model_type))
        class(model) <- c("miQC_grouped_model", class(model))
        attr(model, "group_column") <- group
        return(model)
    }

    model <- .fit_model(metrics, model_type = model_type)

    model
}
