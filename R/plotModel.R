#' plotModel
#'
#' Function to plot quality characteristics of cells in dataset, parameters of
#' compromised and intact distributions, and posterior probability of each cell
#' belonging to the compromised distribution.
#'
#' @param sce (SingleCellExperiment) Input data object.
#'
#' @param model (flexmix) Output of mixtureModel function, which should be
#'   explicitly called first to ensure stability of model parameters.
#'   Default = NULL.
#'
#' @return Returns a ggplot object. Additional plot elements can be added as
#'   ggplot elements (e.g. title, customized formatting, etc).
#'
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom flexmix parameters posterior fitted
#' @importFrom ggplot2 ggplot aes labs geom_point geom_line ylim
#'
#' @export
#'
#' @examples
#' library(scRNAseq)
#' sce <- ZeiselBrainData()
#' mt_genes <- grepl("^mt-",  rownames(sce))
#' feature_ctrls <- list(mito = rownames(sce)[mt_genes])
#' sce <- addPerCellQC(sce, subsets = feature_ctrls, BPPARAM = BiocParallel::MulticoreParam())
#' model <- mixtureModel(sce)
#' sce <- plotModel(sce, model)



plotModel <- function(sce, model = NULL) {
    metrics <- as.data.frame(colData(sce))

    if(is.null(model)) {
        warning("call 'mixtureModel' explicitly to get stable model features")
        model <- mixtureModel(sce)
    }

    fitted_models <- as.data.frame(cbind(metrics$detected, fitted(model)))

    intercept1 <- parameters(model, component = 1)[1]
    intercept2 <- parameters(model, component = 2)[1]
    if (intercept1 > intercept2) {
        compromised_dist <- 1
    } else {
        compromised_dist <- 2
    }

    post <- posterior(model)
    metrics$prob_compromised <- post[, compromised_dist]

    p <- ggplot(metrics, aes(x = detected, y = subsets_mito_percent,
                                colour = prob_compromised)) +
        labs(x = "Unique genes found", y = "Percent reads mitochondrial",
                color = "Probability\ncompromised") +
        geom_point() +
        geom_line(data = fitted_models, inherit.aes = FALSE,
                    aes(x = V1, y = Comp.1), lwd = 2) +
        geom_line(data = fitted_models, inherit.aes = FALSE,
                    aes(x = V1, y = Comp.2), lwd = 2) +
        ylim(0, max(metrics$subsets_mito_percent))

    p
}
