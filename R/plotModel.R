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
#' @importFrom flexmix parameters posterior fitted
#' @importFrom ggplot2 ggplot aes labs geom_point geom_line ylim
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
#' plotModel(sce, model)



plotModel <- function(sce, model = NULL, detected = "detected",
                        subsets_mito_percent = "subsets_mito_percent") {
    metrics <- as.data.frame(colData(sce))

    if (is.null(model)) {
        warning("call 'mixtureModel' explicitly to get stable model features")
        model <- mixtureModel(sce)
    }

    predictions <- fitted(model)
    Comp.1 <- predictions[, 1]
    Comp.2 <- predictions[, 2]
    fitted_model_comp1 <- as.data.frame(cbind(detected = metrics$detected,
                                            Comp.1 = Comp.1))
    fitted_model_comp1 <- subset(fitted_model_comp1, 
                                 fitted_model_comp1$Comp.1 >= 0)
    fitted_model_comp2 <- as.data.frame(cbind(detected = metrics$detected,
                                              Comp.2 = Comp.2))
    fitted_model_comp2 <- subset(fitted_model_comp2, 
                                 fitted_model_comp2$Comp.2 >= 0)

    intercept1 <- parameters(model, component = 1)[1]
    intercept2 <- parameters(model, component = 2)[1]
    if (intercept1 > intercept2) {
        compromised_dist <- 1
    } else {
        compromised_dist <- 2
    }

    post <- posterior(model)
    prob_compromised <- post[, compromised_dist]
    metrics <- cbind(metrics, prob_compromised = prob_compromised)

    p <- ggplot(metrics, aes(x = detected, y = subsets_mito_percent,
                                colour = prob_compromised)) +
        labs(x = "Unique genes found", y = "Percent reads mitochondrial",
                color = "Probability\ncompromised") +
        geom_point() +
        geom_line(data = fitted_model_comp1, inherit.aes = FALSE,
                    aes(x = detected, y = Comp.1), lwd = 2) +
        geom_line(data = fitted_model_comp2, inherit.aes = FALSE,
                    aes(x = detected, y = Comp.2), lwd = 2) +
        ylim(0, max(metrics$subsets_mito_percent))

    p
}
