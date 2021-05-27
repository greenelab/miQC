#' plotMetrics
#'
#' A function to plot the QC parameters used for a miQC model, number of unique
#' genes expressed and percent mitochondrial reads. This function can be run
#' before calling mixtureModel() to assess if miQC is appropriate given the data
#' distribution. See vignette for examples of cases where miQC is and isn't a
#' good choice for filtering.
#'
#' @param sce (SingleCellExperiment) Input data object.
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
#' @param palette (character) Specifies the color to plot cells as. Default is
#'   "#33ADFF".
#'
#' @return Returns a ggplot object. Additional plot elements can be added as
#'   ggplot elements (e.g. title, customized formatting, etc).
#'
#' @importFrom SingleCellExperiment colData
#' @importFrom ggplot2 ggplot aes labs geom_point
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
#' plotMetrics(sce)


plotMetrics <- function(sce, detected = "detected",
                        subsets_mito_percent = "subsets_mito_percent",
                        palette = "#33ADFF") {
    metrics <- as.data.frame(colData(sce))

    p <- ggplot(metrics, aes(x = detected, y = subsets_mito_percent)) +
        labs(x = "Unique genes found", y = "Percent reads mitochondrial") +
        geom_point(colour = palette)

    p
}
