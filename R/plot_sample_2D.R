#' Visualize sample variability in a 2D space
#'
#' @description
#' Performs dimensionality reduction of factor scores for the visualization of sample-level variability
#'
#' @details
#' For a `MOFA2` model it performs a multidimensional scaling (MDS)
#' or uniform manifold approximation and projection (UMAP) of the factor
#' scores. It allows you to color samples based on a covariate available in
#' your meta-data.
#'
#' @param model A MOFA2 model.
#' @param group Boolean flag TRUE/FALSE, to specify if a grouped MOFA model is provided.
#' @param method A string specifying if "UMAP" or "MDS" should be performed
#' @param metadata A data frame containing the annotations of the samples included in the MOFA model.
#' @param sample_id_column A string character that refers to the column in `metadata` where the sample identifier is located.
#' @param color_by A string character that refers to the column in `metadata` where the covariate to be tested is located.
#' @param ... inherited parameters of `uwot::umap()`
#'
#'
#' @return A dataframe in a tidy format containing the manifold and the scatter plot
#' @export
#'
#' @import MOFA2
#' @import tibble
#' @import dplyr
#' @import uwot
#' @import stats
#' @import ggplot2
#'
#' @examples
#' inputs_dir <- base::system.file("extdata", package = "MOFAcellulaR")
#' model <- MOFA2::load_model(file.path(inputs_dir, "testmodel.hdf5"))
#' metadata <- readRDS(file.path(inputs_dir, "testmetadata.rds"))
#' UMAP_embedding <- plot_sample_2D(model = model,
#'                                  group = FALSE,
#'                                  method = "UMAP",
#'                                  metadata = metadata,
#'                                  sample_id_column = "sample",
#'                                  color_by = "batch")
plot_sample_2D <- function(model,
                           group = FALSE,
                           method = "UMAP",
                           metadata,
                           sample_id_column = "sample",
                           color_by,
                           ...) {

  # Get the matrix of factor scores
  if(group){

    factors <- MOFA2::get_factors(model, factors = "all") %>%
      base::do.call(base::rbind, .)

  } else {

    factors <- MOFA2::get_factors(model, factors = "all")[[1]]

  }

  # Make specific reductions
  if(method == "MDS") {
    # Using eucledian
    factors_mds <- stats::cmdscale(stats::dist(factors)) %>%
      base::as.data.frame()

    colnames(factors_mds) <- c("MDS1", "MDS2")

    factors_mds <- factors_mds %>%
      as.data.frame() %>%
      tibble::rownames_to_column(sample_id_column) %>%
      dplyr::left_join(metadata, by = sample_id_column) %>%
      dplyr::rename("color_col" = color_by)

    mds_plt <- ggplot2::ggplot(factors_mds, ggplot2::aes(x = .data$MDS1,
                                       y = .data$MDS2,
                                       color = .data$color_col)) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text = ggplot2::element_text(size =12)) +
      ggplot2::labs(color = "")

    plot(mds_plt)

    return(factors_mds)

  } else if(method == "UMAP") {
    # UMAP of factors
    factors_umap <- uwot::umap(factors, ...) #Inherited parameters
    colnames(factors_umap) <- c("UMAP_1", "UMAP_2")

    factors_umap <- factors_umap %>%
      as.data.frame() %>%
      tibble::rownames_to_column(sample_id_column) %>%
      dplyr::left_join(metadata, by = sample_id_column) %>%
      dplyr::rename("color_col" = color_by)

    umap_plt <- ggplot2::ggplot(factors_umap,
                                ggplot2::aes(x = .data$UMAP_1,
                                             y = .data$UMAP_2,
                                             color = color_col)) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text = ggplot2::element_text(size =12)) +
      ggplot2::labs(color = "")

    plot(umap_plt)

    return(factors_umap)

  }

}




















