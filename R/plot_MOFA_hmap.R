#' Visualize MOFA multicellular model
#'
#' @description
#' Plots a heatmap with factor scores annotated by sample's information and factor's statistics
#'
#' @details
#' This function summarizes a `MOFA2` model by plotting and clustering the factor scores across samples.
#' Additionally, it allows you to annotate each sample with categorical or continous variables. Finally,
#' for each factor, the amount of explained variance captured for each cell type is shown. If provided,
#' summary of the association statistics and sample variables can be provided
#'
#' @param model A MOFA2 model.
#' @param group Boolean flag TRUE/FALSE, to specify if a grouped MOFA model is provided.
#' @param metadata A data frame containing the annotations of the samples included in the MOFA model.
#' @param sample_id_column A string character that refers to the column in `metadata` where the sample identifier is located.
#' @param sample_anns A vector containing strings that refer to the columns in `metadata` to be used to annotate samples
#' @param assoc_list A named list collecting results of `MOFAcellulaR::get_associations()`
#' @param col_rows A named list of lists containing at the first index, all sample_anns used,
#' and at the second index, all levels of an annotation within the first index.
#' As required by `ComplexHeatmap::rowAnnotation` col parameter.
#'
#' @return A dataframe in a tidy format containing the manifold and the scatter plot
#' @export
#'
#' @import MOFA2
#' @import tibble
#' @import dplyr
#' @import purrr
#' @import tidyr
#' @import circlize
#' @import grDevices
#' @import ComplexHeatmap
#' @import grid
#'
#' @examples
#' inputs_dir <- base::system.file("extdata", package = "MOFAcellulaR")
#' model <- MOFA2::load_model(file.path(inputs_dir, "testmodel.hdf5"))
#' metadata <- readRDS(file.path(inputs_dir, "testmetadata.rds"))
#' metadata$var <- stats::rnorm(nrow(metadata))
#'
#' categorical_assoc <- get_associations(model = model,
#'                                       metadata = metadata,
#'                                       sample_id_column = "sample",
#'                                       test_variable = "patient_group",
#'                                       test_type = "categorical",
#'                                       group = FALSE)
#'
#' continuous_assoc <- get_associations(model = model,
#'                                      metadata = metadata,
#'                                      sample_id_column = "sample",
#'                                      test_variable = "var",
#'                                      test_type = "continous",
#'                                      group = FALSE)
#'
#'
#' assoc_list = list("categorical" = categorical_assoc, "continous" = continuous_assoc)
#'
#' plot_MOFA_hmap(model = model,
#'                group = FALSE,
#'                metadata = metadata,
#'                sample_id_column = "sample",
#'                sample_anns = c("patient_group", "batch", "var"),
#'                assoc_list = assoc_list)

plot_MOFA_hmap <- function(model,
                           group = FALSE,
                           metadata,
                           sample_id_column = "sample",
                           sample_anns,
                           assoc_list = NULL,
                           col_rows = NULL) {

  # Build general aesthetics
  # aesthetics definition of the borders
  ht_opt$ROW_ANNO_PADDING <- grid::unit(2.5, "mm")
  ht_opt$COLUMN_ANNO_PADDING <- grid::unit(2.5, "mm")

  # Get the factor matrix to plot
  factor_matrix <- MOFA2::get_factors(model,
                                      factors = "all") %>%
    do.call(rbind, .)

  # Gradient colors --------------------------------------------------------------
  # For factor scores
  max_fact <- factor_matrix %>%
    max() %>%
    abs()

  col_fun_fact <- circlize::colorRamp2(seq((max_fact + 0.5) * -1,
                                           max_fact + 0.5,
                                           length = 50),
                                       grDevices::hcl.colors(50,"Green-Brown",rev = T))

  # 1. Row annotations (patient categories)
  # Here we find the right order in the meta-data -----------------------------------
  row_anns <- MOFA2::get_factors(model, factors = "all") %>%
    do.call(rbind, .) %>%
    rownames() %>%
    tibble::enframe(value = sample_id_column) %>%
    dplyr::select_at(sample_id_column) %>%
    dplyr::left_join(metadata, by = sample_id_column) %>%
    dplyr::select_at(sample_anns) %>%
    as.data.frame()

  # Patient annotations -------------------------------------------------------------
  if(is.null(col_rows)) {

    row_ha <- ComplexHeatmap::rowAnnotation(df = as.data.frame(row_anns),
                                            gap = grid::unit(2.5, "mm"),
                                            border = TRUE)

  } else {

    row_ha <- ComplexHeatmap::rowAnnotation(df = as.data.frame(row_anns),
                                            gap = grid::unit(2.5, "mm"),
                                            border = TRUE,
                                            col = col_rows)

  }

  # 2. Column annotations (R2 and associations)

  # 2.1 R2
  # Add explain variance per view, it must be a matrix, with
  # factors in row (ordered) and views in columns
  r2_list <- model@cache$variance_explained$r2_per_factor

  # If you have a grouped model, then we rename the views with the grouping
  if(group) {

    names_groups <- names(r2_list)

    r2_list <- purrr::map2(r2_list, names_groups, function(dat, nam) {

      dat_mod <- dat
      colnames(dat_mod) <- paste0(nam,"_", colnames(dat_mod))
      return(dat_mod)

    })

  }

  # Finally, homogeneous processing
  r2_per_factor <- r2_list %>%
    do.call(base::cbind, .)

  # Gradient colors --------------------------------------------------------------
  # For R2
  max_r2 <- r2_per_factor %>%
    max()

  if(max_r2 < 90) { # Add tolerance buffer
    col_fun_r2 <- circlize::colorRamp2(seq(0, max_r2 + 5, length = 50),
                                       grDevices::hcl.colors(50,"Oranges",rev = T))
  } else {
    col_fun_r2 <- circlize::colorRamp2(seq(0, 100, length = 50),
                                       grDevices::hcl.colors(50,"Oranges",rev = T))
  }


  # 2.2 Associations with covariates
  # If you provide a list of associations...
  if(!is.null(assoc_list)) {

    # This is a matrix, with factors in rows and p-values
    # for tested covariates in columns
    assoc_pvals <- assoc_list %>%
      tibble::enframe(name = "test") %>%
      tidyr::unnest(c(value)) %>%
      dplyr::mutate(log_adjpval = -log10(.data$adj_pvalue)) %>%
      dplyr::select(test, Factor, log_adjpval) %>%
      tidyr::pivot_wider(names_from = test,
                         values_from = log_adjpval) %>%
      dplyr::select(-Factor) %>%
      as.matrix()

    # Define aesthetics so as not to repeat the calculation
    # Association p-values
    col_fun_assoc <-  circlize::colorRamp2(seq(0, max(assoc_pvals) + 0.5, length = 20),
                                           grDevices::hcl.colors(20,"Purples",rev = T))


    column_ha <- ComplexHeatmap::HeatmapAnnotation("R2" = r2_per_factor,
                                  "pvalue" = assoc_pvals,
                                  gap = grid::unit(2.5, "mm"),
                                  border = TRUE,
                                  col = list(R2 = col_fun_r2,
                                             pvalue = col_fun_assoc))

  } else { # Just add R2 info

    column_ha <- ComplexHeatmap::HeatmapAnnotation("R2" = r2_per_factor,
                                   gap = grid::unit(2.5, "mm"),
                                   border = TRUE,
                                   col = list(R2 = col_fun_r2))

  }

  # Build the final heatmap
  scores_hmap <- ComplexHeatmap::Heatmap(factor_matrix,
                         name = "factor_scores",
                         right_annotation = row_ha,
                         top_annotation = column_ha,
                         cluster_columns = FALSE,
                         show_row_dend = TRUE,
                         show_row_names = FALSE,
                         border = TRUE,
                         gap = grid::unit(2.5, "mm"),
                         col = col_fun_fact)

  ComplexHeatmap::draw(scores_hmap)


  }


















