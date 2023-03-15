#' Indentify highly variable genes
#'
#' @description
#' Identifies highly variable features from a log-normalized count matrix or filters matrices by a list of genes
#' provided by the user.
#'
#' @details
#'
#' This function estimates highly variable genes per cell type using `scran::getTopHVGs`. Alternatively, this function
#' allows the user to provide the features to be used in each cell type. If prior genes are used, for cell types
#' where this information is missing, highly variable genes will be calculated
#'
#' @param pb_dat_list List of SummarizedExperiment generated from `MOFAcellulaR::filt_profiles()`
#' @param prior_hvg NULL by default. Alternatively, a named list with a character vector containing features to select.
#' @param var.threshold Numeric. Inherited from `scran::getTopHVGs()`. Minimum threshold on the metric of variation
#' @return A named list of SummarizedExperiments per cell type provided with filtered normalized log transformed data in their `logcounts` assay
#' @export
#'
#' @import SummarizedExperiment
#' @import S4Vectors
#' @import edgeR
#' @import purrr
#' @import scran
#'
#' @examples
#' inputs_dir <- base::system.file("extdata", package = "MOFAcellulaR")
#' load(file.path(inputs_dir, "testpbcounts.rda"))
#' load(file.path(inputs_dir, "testcoldata.rda"))
#'
#' pb_obj <- create_init_exp(counts = testpbcounts,
#'                           coldata = testcoldata)
#'
#' ct_list <- filt_profiles(pb_dat = pb_obj,
#'                          cts = c("Fib","CM"),
#'                          ncells = 5,
#'                          counts_col = "cell_counts",
#'                          ct_col = "cell_type")
#'
#' ct_list <- filt_gex_byexpr(pb_dat_list = ct_list,
#'                            min.count = 5,
#'                            min.prop = 0.25)
#'
#' ct_list <- tmm_trns(pb_dat_list = ct_list,
#'                     scale_factor = 1000000)
#'
#' ct_list <- filt_gex_byhvg(pb_dat_list = ct_list,
#'                           prior_hvg = NULL,
#'                           var.threshold = 0)
#'
filt_gex_byhvg <- function(pb_dat_list,
                           prior_hvg = NULL,
                           var.threshold = 1) {
# If there is no prior, calculate hvgs for every member of the list
  if(is.null(prior_hvg)) {

    pb_dat_red <- purrr::map(pb_dat_list, function(x) {
      hvg <- scran::getTopHVGs(x,var.threshold = var.threshold)
      return(x[hvg, ])
    })

    return(pb_dat_red)

  } else {

    cts_in_data <- purrr::set_names(names(pb_dat_list))
    cts_in_prior <- purrr::set_names(names(prior_hvg))

    in_cts <- cts_in_data[cts_in_data %in% cts_in_prior]
    out_cts <- cts_in_data[!cts_in_data %in% cts_in_prior]

    in_cts_data <- pb_dat_list[in_cts]
# For cell types with prior, filter with available genes
    for(ct in in_cts) {
      ct_genes <- in_cts_data[[ct]] %>% rownames()
      ct_genes <- ct_genes[ct_genes %in% prior_hvg[[ct]]]
      in_cts_data[[ct]] <- in_cts_data[[ct]][ct_genes,]
    }

    if(length(out_cts) == 0) {

      return(in_cts_data)

    } else {
# For the rest, calculate hvgs
      out_cts_data <- pb_dat_list[out_cts]

      out_cts_data <- purrr::map(out_cts_data, function(x) {
        hvg <- scran::getTopHVGs(x,var.threshold = var.threshold)
        return(x[hvg, ])
      })

      return(c(in_cts_data, out_cts_data))

    }
  }
}

#' Filter background expression of marker genes
#'
#' @description
#' For a collection of matrices, we exclude features that are considered background based on prior
#' knowledge of marker genes
#'
#' @details
#' Performs filtering of highly variable genes (after data transformation).
#' This is based on marker genes. The assumption is that background gene expression can be traced
#' by expression of cell type marker genes in cell types which shouldn't express the gene.
#' Marker genes will be only kept in the matrix if they are expressed in the expected cell type
#'
#' @param pb_dat_list List of SummarizedExperiment generated from `MOFAcellulaR::filt_profiles()`
#' @param prior_mrks A named list providing marker genes per cell type. Names should be identical as in `pb_dat_list`
#' @return A named list of SummarizedExperiments per cell type provided with filtered normalized log transformed data in their `logcounts` assay
#' @export
#'
#' @import SummarizedExperiment
#' @import S4Vectors
#' @import edgeR
#' @import purrr
#' @import scran
#' @import tidyr
#' @import tibble
#' @import dplyr
#'
#' @examples
#' inputs_dir <- base::system.file("extdata", package = "MOFAcellulaR")
#' load(file.path(inputs_dir, "testpbcounts.rda"))
#' load(file.path(inputs_dir, "testcoldata.rda"))
#'
#' pb_obj <- create_init_exp(counts = testpbcounts,
#'                           coldata = testcoldata)
#'
#' ct_list <- filt_profiles(pb_dat = pb_obj,
#'                          cts = c("Fib","CM"),
#'                          ncells = 5,
#'                          counts_col = "cell_counts",
#'                          ct_col = "cell_type")
#'
#' ct_list <- filt_gex_byexpr(pb_dat_list = ct_list,
#'                            min.count = 5,
#'                            min.prop = 0.25)
#'
#' ct_list <- tmm_trns(pb_dat_list = ct_list,
#'                     scale_factor = 1000000)
#'
#' ct_list <- filt_gex_byhvg(pb_dat_list = ct_list,
#'                           prior_hvg = NULL,
#'                           var.threshold = 0)
#'
#' prior_hvg_test <- list("CM" = c("TTN"),
#'                        "Fib" = c("POSTN"))
#'
#' ct_list <- filt_gex_bybckgrnd(pb_dat_list = ct_list,
#'                               prior_mrks = prior_hvg_test)
#'
filt_gex_bybckgrnd <- function(pb_dat_list, prior_mrks) {

  # Current genes per view
  ct_genes <- purrr::map(pb_dat_list, rownames) %>%
    tibble::enframe("view","gene") %>%
    tidyr::unnest(.data$gene)

  prior_mrks_df <- prior_mrks %>%
    tibble::enframe("view_origin","gene") %>%
    tidyr::unnest(.data$gene) %>%
    dplyr::mutate(marker_gene = TRUE)

  # Here are genes that aren't cell type markers
  ok_genes <- ct_genes %>%
    dplyr::left_join(prior_mrks_df, by = "gene") %>%
    dplyr::filter(is.na(.data$marker_gene)) %>%
    dplyr::select_at(c("view", "gene"))

  # Here are genes selected as HVG that are marker
  # genes, we will keep only genes if they appear
  # in the right cell
  not_bckground_genes <- ct_genes %>%
    dplyr::left_join(prior_mrks_df, by = "gene") %>%
    stats::na.omit() %>%
    tidyr::unnest(c()) %>%
    dplyr::filter(.data$view == .data$view_origin) %>%
    dplyr::select_at(c("view", "gene"))

  clean_hvgs <- dplyr::bind_rows(ok_genes,
                          not_bckground_genes) %>%
    dplyr::group_by(view) %>%
    tidyr::nest() %>%
    dplyr::mutate(data = map(.data$data, ~.x[[1]])) %>%
    tibble::deframe()

  pb_dat_list <- pb_dat_list %>%
    filt_gex_byhvg(pb_dat_list = .,
                   prior_hvg = clean_hvgs,
                   var.threshold = NULL)

  return(pb_dat_list)
}
