#' Filter pseudobulk profiles
#'
#' @description
#' Filter pseudobulk profiles of specific cell types based on the number of cells from which they were generated.
#'
#' @details
#'
#' This function assumes that you have a SummarizedExperiment object
#' with information in `colData(object)` specifying the number of cells used
#' for each profile and the cell-type grouping this profiles. The function then will
#' select only the cell-types provided by the user and filter the profiles
#' with less cells as the ones specified.
#'
#' @param pb_dat SummarizedExperiment generated from `MOFAcellulaR::create_init_exp()`
#' @param cts A vector containing the names of cells to be used in the analysis
#' @param ncells Number of minimum cells of each pseudobulk profile
#' @param counts_col String pointing to the column in `colData(pb_dat)` where the number of cells per pseudobulk was stored
#' @param ct_col String pointing to the column in `colData(pb_dat)` where the cell-type category per pseudobulk was stored
#'
#' @return A named list of SummarizedExperiments per cell type provided with filtered pseudobulk profiles
#' @export
#'
#' @import SummarizedExperiment
#' @import S4Vectors
#' @import purrr
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
filt_profiles <- function(pb_dat,
                          cts,ncells = 50,
                          counts_col = "cell_counts",
                          ct_col = "cell_type"){

  # by n of cells

  ix <- base::which(SummarizedExperiment::colData(pb_dat)[,counts_col] >= ncells)
  pb_dat <- pb_dat[, ix]

  # by views of interest

  if(is.null(cts)) {

    cts <- purrr::set_names(SummarizedExperiment::colData(pb_dat)[,ct_col] %>%
                       unique())

  } else {

    cts <- purrr::set_names(cts)

  }

  pb_dat_list <- purrr::map(cts, function(ctype) {

    ix <- base::which(SummarizedExperiment::colData(pb_dat)[,ct_col] == ctype)

    return(pb_dat[,ix])

  })

  return(pb_dat_list)

}

#' Filter lowly expressed genes from pseudobulk profiles
#'
#' @description
#' Filter lowly expressed genes from pseudobulk profiles using `edgeR`.
#'
#' @details
#'
#' This function wraps `edgeR::filterByExpr()` to be applied to lists of SummarizedExperiments.It assumes
#' that all samples are part of the same group.
#'
#' @param pb_dat_list List of SummarizedExperiment generated from `MOFAcellulaR::filt_profiles()`
#' @param min.count Numeric, minimum counts per sample to be considered. Check `?edgeR::filterByExpr()` for details.
#' @param min.prop Numeric, minimum proportion of samples containing the minimum counts. Check `?edgeR::filterByExpr()` for details.
#'
#' @return A named list of SummarizedExperiments per cell type provided with filtered pseudobulk profiles
#' @export
#'
#' @import SummarizedExperiment
#' @import S4Vectors
#' @import edgeR
#' @import purrr
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
filt_gex_byexpr <- function(pb_dat_list, min.count, min.prop) {

  pb_dat_red <- purrr::map(pb_dat_list, function(x) {

    useful_genes <- edgeR::filterByExpr(x,
                                        min.count = min.count,
                                        min.prop = min.prop)

    return(x[useful_genes, ])
  })

  return(pb_dat_red)

}

#' Filter views with not enough features
#'
#' @description
#' Filter complete views based on their number of profiled features
#'
#' @details
#'
#' This function allows the user to control the number of minimum features per view
#'
#' @param pb_dat_list List of SummarizedExperiment generated from `MOFAcellulaR::filt_profiles()`
#' @param ngenes Numeric, minimum number of features per view.
#'
#' @return A named list of SummarizedExperiments per cell type provided with filtered pseudobulk profiles
#' @export
#'
#' @import SummarizedExperiment
#' @import S4Vectors
#' @import purrr
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
#' ct_list <- filt_views_bygenes(pb_dat_list = ct_list,
#'                               ngenes = 15)
#'
filt_views_bygenes <- function(pb_dat_list,
                               ngenes) {

  view_count <- pb_dat_list %>%
    purrr::map_dbl(~ .x %>% nrow())

  view_count <- view_count[view_count >= ngenes]

  views <- names(view_count)

  return(pb_dat_list[views])

}

#' Filter views with not enough samples
#'
#' @description
#' Filter complete views based on their number of profiled samples
#'
#' @details
#'
#' This function allows the user to control the number of minimum samples per view
#'
#' @param pb_dat_list List of SummarizedExperiment generated from `MOFAcellulaR::filt_profiles()`
#' @param nsamples Numeric, minimum number of samples per view.
#'
#' @return A named list of SummarizedExperiments per cell type provided with filtered pseudobulk profiles
#' @export
#'
#' @import SummarizedExperiment
#' @import S4Vectors
#' @import purrr
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
#' ct_list <- filt_views_bysamples(pb_dat_list = ct_list,
#'                               nsamples = 2)
#'
filt_views_bysamples <- function(pb_dat_list,
                                 nsamples) {

  view_count <- pb_dat_list %>%
    purrr::map_dbl(~ .x %>% ncol())

  view_count <- view_count[view_count >= nsamples]

  views <- names(view_count)

  return(pb_dat_list[views])

}

#' Filter samples that have a limited coverage of features in a view
#'
#' @description
#' Excludes samples from views with a coverage lower than the one stated
#'
#' @details
#'
#' This function allows the user to exclude samples that have a limited number of features that could create potential outliers
#'
#' @param pb_dat_list List of SummarizedExperiment generated from `MOFAcellulaR::filt_profiles()`
#' @param prop_coverage Numeric, minimum proportion of coverage of features (different from 0).
#'
#' @return A named list of SummarizedExperiments per cell type provided with filtered pseudobulk profiles
#' @export
#'
#' @import SummarizedExperiment
#' @import S4Vectors
#' @import purrr
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
#' ct_list <- filt_views_bygenes(pb_dat_list = ct_list,
#'                               ngenes = 15)
#'
#' ct_list <- filt_samples_bycov(pb_dat_list = ct_list,
#'                               prop_coverage = 0.9)
#'
#'
filt_samples_bycov <- function(pb_dat_list, prop_coverage = 0.9) {

  pb_dat_list_new <- purrr::map(pb_dat_list, function(pb_obj) {

    mat <- SummarizedExperiment::assay(pb_obj, "counts")

    sample_vect <- (mat != 0) %>% colSums(.)/nrow(mat)

    sel_samples <- names(sample_vect[which(sample_vect >= prop_coverage)])

    return(pb_obj[,sel_samples])

  })

  return(pb_dat_list_new)
}
