#' TMM normalization of single cell data sets
#'
#' @description
#' Performs for a list of SummarizedExperiments, TMM normalization scaled by a factor specified by the user.
#'
#' @details
#'
#' This function estimates TMM normalization factors and normalizes a list of gene count matrices, data
#' is additionally scaled using a factor specified by the user and `log1p()` transformed
#'
#' @param pb_dat_list List of SummarizedExperiment generated from `MOFAcellulaR::filt_profiles()`
#' @param scale_factor Numeric. Scaled counts are multiplied by this factor before log transformation
#' @return A named list of SummarizedExperiments per cell type provided with normalized log transformed data in their `logcounts` assay
#' @export
#'
#' @import SummarizedExperiment
#' @import S4Vectors
#' @import edgeR
#' @import purrr
#'
#' @examples
#' inputs_dir <- base::system.file("data", package = "MOFAcellulaR")
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
tmm_trns <- function(pb_dat_list, scale_factor = 1000000) {

  pb_dat_red <- purrr::map(pb_dat_list, function(x) {
    all_nf <- edgeR::calcNormFactors(x, method = "TMM")
    sfs <- all_nf$samples$lib.size * all_nf$samples$norm.factors
    pb <- base::sweep(assay(x, "counts"), MARGIN = 2, sfs, FUN = "/")
    SummarizedExperiment::assay(x, "logcounts") <- base::log1p(pb * scale_factor)

    return(x)

  })

  return(pb_dat_red)

}
