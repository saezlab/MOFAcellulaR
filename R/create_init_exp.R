#' Single-cell SummarizedExperiment object
#'
#' @description
#' Creates a SummarizedExperiment object necessary to build a multi-view representation.
#'
#' @details
#'
#' This function is the first step for a multicellular factor analysis.
#' It collects in a single object the pseudobulk counts of a single cell experiment
#' and its annotations.
#'
#' @param counts Named numeric matrix with features in rows and samples in columns.
#' @param coldata A data frame containing the annotations of the samples.
#'
#' @return SummarizedExperiment with provided data
#' @export
#'
#' @import SummarizedExperiment
#' @import S4Vectors
#'
#' @examples
#' inputs_dir <- base::system.file("extdata", package = "MOFAcellulaR")
#' load(file.path(inputs_dir, "testpbcounts.rda"))
#' load(file.path(inputs_dir, "testcoldata.rda"))
#' pb_obj <- create_init_exp(counts = testpbcounts,  coldata = testcoldata)
create_init_exp <- function(counts, coldata) {

  pb_dat <- SummarizedExperiment::SummarizedExperiment(assays = list("counts" = counts),
                                                       colData = S4Vectors::DataFrame(coldata))

  return(pb_dat)
}

#' Create MOFA-ready dataframe
#'
#' @description
#' Creates from a list of SummarizedExperiments a multi-view representation for MOFA
#'
#' @details
#'
#' This function is the last data preparation step for a multicellular factor analysis.
#' It collects a collection of cell-type-specific SummarizedExperiments into a
#' single data frame ready to be used in MOFA. Features are modified
#' so as to reflect their cell type of origin.
#'
#' @param pb_dat_list List of SummarizedExperiment generated from `MOFAcellulaR::filt_profiles()`
#'
#' @return Data frame in a multiview representation
#' @export
#'
#' @import SummarizedExperiment
#' @import S4Vectors
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
#' pb_obj <- create_init_exp(counts = testpbcounts,  coldata = testcoldata)
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
#' multiview_dat <- pb_dat2MOFA(pb_dat_list = ct_list)
pb_dat2MOFA <- function(pb_dat_list) {

  pb_red <- purrr::map(pb_dat_list, function(x) {

    dat <- SummarizedExperiment::assay(x, "logcounts")

    colnames(dat) <- SummarizedExperiment::colData(x)[,"donor_id"]

    dat %>%
      base::as.data.frame() %>%
      tibble::rownames_to_column("feature") %>%
      tidyr::pivot_longer(-.data$feature, names_to = "sample", values_to = "value")

  }) %>%
    tibble::enframe(name = "view") %>%
    tidyr::unnest(cols = c(value)) %>%
    dplyr::mutate(feature = paste0(.data$view, "_", .data$feature))

  return(pb_red)

}


#' Center view-data
#'
#' @description
#' Centers each element of a list of SummarizedExperiments
#'
#' @details
#'
#' Given that the MOFA model in general uses centered data, when interested in projecting new
#' data to a new manifold, it is needed to perform centering.
#'
#' @param pb_dat_list List of SummarizedExperiment generated from `MOFAcellulaR::filt_profiles()`
#'
#' @return A named list of SummarizedExperiments per cell type provided with centered pseudobulk profiles
#' @export
#'
#' @import SummarizedExperiment
#' @import S4Vectors
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
#' pb_obj <- create_init_exp(counts = testpbcounts,  coldata = testcoldata)
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
#' ct_list <- center_views(pb_dat_list = ct_list)
#'
#' multiview_dat <- pb_dat2MOFA(pb_dat_list = ct_list)
center_views <- function(pb_dat_list) {

  pb_red <- purrr::map(pb_dat_list, function(x) {

    dat <- SummarizedExperiment::assay(x, "logcounts") %>% t() # genes in columns now

    scaled_dat <- base::scale(dat, scale = FALSE) %>% t() # genes in rows now

    SummarizedExperiment::assay(x, "logcounts") <- scaled_dat #replace assay

    return(x)

  })

  return(pb_red)

}
