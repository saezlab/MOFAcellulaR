#' Extract factor scores from a model in a tidy format with meta data
#'
#' @description
#' Generates a tidy dataframe for factors of interest useful for plotting functions.
#'
#' @details
#'
#' This function simplifies the extraction of factor scores from MOFA models.
#'
#' @param model A MOFA2 model.
#' @param metadata A data frame containing the annotations of the samples included in the MOFA model.
#' @param factor A string character to define the factor to be extracted. Alternatively, "all" to get all factors.
#' @param group Boolean flag TRUE/FALSE, to specify if a grouped MOFA model is provided.
#' @param sample_id_column A string character that refers to the column in `metadata` where the sample identifier is located.
#'
#'
#' @return A dataframe in a tidy format containing the factor scores per sample together with metadata
#' @export
#'
#' @import MOFA2
#' @import tibble
#' @import dplyr
#'
#' @examples
#' inputs_dir <- base::system.file("extdata", package = "MOFAcellulaR")
#' model <- MOFA2::load_model(file.path(inputs_dir, "testmodel.hdf5"))
#' metadata <- readRDS(file.path(inputs_dir, "testmetadata.rds"))
#' all_factors <- get_tidy_factors(model = model,
#'                                 metadata = metadata,
#'                                 factor = "all",
#'                                 sample_id_column = "sample")
#' Factor3 <- get_tidy_factors(model = model,
#'                             metadata = metadata,
#'                             factor = "Factor3",
#'                             sample_id_column = "sample")
get_tidy_factors <-  function(model,
                              metadata,
                              factor,
                              group = FALSE,
                              sample_id_column = "sample") {

  if(group) {

    if(factor == "all") {
      #Here we return all factors and bind the distinct groups
      factor_scores <- MOFA2::get_factors(model, factors = "all") %>%
        base::do.call(base::rbind, .)

    } else {
      #Here we return specific factors
      factor_scores <- MOFA2::get_factors(model, factors = "all") %>%
        base::do.call(base::rbind, .)

      factor_scores <- factor_scores[, factor, drop= F]

    }

  } else {

    if(factor == "all") {
      #Here we return all factors
      factor_scores <- MOFA2::get_factors(model, factors = "all")[[1]]
    } else {
      #Here we return specific factors
      factor_scores <- MOFA2::get_factors(model, factors = "all")[[1]][,factor, drop= F]
    }

  }

  # Merge with provided meta_data
  factor_scores <- factor_scores %>%
    base::as.data.frame() %>%
    tibble::rownames_to_column(sample_id_column) %>%
    dplyr::left_join(metadata,
              by = sample_id_column)  %>%
    tidyr::pivot_longer(-colnames(metadata),
                        names_to = "Factor") %>%
    dplyr::rename("sample" = sample_id_column)

  return(factor_scores)

}

#' Extract feature weights scores from a model in a tidy format for a factor of interest
#'
#' @description
#' Generates a tidy dataframe of feature weights for factors of interest useful for plotting function and enrichment analysis
#'
#' @details
#'
#' This function simplifies the extraction of factor feature weights from MOFA models.
#'
#' @param model A MOFA2 model.
#' @param factor A string character to define the factor to be extracted.
#'
#'
#' @return A dataframe in a tidy format containing the factor feature weights
#' @export
#'
#' @import MOFA2
#' @import dplyr
#' @import purrr
#'
#' @examples
#' inputs_dir <- base::system.file("extdata", package = "MOFAcellulaR")
#' model <- MOFA2::load_model(file.path(inputs_dir, "testmodel.hdf5"))
#' gene_weights <- get_geneweights(model = model,
#'                                 factor = "Factor1")
get_geneweights <- function(model, factor) {

  factor_loadings <- MOFA2::get_weights(model, as.data.frame = T) %>%
    base::as.data.frame() %>%
    dplyr::mutate(feature = purrr::map2_chr(view, feature, function(v, f) {

      gsub(paste0(v, "_"), "",f)

    })) %>%
    dplyr::rename("ctype" = view) %>%
    dplyr::rename("factors" = factor) %>%
    dplyr::filter(factors %in% factor) %>%
    dplyr::select(-factors)

  return(factor_loadings)

}
