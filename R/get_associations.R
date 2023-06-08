#' Associate factors to covariates of interest
#'
#' @description
#' Performs Analysis of Variance (ANOVA) or linear models to associate
#' factor scores with covariates of interest provided by the user
#'
#' @details
#' Given a covariate of interest and a defined test, this function tests for associations with factor scores. For
#' categorical tests, ANOVAs are fitted, while for continous variables, linear models. P-values are corrected
#' using the Benjamini-Hochberg procedure.
#'
#' @param model A MOFA2 model.
#' @param metadata A data frame containing the annotations of the samples included in the MOFA model.
#' @param sample_id_column A string character that refers to the column in `metadata` where the sample identifier is located.
#' @param test_variable A string character that refers to the column in `metadata` where the covariate to be tested is located.
#' @param test_type A string character ("categorical", "continuous").
#' @param categorical_type A string character ("parametric", "not_parametric"), only applies for categorical data.
#' @param group Boolean flag TRUE/FALSE, to specify if a grouped MOFA model is provided.
#'
#'
#' @return A dataframe in a tidy format containing the p-values of the association tests per factor
#' @export
#'
#' @import MOFA2
#' @import tibble
#' @import dplyr
#' @import broom
#' @import stats
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
#'                                       metadata = metadata,
#'                                       sample_id_column = "sample",
#'                                       test_variable = "var",
#'                                       test_type = "continuous",
#'                                       group = FALSE)
get_associations <- function(model,
                             metadata,
                             sample_id_column,
                             test_variable,
                             test_type = "categorical",
                             categorical_type = "parametric",
                             group = FALSE) {

  factors <- get_tidy_factors(model = model,
                              metadata = metadata,
                              factor = "all",
                              group = group,
                              sample_id_column = sample_id_column)

  # Get factors associated with patient group
  factors <- factors %>%
    dplyr::select_at(c("sample", test_variable, "Factor", "value")) %>%
    na.omit() %>% # this allows to make subsets of data
    dplyr::group_by(Factor) %>%
    tidyr::nest() %>%
    dplyr::mutate(pvalue = purrr::map(.data$data, function(dat) {

      # Fit ANOVAs if testing variable is categorical
      if(test_type == "categorical") {

        if(categorical_type == "parametric") {

          factor_aov <- stats::aov(as.formula(paste0("value ~ ", test_variable)), data = dat) %>%
            broom::tidy() %>%
            dplyr::filter(term == test_variable) %>%
            dplyr::select_at(c("term", "p.value"))

          return(factor_aov)

        } else {

          factor_kw <- stats::kruskal.test(as.formula(paste0("value ~ ", test_variable)), data = dat) %>%
            broom::tidy() %>%
            dplyr::mutate(term = test_variable) %>%
            dplyr::select_at(c("term", "p.value"))

          return(factor_kw)

        }

      } else {
        # Fit linear model if testing variable is continous
        factor_lm <- stats::lm(as.formula(paste0("value ~ ", test_variable)), data = dat) %>%
          broom::tidy() %>%
          dplyr::filter(term == test_variable) %>%
          dplyr::select_at(c("term", "p.value"))

        return(factor_lm)

      }

    })) %>%
    tidyr::unnest(pvalue) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(adj_pvalue = stats::p.adjust(p.value))

  expl_var <- factors %>%
    dplyr::select_at(c("Factor", "term", "p.value", "adj_pvalue"))

  return(expl_var)

}
