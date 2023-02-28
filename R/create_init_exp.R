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
#' @examples
#' x <- matrix(rep(1,4),nrow = 2, ncol = 2)
#' rownames(x) <- c("featA", "featB")
#' colnames(x) <- c("Sample1", "Sample2")
#' y <- data.frame(cell_type = c("A","A"),
#' donor_id = c("Sample1", "Sample2"),
#' cell_counts = c(30,30))
#' se <- create_init_exp(x, y)
create_init_exp <- function(counts, coldata) {

  pb_dat <- SummarizedExperiment::SummarizedExperiment(assays = list("counts" = counts),
                                                       colData = S4Vectors::DataFrame(coldata))

  return(pb_dat)
}
