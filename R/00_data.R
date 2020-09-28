#' fg_data_fca
#'
#' @format A list containing the following elements
#'  derived from the flowCAP-II AML data set for cell populations up to layer 3.
#' \itemize{
#'   \item{\code{count}: A numeric sample x cell population node
#'   matrix with cell count values.}
#'   \item{\code{meta}: A data frame containing meta information on
#'    samples in \code{count}; it contains columns:}
#'   \itemize{
#'       \item{\code{class}: a string indicating whether a sample is
#'        from a "control" or "aml" subject.}
#'       \item{\code{id}: a string containing sample id's.}
#'       \item{\code{train}: a logical variable indicating whether
#'        a sample is from the train or test set.}
#'       \item{\code{subject}: a numeric variable containing the
#'        id of the subject from whom the sample came from.}
#'       \item{\code{tube}: the tube or panel number;
#'        all samples in this data set is
#'        analyzed under the 6th panel.}
#'   }
#' }
#' @source {
#'   \insertRef{aghaeepour2013critical}{flowGraph}
#' }
"fg_data_fca"

#' fg_data_pos15
#'
#' @format A list containing the following elements
#'  for a positive control data set with markers A, B, C, D. This is a positive
#'  control data set where node A+B+C+ increased by 50%.
#' \itemize{
#'   \item{\code{count}: A numeric sample x cell population node
#'   matrix with cell count values}
#'   \item{\code{meta}: A data frame containing meta information on
#'    samples in \code{count}; it contains columns:}
#'   \itemize{
#'       \item{\code{id}: a string containing sample id's.}
#'       \item{\code{class}: a string indicating whether a sample is
#'        from a "control" or "exp" (experiment) subject.}
#'   }
#' }
"fg_data_pos15"


#' fg_data_pos30
#'
#' @format A list containing the following elements
#'  for a positive control data set with
#'  markers A, B, C, D; note it was made with two and three thresholds
#'  for markers A and B to test functions with multiple thresholds
#'  (this is a positive
#'  control data set where nodes A+..B+..C+ increased by 50%).
#' \itemize{
#'   \item{\code{count}: A numeric sample x cell population node
#'   matrix with cell count values}
#'   \item{\code{meta}: A data frame containing meta information on
#'   samples in \code{count}; it contains columns:}
#'   \itemize{
#'       \item{\code{id}: a string containing sample id's.}
#'       \item{\code{class}: a string indicating whether a sample is
#'        from a "control" or "exp" (experiment) subject.}
#'   }
#' }
"fg_data_pos30"


# # load data
# .onLoad <- function(libname, pkgname) {
#     data("fg_data_fca", "fg_data_pos15", "fg_data_pos30",
#          package=pkgname, envir=paresave(nt.env(environment()))
# }
