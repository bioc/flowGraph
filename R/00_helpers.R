## HELPERS ---------------------------------------------------------------------



#' #' @title Prepares number of cores for parallel backend.
#' #' @description Prepares number of cores needed for parallel backend.
#' #' @param no_cores An integer indicating how many cores to parallelize on.
#' #' @return An integer indicating how many cores to parallelize on.
#' #' @details Given the number of cores a user wishes to parallelize processes on
#' #'  \code{no_cores}, \code{ncores} ensures this value does not
#' #'  exceed the actual number of cores the user's computer contains.
#' #' @examples
#' #'
#' #'  # NOT EXPORTED
#' #'  no_cores <- 100
#' #'  flowGraph:::ncores(no_cores)
#' #'
#' #' @seealso
#' #'  \code{\link[parallel]{detectCores}}
#' #' @rdname ncores
#' #' @importFrom parallel detectCores
#' ncores <- function(no_cores=1)
#'     max(1, min(no_cores, parallel::detectCores()-1))


#' @title Formats time into string.
#' @description Formats time into a string HH:MM:SS given time zone.
#' @param time A time variable of class \code{POSIXct}, \code{POSIXt}.
#' @return Time formatted as a string; used in \code{time_output} function.
#' @examples
#'  # NOT EXPORTED
#'  flowGraph:::tstr(Sys.time())
#'
#' @rdname tstr
tstr <- function(time) format(.POSIXct(time), "%H:%M:%S")


#' @title Outputs elapsed time.
#' @description Given a time, prints the time elapsed from that time until now.
#' @param start A time variable of class \code{POSIXct}, \code{POSIXt}.
#' @param msg A string with a message to print out after the elapsed time.
#' @return Prints to console, the time from which process
#'  started \code{start} - ended, and > time elapsed from
#'  \code{start} until now.
#' @examples
#'
#'  start <- Sys.time()
#'  flowGraph:::time_output(start,'start - now > time elapsed')
#'
#' @rdname time_output
# #' @export
time_output <- function(start, msg="") {
    start <- as.POSIXct(start)
    end <- Sys.time()
    time_elapsed <- difftime(end, start, units="secs")
    message(msg, ifelse(msg == "", "", ": "),
            tstr(start), "-", tstr(end), " > ", tstr(time_elapsed))
}


#' @title Prepares parallel loop indices.
#' @description \code{loop_ind_f} is a helper function that splits
#'  a vector of loop indices into a list of multiple loop indices
#'  for use in parallel processes within the flowGraph package.
#' @param x A vector of loop indices.
#' @param n An integer, or the number of vectors to split \code{x} into.
#' @return list of \code{n} vectors with elements from \code{x}.
#' @examples
#'
#'  old_loop_inds <- 1:10
#'  no_cores <- 5
#'
#'  new_loop_inds <- flowGraph:::loop_ind_f(old_loop_inds, no_cores)
#'  # future::plan(future::multiprocess)
#'  # example_indices <- furrr::future_map(new_loop_inds, function(ii) {
#'  #     purrr::map(ii, function(i) i )
#'  # s})
#'
#' @rdname loop_ind_f
# #' @export
loop_ind_f <- function(x, n) {
    if (n == 1) return(base::list(x))
    return(base::split(x, ceiling(seq_along(x)/ceiling(base::length(x)/n))))
}


#' @title Summarizes a numeric matrix.
#' @description Summarizes a numeric matrix.
#' @param m A numeric matrix.
#' @param feat_type Name of the matrix \code{m}.
#' @return A data frame containing one row summarizing \code{m};
#'  see \code{\link[flowGraph]{fg_get_feature_desc}}.
#' @examples
#'
#'  summary_table(matrix(rnorm(12),nrow=3), feat_type='random')
#'
#' @rdname summary_table
#' @export
summary_table <- function(m, feat_type="") {
    m <- as.matrix(m)
    base::data.frame(feat=feat_type, nrow=base::nrow(m),
                     ncol=base::ncol(m), inf=sum(is.infinite(m)),
                     neginf=sum(m == -Inf), na=sum(base::is.na(m)),
                     nan=sum(is.nan(m)), neg=sum(m < 0),
                     pos=sum(m > 0), zero=sum(m == 0), max=max(m[is.finite(m)]),
                     min=min(m[is.finite(m)]))
}


#' @title Normalizes matrix values by class.
#' @description Used only in the \code{\link[flowGraph]{fg_feat_mean_class}}
#'  function; for each class in the \code{classes} vector,
#'  \code{meandiff} takes the column mean
#'  of the rows in the given matrix associated with that class;
#'  it then takes the difference point by point between these means and
#'  the original rows for that class.
#' @param m0 A numeric matrix.
#' @param classes A vector whose length is equal to the number of
#'  rows in the given matrix.
#' @return A numeric matrix whose dimensions equate to that of the input
#'  and whose values are normalized per class.
#' @examples
#'
#'  classes <- append(rep('apples',4), rep('oranges',3))
#'  m0 <- matrix(rnorm(35), nrow=7)
#'  m <- flowGraph:::mean_diff(m0, classes)
#'
#' @seealso
#'  \code{\link[flowGraph]{fg_feat_mean_class}}
#' @rdname meandiff
# #' @export
#' @importFrom purrr map
mean_diff <- function(m0, classes) {
    m <- m0
    for (pi in base::unique(classes)) {
        pii <- classes == pi
        m_ <- m[pii, , drop=FALSE]
        m_m <- base::colMeans(m_)
        m[pii, ] <- base::do.call(rbind,
                                  purrr::map(base::seq_len(base::nrow(m_)),
                                             function(i) m_[i, ] - m_m))
    }
    base::dimnames(m) <- base::dimnames(m0)
    return(m)
}



