#' @title Adds a feature.
#' @description Adds a feature created using \code{feat_fun} from \code{fg} OR
#'  \code{m} into a given flowGraph object. Only use this function if
#'  you cannot generate the desired features using the existing flowGraph
#'  functions starting with \code{fg_feat_<feature name>}.
#' @param fg flowGraph object.
#' @param type A string specifying the type of the feature being
#'  added i.e. 'node' or 'edge'.
#' @param feature A string indicating the unique name of the feature added.
#' @param m A numeric matrix with feature values; it should contain the
#'  same sample id's on row names as in \code{fg_get_meta(fg)$id}
#'  and node or edge names
#'  as column names (i.e. if \code{m} is a node feature, it would have the same
#'  column names as those in \code{fg_get_graph(fg)$v$phenotype};
#'  if it is an edge
#'  feature, its column names should be the same as
#'  \code{paste0(fg_get_graph(fg)$e$from, '_', fg_get_graph(fg)$e$to)}).
#' @param feat_fun A function that ouputs a feature matrix as in \code{m} given
#'  \code{fg} and other optional parameters.
#' @param overwrite A logical variable indicating whether or not the function
#'  should replace the existing feature with the same name if
#'  one is already in \code{fg}.
#' @param ... Other parameters that would be used as input into \code{feat_fun}.
#' @return flowGraph object.
#' @details \code{fg_add_feature} adds the given new feature matrix to the
#'  given flowGraph object \code{fg} updating slots
#'  \code{feat} and \code{feat_desc}.
#'  See \code{\link[flowGraph]{flowGraph-class}}
#'  slot \code{feat} and \code{feat_desc} for what should be in these slots.
#'  We do not recommend users to directly use this method unless there is
#'  a clear understanding on how the row and column names should be specified.
#'  Instead, we recommend users to use the functions listed in the "See also"
#'  sections prefixed with "fg_feat_".
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=no_cores)
#'  fg_get_feature_desc(fg)
#'
#'  fg <- fg_add_feature(fg, type="node", feature="count_copy",
#'                       m=fg_data_pos30$count)
#'  fg_get_feature_desc(fg)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_feat_node_prop}}
#'  \code{\link[flowGraph]{fg_feat_node_specenr}}
#'  \code{\link[flowGraph]{fg_get_feature}}
#'  \code{\link[flowGraph]{fg_rm_feature}}
#'  \code{\link[flowGraph]{fg_get_feature_desc}}
#' @rdname fg_add_feature
#' @export
fg_add_feature <- function(
    fg, type="node", feature, m=NULL,
    feat_fun=NULL, overwrite=FALSE, ...
) {
    type <- match.arg(type, c("node", "edge"))
    exists_ <- feature %in% base::names(fg@feat[[type]])
    if (exists_) {
        if (!overwrite) {
            message("skipped")
            return(fg)
        }
        f_ind <- base::names(fg@feat[[type]]) != feature
        fg@feat_desc[[type]] <- fg@feat_desc[[type]][f_ind,,drop=FALSE]
        fg@feat[[type]] <- fg@feat[[type]][f_ind]
    }

    if (base::is.null(m)) {
        if (base::is.null(feat_fun))
            stop("please provide a feature matrix or a function to create one")
        m <- feat_fun(fg, ...)
    }

    fg@feat_desc[[type]] <-
        rbind(fg@feat_desc[[type]], summary_table(m, feature))
    fg@feat[[type]][[feature]] <- m

    return(fg)
}


#' @title Removes a feature.
#' @description Removes a feature from a given flowGraph object.
#' @param fg flowGraph object.
#' @param type A string specifying the type of the feature being
#'  removed i.e. 'node' or 'edge'.
#' @param feature A string indicating the unique name of the feature removed;
#'  note we cannot remove the 'node' 'count' feature type.
#' @return flowGraph object with specified feature removed.
#' @details \code{fg_rm_feature} removes a specified feature matrix from the
#'  given flowGraph object \code{fg} updating slots
#'  \code{feat} and \code{feat_desc}.
#'  See \code{\link[flowGraph]{flowGraph-class}}
#'  slot \code{feat} and \code{feat_desc} for what should be in these slots.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  no_cores=no_cores)
#'  fg_get_feature_desc(fg)
#'
#'  fg <- fg_rm_feature(fg, type="node", feature="prop")
#'  fg_get_feature_desc(fg)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_add_feature}}
#'  \code{\link[flowGraph]{fg_get_feature}}
#'  \code{\link[flowGraph]{fg_get_feature_desc}}
#'  \code{\link[flowGraph]{fg_rm_summary}}
#' @rdname fg_rm_feature
#' @export
fg_rm_feature <- function(fg, type="node", feature=NULL) {
    type <- match.arg(type, c("node", "edge"))
    if (feature=="count" & type=="node")
        stop("cannot remove the count node feature from a flowGraph object.")
    ft_ind <- which(names(fg@feat[[type]]) == feature)
    fg@feat[[type]][[feature]] <- NULL
    if (base::length(ft_ind)==1) {
        fg@feat_desc[[type]] <- fg@feat_desc[[type]][-ft_ind,, drop=FALSE]
    } else {
        warning("feature not found, nothing was dropped")
    }
    # don't need the drop part, but just in case.
    return(fg)
}


#' @title Adds a feature summary.
#' @description Adds a feature summary into a given flowGraph object.
#'  Only use this function if your summary statistic cannot be calcuated
#'  using the \code{\link[flowGraph]{fg_summary}} function.
#' @param fg flowGraph object.
#' @param type A string indicating feature type the summary was created for;
#'  'node' or 'edge'.
#' @param summary_meta The user must provide \code{type} and
#'  \code{summary_meta}.
#'
#'  \code{summary_meta} is a list containing
#'  \code{feature} (feature name), \code{test_name} (summary statistic name),
#'  \code{class} (class), \code{label1}, and \code{label2} (class labels compared).
#'  See \code{\link[flowGraph]{fg_get_summary_desc}} for details.
#' @param p A list containing summary values; this list contains elements:
#'  \code{values} (a vector containing summary statistics e.g. p-values;
#'   this vector should be named by their associated phenotype or edge name),
#'  \code{test_custom} (a function of the statistical test used), and
#'  \code{adjust_custom} (a function of the p-value correction method used).
#'  This list must contain the \code{values} element.
#' @param summ_fun  A function that ouputs a feature summary matrix
#'  as in \code{p} given \code{fg} and other optional parameters.
#' @param overwrite A logical variable indicating whether or not the function
#'  should replace the existing feature summary with the
#'  same name if one is already in \code{fg}.
#' @param ... Other parameters that would be used as input into \code{summ_fun}.
#' @return flowGraph object.
#' @details \code{fg_add_summary} adds the given feature summary list \code{p}
#'  or the output of the given function \code{summ_fun} to the
#'  given flowGraph object \code{fg} updating slots
#'  \code{summary} and \code{summary_desc}.
#'  See \code{\link[flowGraph]{flowGraph-class}}
#'  slot \code{summary} and \code{summary_desc}
#'  for what should be in these slots. We do not recommend users directly use
#'  this function unless what is required is duly in the above slots is
#'  well understood --- note these slots are used in plotting functions
#'  e.g. \code{\link[flowGraph]{fg_plot}}. We instead recommend users to use
#'  the \code{\link[flowGraph]{fg_summary}} function.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  no_cores=no_cores)
#'
#'  # get samples that we are going to compare
#'  m <- fg_get_feature(fg, type="node", feature="prop")
#'  m1_ <- m[fg_data_pos30$meta$class=="control",,drop=FALSE]
#'  m2_ <- m[fg_data_pos30$meta$class=="exp",,drop=FALSE]
#'
#'  # define test or summary function to conduct comparison
#'  test_custom <- function(x,y)
#'      tryCatch(stats::t.test(x,y)$p.value, error=function(e) 1)
#'  values_p <- sapply(seq_len(ncol(m)), function(j)
#'      test_custom(m1_[,j], m2_[,j]) )
#'  values_p <- p.adjust(values_p , method="BY")
#'
#'  # the user can choose to fill either parameter "p" or "summ_fun",
#'  # the latter of which must output a list with the same elements as "p".
#'  # see documentation for ?flowGraph-class, slot "summary" for
#'  # details on what should be in "p".
#'  p <- list(values=values_p, test_fun=test_custom, adjust_fun="BY")
#'  fg <- fg_add_summary(fg, type="node", summary_meta=list(
#'       feature="prop", test_name="wilcox_BY",
#'       class="class", label1="control", label2="exp"), p=p)
#'
#'  fg_get_summary_desc(fg)
#'
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_summary}}
#'  \code{\link[flowGraph]{fg_get_summary}}
#'  \code{\link[flowGraph]{fg_rm_summary}}
#'  \code{\link[flowGraph]{fg_get_summary_desc}}
#'  \code{\link[flowGraph]{fg_add_feature}}
#' @rdname fg_add_summary
#' @export
#' @importFrom purrr map_lgl compact
fg_add_summary <- function(
    fg, type="node", summary_meta=NULL,
    p=NULL, summ_fun=NULL, overwrite=FALSE, ...
) {
    type <- match.arg(type, c("node", "edge"))
    options(stringsAsFactors=FALSE)
    try({
        index <- flowGraph:::fg_get_summary_index(
            fg, type=type, summary_meta=summary_meta)
        if (!overwrite) {
            message("summary exists, skipped")
            return(fg)
        }
        fg <- fg_rm_summary(fg, type=type, index=index)
    }, silent=TRUE)

    if (base::is.null(p)) {
        if (base::is.null(summ_fun))
            stop("provide summary statistic values or a function to create one")
        p <- summ_fun(fg, type=type, ...)  # list(values, test, adjust)
    }

    # legacy, left here just in case
    if ("m1"%in%names(p)) p[["m1"]] <- NULL
    if ("m2"%in%names(p)) p[["m2"]] <- NULL
    p <- purrr::compact(p)

    sm <- data.frame(matrix(summary_meta, nrow=1))
    colnames(sm) <- c("feat", "test_name","class","label1", "label2")
    if (base::length(fg@summary_desc[[type]])==0) {
        fg@summary_desc[[type]] <- sm
    } else {
        fg@summary_desc[[type]] <- rbind(fg@summary_desc[[type]], sm)
        rownames(fg@summary_desc[[type]]) <- NULL
    }

    if (base::is.null(fg@summary[[type]])) fg@summary[[type]] <- list()
    fg@summary[[type]][[nrow(fg@summary_desc[[type]])]] <- p

    return(fg)
}


#' @title Removes a feature summary.
#' @description Removes a feature summary from a given flowGraph object;
#'  while \code{fg} is required, the user can choose to input parameters
#'  \code{summary_meta}, \code{index}, or all of \code{type},
#'  \code{feat}, \code{test_name}, \code{class}, \code{label1},
#'   and \code{label2}.
#'  See \code{\link[flowGraph]{fg_get_summary_desc}} for details.
#' @param fg flowGraph object.
#' @param type A string indicating feature type the summary was created for;
#'  'node' or 'edge'.
#' @param index The user must provide \code{type} and
#'  additionally, one of \code{summary_meta} or \code{index}.
#'
#'  \code{index} is an integer indicating the row in
#'  \code{fg_get_summary_desc(<flowGraph>)} of the corresponding type and
#'  summary the user would like to retrieve.
#' @param summary_meta The user must provide \code{type} and
#'  additionally, one of \code{summary_meta} or \code{index}.
#'
#'  \code{summary_meta} is a list containing

#'  \code{feat} (feature name), \code{test_name} (summary statistic name),
#'  \code{class} (class), \code{label1}, and \code{label2} (class labels compared).
#'  See \code{\link[flowGraph]{fg_get_summary_desc}} for details.
#' @return flowGraph object.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=no_cores)
#'
#'  fg <- fg_summary(fg, no_cores=no_cores, class="class", label1="control",
#'                   overwrite=FALSE, test_name="wilcox_byLayer", diminish=FALSE,
#'                   node_features=NULL, edge_features=NULL)
#'  fg_get_summary_desc(fg)
#'
#'  fg <- fg_rm_summary(fg, summary_meta=c(
#'      feature="count",test_name="wilcox_byLayer",
#'      class="class", label1="control", label2="exp"))
#'  fg_get_summary_desc(fg)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_get_summary}}
#'  \code{\link[flowGraph]{fg_add_summary}}
#'  \code{\link[flowGraph]{fg_get_summary_desc}}
#'  \code{\link[flowGraph]{fg_rm_feature}}
#' @rdname fg_rm_summary
#' @export
#' @importFrom purrr map_lgl
fg_rm_summary <- function(fg, type="node", index=NULL, summary_meta=NULL) {
    type <- match.arg(type, c("node", "edge"))
    index <- fg_get_summary_index(fg,type=type, index,summary_meta)
    fg@summary[[type]][[index]] <- NULL
    # don't need the drop part, but just in case.
    fg@summary_desc[[type]] <- fg@summary_desc[[type]][-index,, drop=FALSE]
    return(fg)
}

#' @title Clears all featuresin a flowGraph object.
#' @description Returns a flowGraph object with only the \code{count} feature.
#' @param fg flowGraph object.
#' @return flowGraph object with only the \code{count} \code{node} feature.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  no_cores=no_cores)
#'
#'  fg <- fg_clear_features(fg)
#'  fg_get_summary_desc(fg)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#' @rdname fg_clear_features
#' @export
#' @importFrom purrr compact
fg_clear_features <- function(fg) {
    if (length(fg@feat$node)>1)
        fg@feat$node[!base::names(fg@feat$node)%in%"count"] <- NULL
    fg@feat$node <- purrr::compact(fg@feat$node)
    fg@feat$edge <- fg@feat_desc$edge <- list()
    fg@feat_desc$node = fg@feat_desc$node[fg@feat_desc$node$feat=="count",]
    return(fg)
}


#' @title Removes all summary statistics.
#' @description Removes all summary statistics in a flowGraph object;
#'  we recommend doing this to save space.
#' @param fg flowGraph object.
#' @return flowGraph object with an empty \code{summary} slot.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=no_cores, node_features="count")
#'  fg_get_summary_desc(fg)
#'
#'  fg <- fg_clear_summary(fg)
#'  fg_get_summary_desc(fg)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_summary}}
#' @rdname fg_clear_summary
#' @export
fg_clear_summary <- function(fg) {
    fg@summary <- fg@summary_desc <- base::list()
    fg@etc$actualVSexpect <- base::list()
    return(fg)
}


#' @title Clears all features and feature summaries in a flowGraph object.
#' @description Returns a flowGraph object with only the \code{count} feature
#'  and meta data. This function clears all other features and
#'  feature summaries to save space.
#' @param fg flowGraph object.
#' @return flowGraph object with all summary statistics and feature values
#'  removed except for the node count feature.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  no_cores=no_cores)
#'
#'  fg <- fg_extract_raw(fg)
#'  show(fg)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#' @rdname fg_extract_raw
#' @export
fg_extract_raw <- function(fg) {
    fg <- fg_clear_summary(fg)
    fg <- fg_clear_features(fg)
    return(fg)
}

#' @title Replace marker names.
#' @description Replace marker names in a flowGraph object.
#' @param fg flowGraph object.
#' @param markers_new A string vector of new marker names;
#'  if \code{markers_old} is set to \code{NULL},
#'  each marker in \code{markers_new} should correspond to
#'  each marker in the \code{markers} slot of the \code{flowGraph} object.
#' @param markers_old A string vector of old marker names user wants to replace;
#' these marker names corresponding to those
#' in \code{fg_get_markers(fg)} with the same length as \code{markers_new}.
#' If \code{markers_old=NULL}, \code{markers_new} should be the same length as
#' \code{fg_get_markers(fg)}.
#' @return flowGraph object with marker names replaced.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=no_cores)
#'
#'  fg <- fg_gsub_markers(fg, c("Anew", "Bnew", "Cnew", "Dnew"))
#'  fg_get_feature_desc(fg)
#'
#' @rdname fg_gsub_markers
#' @export
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_gsub_ids}}
#' @importFrom purrr map
fg_gsub_markers <- function(fg, markers_new, markers_old=NULL) {

    if (base::is.null(markers_old)) {
        if (base::length(markers_new) != base::length(fg@markers))
            stop("incorrect number of markers\n")

        markers_old <- fg@markers
    } else {
        if (base::length(markers_old) != base::length(markers_new))
            stop("incorrect number of markers\n")
    }
    for (mi in base::seq_len(base::length(markers_old))) {
        fg@graph$v$phenotype <-
            base::gsub(markers_old[mi], markers_new[mi], fg@graph$v$phenotype)
        fg@graph$e$from <-
            base::gsub(markers_old[mi], markers_new[mi], fg@graph$e$from)
        fg@graph$e$to <-
            base::gsub(markers_old[mi], markers_new[mi], fg@graph$e$to)
        fg@graph$e$marker <-
            base::gsub(markers_old[mi], markers_new[mi], fg@graph$e$marker)

        fg@edge_list$parent <-
            purrr::map(fg@edge_list$parent, function(x)
                base::gsub(markers_old[mi], markers_new[mi], x))
        base::names(fg@edge_list$parent) <-
            base::gsub(markers_old[mi], markers_new[mi],
                       base::names(fg@edge_list$parent))
        fg@edge_list$child <-
            purrr::map(fg@edge_list$child, function(x)
                purrr::map(x, function(y)
                    base::gsub(markers_old[mi], markers_new[mi], y) ))
        base::names(fg@edge_list$child) <-
            base::gsub(markers_old[mi], markers_new[mi],
                       base::names(fg@edge_list$child))

    }
    fg@feat$node <- purrr::map(fg@feat$node, function(x) {
        base::colnames(x) <- fg@graph$v$phenotype
        x
    })
    if (base::length(fg@feat$edge) > 0) {
        ecn <- paste0(fg@graph$e$from, "_", fg@graph$e$to)
        fg@feat$edge <- purrr::map(fg@feat$edge, function(x) {
            base::colnames(x) <- ecn
            x
        })
    }

    fg@markers[match(markers_old, fg@markers)] <- markers_new
    return(fg)
}

#' @title Replace sample id's.
#' @description Replace sample id's in a flowGraph object.
#' @param fg flowGraph object.
#' @param ids_new A string vector of new sample id's; if \code{ids_old} is
#' set to \code{NULL}, each id in \code{ids_new} should correspond to
#' each id in \code{fg_get_meta(fg)$id}.
#' @param ids_old A string vector of old sample id's the user wants to replace;
#' these marker names corresponding to those
#' in \code{fg_get_meta(fg)$id} with the same length as \code{ids_new}.
#' If \code{ids_old=NULL}, \code{ids_new} should be the same length as
#' \code{fg_get_meta(fg)$id}.
#' @return flowGraph object with sample id's replaced.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=no_cores)
#'
#'  fg <- fg_gsub_ids(fg, ids_new=paste0(fg_get_meta(fg)$id, "_new"))
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_get_feature_desc}}
#'  \code{\link[flowGraph]{fg_gsub_markers}}
#' @rdname fg_gsub_ids
#' @export
#' @importFrom purrr map
fg_gsub_ids <- function(fg, ids_new, ids_old=NULL) {

    if (base::is.null(ids_old)) {
        if (base::length(ids_new) != base::nrow(fg@meta))
            stop("incorrect number of ids\n")

        ids_old <- fg@meta$id
    } else if (base::length(ids_old) != base::length(ids_new)) {
        stop("incorrect number of ids\n")
    }

    ids_ind <- base::match(ids_old, fg@meta$id)
    fg@meta$id[ids_ind] <- ids_new
    fg@feat$node <- purrr::map(fg@feat$node, function(x) {
        base::rownames(x)[ids_ind] <- ids_new
        x
    })
    if (base::length(fg@feat$edge) > 0) {
        fg@feat$edge <- purrr::map(fg@feat$edge, function(x) {
            base::rownames(x)[ids_ind] <- ids_new
            x
        })
    }

    fg@feat_desc <- fg_get_feature_desc(fg, re_calc=TRUE)
    return(fg)
}


#' @title Merges the samples from two flowGraph objects.
#' @description Merges the samples from two flowGraph objects together;
#'  we recommend removing all summary statistics from the new flowGraph object
#'  as those won't be adjusted: \code{\link[flowGraph]{fg_clear_summary}}.
#' @param fg1 flowGraph object.
#' @param fg2 flowGraph object.
#' @return flowGraph object.
#' @details Appends the samples from \code{fg2} onto those in \code{fg1}.
#'  This function requires that the two flowGraph objects must have the
#'  same phenotypes. Therefore, we recommend users to use,
#'  instead, \code{\link[flowGraph]{fg_merge}}.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg0 <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=no_cores)
#'
#'  fg1 <- fg_extract_samples(fg0, fg_get_meta(fg0)$id[1:5])
#'  fg2 <- fg_extract_samples(fg0, fg_get_meta(fg0)$id[4:7])
#'  fg <- fg_merge_samples(fg1, fg2)
#'  fg_get_feature_desc(fg)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_get_feature_desc}}
#'  \code{\link[flowGraph]{fg_merge}}
#'  \code{\link[flowGraph]{fg_extract_samples}}
#' @rdname fg_merge_samples
#' @export
#' @importFrom purrr map
fg_merge_samples <- function(fg1, fg2) {
    fg <- fg1

    cnames = intersect(colnames(fg1@meta), colnames(fg2@meta))

    id2 <- setdiff(fg2@meta$id, fg1@meta$id)
    meta2 <- fg2@meta[fg2@meta$id%in%id2,cnames,drop=FALSE]
    meta1 <- fg1@meta[,cnames,drop=FALSE]

    fg@meta <- base::rbind(meta1, meta2)
    nfs <- base::intersect(base::names(fg1@feat$node),
                           base::names(fg2@feat$node))
    fg@feat$node=purrr::map(nfs, function(xi) {
        a <- fg1@feat$node[[xi]]
        b <- fg2@feat$node[[xi]][id2, , drop=FALSE]
        abcol <- base::intersect(base::colnames(a), base::colnames(b))
        ab <- base::rbind(a[,base::match(abcol, base::colnames(a)), drop=FALSE],
                          b[base::match(base::setdiff(base::rownames(b),
                                                      base::rownames(a)),
                                        base::rownames(b)),
                            base::match(abcol, base::colnames(b)), drop=FALSE])
    })
    base::names(fg@feat$node) <- nfs
    if (base::length(fg@feat$edge) > 0) {
        efs <- base::intersect(base::names(fg1@feat$edge),
                               base::names(fg2@feat$edge))
        fg@feat$edge <- purrr::map(efs, function(xi) {
            a <- fg1@feat$edge[[xi]]
            b <- fg2@feat$edge[[xi]][id2,, drop=FALSE]
            abcol <- base::intersect(base::colnames(a),
                                     base::colnames(b))
            ab <- base::rbind(
                a[, base::match(abcol, base::colnames(a)), drop=FALSE],
                b[base::match(base::setdiff(base::rownames(b),
                                            base::rownames(a)),
                              base::rownames(b)),
                  base::match(abcol, base::colnames(b)),drop=FALSE])
        })
        base::names(fg@feat$edge) <- efs
    }

    fg@feat_desc <- fg_get_feature_desc(fg, re_calc=TRUE)
    fg <- fg_clear_summary(fg)
    warning("merging two flowGraph objects will clear the summary statistics")
    return(fg)
}


#' @title Extracts a set of samples from a flowGraph object.
#' @description Extracts or removes a specified set of samples from
#'  a flowGraph object.
#' @param fg flowGraph object.
#' @param sample_ids A string vector of sample id's that the user wants to
#'  keep in \code{fg}.
#' @param rm_summary A logical indicating whether or not to clear summary.
#' @return flowGraph object.
#' @details The summaries in \code{fg} will not be modified;
#'  we recommend the user recalculates them.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg0 <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=no_cores)
#'  fg_get_feature_desc(fg0)
#'
#'  fg <- fg_extract_samples(fg0, fg_get_meta(fg0)$id[1:5])
#'  fg_get_feature_desc(fg)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_get_feature_desc}}
#'  \code{\link[flowGraph]{fg_merge}}
#'  \code{\link[flowGraph]{fg_extract_phenotypes}}
#' @rdname fg_extract_samples
#' @export
#' @importFrom purrr map
fg_extract_samples <- function(fg, sample_ids, rm_summary=TRUE) {
    if (!any(sample_ids %in% fg@meta$id))
        stop("please provide valid sample id's; see @meta$id\n")


    id_inds <- base::match(sample_ids, fg@meta$id)

    fg@meta <- fg@meta[id_inds, , drop=FALSE]
    fg@feat$node <- purrr::map(fg@feat$node, function(x)
        x[id_inds, , drop=FALSE])
    if (base::length(fg@feat$edge) > 0)
        fg@feat$edge <- purrr::map(fg@feat$edge, function(x)
            x[id_inds, , drop=FALSE])

    fg@feat_desc <- fg_get_feature_desc(fg, re_calc=TRUE)

    if (!base::is.null(fg@summary$node) | !base::is.null(fg@summary$edge))
        if (rm_summary) {
            warning("subsetting samples mean that summary statistics will no longer be valide, removing summary statistics.")
            fg <- fg_clear_summary(fg)
        }
    return(fg)
}


#' @title Extracts a set of phenotypes from a flowGraph object.
#' @description Extracts or removes a specified set of
#'  phenotypes from a flowGraph object.
#' @param fg flowGraph object.
#' @param phenotypes A string vector of phenotype or
#'  cell population name labels.
#' @return flowGraph object.
#' @details The \code{summary} in \code{fg} will not be modified;
#'  we recommend users recalculate them.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg0 <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=no_cores)
#'  fg_get_feature_desc(fg0)
#'
#'  fg <- fg_extract_phenotypes(fg0, fg_get_graph(fg0)$v$phenotype[1:10])
#'  fg_get_feature_desc(fg)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_get_feature_desc}}
#'  \code{\link[flowGraph]{fg_merge}}
#'  \code{\link[flowGraph]{fg_extract_samples}}
#'  \code{\link[flowGraph]{fg_merge_samples}}
#' @rdname fg_extract_phenotypes
#' @export
#' @importFrom purrr map compact
fg_extract_phenotypes <- function(fg, phenotypes) {

    if (!any(phenotypes %in% fg@graph$v$phenotype))
        stop("please provide valid phenotypes; see @graph$v\n")

    id_inds <- fg@graph$v$phenotype %in% phenotypes
    id_inds_ <- fg@graph$e$from %in% phenotypes &
        fg@graph$e$to %in% phenotypes

    fg@graph$v <- fg@graph$v[id_inds,, drop=FALSE]
    fg@graph$e <- fg@graph$e[id_inds_,, drop=FALSE]
    fg <- fg_set_layout(fg, fg@plot_layout)

    fg@feat$node <- purrr::map(fg@feat$node, function(x)
        x[, id_inds, drop=FALSE])
    if (base::length(fg@feat$edge) > 0)
        fg@feat$edge <- purrr::map(fg@feat$edge, function(x)
            x[, id_inds_, drop=FALSE])

    fg@edge_list$child <-
        fg@edge_list$child[phenotypes[phenotypes %in%
                                          base::names(fg@edge_list$child)]]
    fg@edge_list$child <- purrr::map(fg@edge_list$child, function(x)
        purrr::map(x, function(y) {
            a <- y[y %in% phenotypes]
            if (base::length(a) == 0)
                return(NULL)
            a
        }))

    fg@edge_list$child <- purrr::compact(fg@edge_list$child)
    fg@edge_list$parent <-
        fg@edge_list$parent[phenotypes[phenotypes %in%
                                           base::names(fg@edge_list$parent)]]
    fg@edge_list$parent <- purrr::map(fg@edge_list$parent, function(x) {
        a <- x[x %in% phenotypes]
        if (base::length(a) == 0)
            return(NULL)
        a
    })
    fg@edge_list$parent <- purrr::compact(fg@edge_list$parent)

    if (length(fg@summary)>0) {
        if (length(fg@summary$node)>0)
            fg@summary$node <- purrr::map(fg@summary$node, function(x) {
                x$values <- x$values[id_inds]
                x
            })
        if (length(fg@summary$edge)>0)
            fg@summary$edge <- purrr::map(fg@summary$edge, function(x) {
                x$values <- x$values[id_inds_]
                x
            })
    }

    fg@feat_desc <- fg_get_feature_desc(fg, re_calc=TRUE)
    return(fg)
}


#' @title Merges two flowGraph objects together.
#' @description Merges two flowGraph objects together.
#' @param fg1 flowGraph object.
#' @param fg2 flowGraph object.
#' @param method_sample A string indicating how samples from flowGraph objects
#'  should be merged:
#'  \itemize{
#'    \item{\code{union}: keep all samples from both flowGraph objects;
#'     in this case \code{method_phenotype} must be \code{intersect}.}
#'    \item{\code{intersect}: keep only samples that exist
#'     in both \code{fg1} and \code{fg2}.}
#'    \item{\code{setdiff}: keep only samples that exist
#'     in \code{fg1} and not in \code{fg2}.}
#'    \item{\code{none}: keep all samples in \code{fg1}.}
#'  }
#' @param method_phenotype  A string indicating how phenotypes from
#' flowGraph objects should be merged:
#'  \itemize{
#'   \item{\code{intersect}: keep only phenotypes that exist in both
#'    \code{fg1} and \code{fg2}.}
#'   \item{\code{setdiff}: keep only phenotypes that exist in
#'    \code{fg1} and not in \code{fg2}.}
#'   \item{\code{none}: keep all phenotypes in \code{fg1}.}
#'  }
#' @return flowGraph object.
#' @details \code{fg_merge} is a generic function that merges the samples and
#'  phenotypes of two flowGraph objects.
#'  Note that if \code{method_sample="union"}
#'  then \code{method_phenotype} must be set to "intersect".
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg0 <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=no_cores)
#'
#'  fg1 <- fg_extract_samples(fg0, fg_get_meta(fg0)$id[1:5])
#'  fg2 <- fg_extract_samples(fg0, fg_get_meta(fg0)$id[4:7])
#'  fg <- fg_merge(fg1, fg2, method_sample="intersect",
#'                           method_phenotype="intersect")
#'  fg_get_feature_desc(fg)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_extract_samples}}
#'  \code{\link[flowGraph]{fg_extract_phenotypes}}
#'  \code{\link[flowGraph]{fg_merge_samples}}
#' @rdname fg_merge
#' @export
fg_merge <- function(
    fg1, fg2, method_sample=c("union", "intersect", "setdiff", "none"),
    method_phenotype=c("intersect", "setdiff", "none")
) {
    method_sample <- match.arg(method_sample)
    method_phenotype <- match.arg(method_phenotype)

    if (method_sample == "union" & method_phenotype != "intersect")
        stop("if method_sample=union,
             then we must set method_phenotype=intersect;\n
             otherwise, features become difficult to compare.")
    if (method_sample == "union")
        return(fg_merge_samples(fg1, fg2))
    if (method_sample == "none" & method_phenotype == "none")
        return(fg1)

    sample_id1 <- fg1@meta$id
    sample_id2 <- fg2@meta$id
    sample_id1_int <- base::intersect(sample_id1, sample_id2)
    sample_id1_new <- base::setdiff(sample_id1, sample_id2)
    sample_id1_uni <- base::union(sample_id1, sample_id2)

    phen1 <- fg1@graph$v$phenotype
    phen2 <- fg2@graph$v$phenotype
    phen1_int <- base::intersect(phen1, phen2)
    phen1_new <- base::setdiff(phen1, phen2)
    phen1_uni <- base::union(phen1, phen2)

    if (method_sample == "none")
        sample1 <- fg1

    if (method_sample == "intersect") {
        if (base::length(sample_id1_int) == 0)
            stop("no intersecting samples")
        sample1 <- fg_extract_samples(fg1, sample_id1_int)
    }

    if (method_sample == "setdiff") {
        if (base::length(sample_id1_new) == 0)
            stop("no setdiff samples")
        sample1 <- fg_extract_samples(fg1, sample_id1_new)
    }

    if (method_phenotype == "none")
        return(sample1)

    if (method_phenotype == "intersect") {
        if (base::length(phen1_int) == 0)
            stop("no intersecting phenotypes")
        if (all(phen1_int %in% phen2))
            return(sample1)
        return(fg_extract_phenotypes(sample1, phen1_int))
    }

    if (method_phenotype == "setdiff") {
        if (base::length(phen1_new) == 0)
            stop("no setdiff phenotypes")
        if (all(phen1_new %in% phen1))
            return(sample1)
        return(fg_extract_phenotypes(sample1, phen1_new))
    }

}


#' @title Replaces sample meta.
#' @description Replaces sample meta in a given flowGraph object.
#' @param fg flowGraph object.
#' @param meta A data frame containing meta data; see details in
#'  \code{\link[flowGraph]{flowGraph-class}}.
#' @return A flowGraph object with an updated sample meta.
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_get_meta}}
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=no_cores)
#'  head(fg_get_meta(fg))
#'
#'  new_df <- fg_data_pos30$meta
#'  new_df$id[1] <- "newID"
#'
#'  fg <- fg_replace_meta(fg, new_df)
#'  head(fg_get_meta(fg))
#'
#' @rdname fg_replace_meta
#' @export
#' @importFrom purrr map
fg_replace_meta <- function(fg, meta) {
    if (nrow(meta)!=nrow(fg@meta))
        stop("meta must have same number of rows as the number of samples")
    if (!"id"%in%colnames(meta)) stop("meta must have id column")
    fg@meta <- meta
    fg@feat$node <- purrr::map(fg@feat$node, function(x) {
        base::rownames(x) <- meta$id
        x
    })
    if (base::length(fg@feat$edge) > 0)
        fg@feat$edge <- purrr::map(fg@feat$edge, function(x) {
            base::rownames(x) <- meta$id
            x
        })

    return(fg)
}



