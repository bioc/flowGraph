#' @title Calculates feature summary statistics.
#' @description Calculates feature summary statistics for flowGraph features;
#'  users can choose from a list of statistical significance tests/adjustments
#'  or define custom summary functions.
#'  For special cases, see example in function
#'  \code{\link[flowGraph]{fg_add_summary}} on
#'  how to manually calculate summary statistics without using this function.
#' @param fg flowGraph object.
#' @param no_cores An integer indicating how many cores to parallelize on.
#' @param class A string corresponding to the column name or index
#'  of \code{fg_get_meta(fg)} whose values represent
#'  the class label of each sample.
#' @param label1 A string from the \code{class} column of the
#'     \code{meta} slot indicating one of the labels compared to create
#'     the summary statistic. If you would like to compare all other class
#'     labels against one label, set \code{label2} to \code{NULL}.
#'     If you would like to compare all labels against
#'     all labels, set \code{label1} and \code{label2} to \code{NULL}.
#' @param label2 A string from the \code{class} column of the
#'     \code{meta} slot indicating one of the labels compared to create
#'     the summary statistic.
#' @param class_labels A list of vectors, each containing two strings
#'  represeting labels to compare; this parameter is an alternative to
#'  parameters \code{label1} and \code{label2} that supports multiple label
#'  pairings.
#' @param node_features A string vector indicating which node feature(s)
#'  to perform summary statistics on; set to \code{NULL} or \code{"NONE"}
#'  and the function will perform summary statistics on all or no node features.
#' @param edge_features A string vector indicating which edge feature(s)
#'  to perform summary statistics on; set to \code{NULL} or \code{"NONE"}
#'  and the function will perform summary statistics on all or no edge features.
#' @param test_name A string with the name of the test you are performing.
#' @param diminish A logical variable indicating whether to use
#'  diminishing summary statistics;
#'  if \code{TRUE}, a summary statistic for a node or edge will only be done if
#'  at least one of its parent node or edge is significant. Otherwise, the test
#'  will be performed on all nodes or edges.
#' @param p_thres A double indicating the summary statistic threshold;
#'  if the result of a statistical test is greater than \code{p_thres}, then it
#'  is insignificant.
#' @param p_rate A double; if \code{diminish=TRUE}, then \code{p_rate}
#'  needs to be specified.
#'  to determine whether or not a node or edge's parent is significant,
#'  we use \code{p_thres}. However, the higher the layer on which a node
#'  resides or to which an edge points to,
#'  the less stringent this \code{p_thres} should be.
#'  Therefore, we set \code{p_thres}
#'  as the threshold for the parent node or edge of the last
#'  layer and multiply \code{p_thres}
#'  by \code{p_rate} for each increasing layer e.g. given
#'  default values and 4 layers,
#'  the thresholds for layers 1 through 4 would be .4, .2, .1, and .05.
#' @param test_custom A function or a string indicating the statistical test
#'  to use. If a string is provided, it should be one of
#'  \code{c("t","wilcox","ks","var","chisq")};
#'  these correspond to statistical tests \code{stats::t.test},
#'  \code{stats::wilcox.test}, and so on.
#'  If a function is provided, it should take as input two numeric vectors and
#'  output a numeric variable.
#' @param effect_size A logical variable indicating whether or not to calculate
#'  effect size statistic (cohen's d) for this set of class labels;
#'  later used for plotting.
#' @param adjust0 A logical variable indicating whether or not to calculate the
#'  minimum percentage of values from samples of each class label that falls
#'  within the range of \code{adjust0_lim}. This is only done for SpecEnr values
#'  as p-values become unstable when comparing near 0 values.
#' @param adjust0_lim A vector of two numeric values indicating a range around
#'  0, default set to -0.1 and 0.1.
#' @param btwn A logical variable indicating whether or not to calculate the
#'  \code{btwn} data frame given in the \code{fg_get_summary} function.
#' @param btwn_test_custom Same as \code{test_custom} but for \code{btwn}.
#' @param save_functions A logical variable indicating whether to save test and adjust functions.
#' @param overwrite A logical variable indicating whether to
#'  overwrite the existing summary statistics if it exists.
#' @return flowGraph object containing claculated summary statistics.
#' @details \code{fg_summary} calculates a summary statistic as specified
#'  by the user in parameters \code{test_name}, \code{diminish}
#'  (\code{p_thres}, \code{p_rate}), and \code{test_custom}.
#'  The test is done for a node or edge feature of interest within a given
#'  flowGraph object as specified by parameters \code{node_features},
#'  \code{edge_features}. It then returns information on the summary statistic
#'  inside the same flowGraph object and returns it to the user.
#'  See \code{\link[flowGraph]{flowGraph-class}} slot \code{summary}
#'  for details on the contents.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=no_cores)
#'  fg_get_summary_desc(fg)
#'
#'  fg <- fg_summary(fg, no_cores=no_cores, class="class", label1="control",
#'                   overwrite=FALSE, test_name="t", diminish=FALSE,
#'                   node_features="count", edge_features="NONE")
#'  fg_get_summary_desc(fg)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_clear_summary}}
#' @rdname fg_summary
#' @export
#' @importFrom purrr map_lgl map map_dbl
#' @importFrom stats t.test wilcox.test ks.test var.test chisq.test p.adjust
fg_summary <- function(
    fg, no_cores=1, class="class", label1=NULL, label2=NULL, class_labels=NULL,
    node_features="SpecEnr", edge_features="NONE",
    # name your statistical test / adjustment method;
    # if same name, overwrite=T will overwrite, overwrite=F will return nothing
    test_name="t_diminish",
    # don't test if all parents are insignificant, stricter the lower the layer
    diminish=TRUE,
    p_thres=.05, p_rate=2, # only used if diminish=TRUE
    test_custom="t",

    effect_size=TRUE, adjust0=TRUE, adjust0_lim=c(-.1, .1),
    btwn=TRUE, btwn_test_custom="t", save_functions=FALSE,
    overwrite=FALSE
) {
    # if label1 is NULL, compare all
    # if node/edge_features is NULL, do all

    fg_meta <- fg_get_meta

    if (!class%in%colnames(fg_meta)) {
        print(class)
        stop(paste0("invalid class, choose one from fg_get_meta(fg): did not find ", fg_meta))
    }

    test_name <- gsub("[.]","_",test_name)
    classes_ <- fg_meta[,class]
    classes <- unique(classes_)
    if (length(classes)<2)
        stop("please provide a valid class from columns of fg_get_meta(fg).")

    # class list lists what classes to compare
    if (!is.null(label1))
        if (!label1%in%fg_meta[,class] & !is.null(class_labels)) {
            label1 <- NULL
            stop("we will proceed to compare all class labels with each other.")
        }
    if (!is.null(class_labels)) {
        if (!all(purrr::map_lgl(class_labels, function(x) x%in%classes_)))
            stop("please provide valid class labels class_labels.")
        class_list <- list(class_labels)
    } else if (!is.null(label1) & is.null(label2)) {
        class_list <- purrr::map(classes[classes!=label1],
                               function(x) c(label1, x))
    } else if (!is.null(label1) & !is.null(label2)) {
        class_list <- list(a=c(label1, label2))
    } else {
        class_list <- purrr::map(2:length(classes), function(i)
            purrr::map(seq_len(i-1), function(j) c(classes[i], classes[j])))
        class_list <- unlist(class_list, recursive=FALSE)
    }

    # test functions
    if (!is.function(test_custom)) {
        if (test_custom=="t") {
            test_custom <- function(x,y)
                tryCatch(stats::t.test(x,y)$p.value, error=function(e) 1)
        } else if (test_custom=="wilcox") {
            test_custom <- function(x,y)
                tryCatch(stats::wilcox.test(x,y)$p.value, error=function(e) 1)
        } else if (test_custom=="ks") {
            test_custom <- function(x,y)
                tryCatch(stats::ks.test(x,y)$p.value, error=function(e) 1)
        } else if (test_custom=="var") {
            test_custom <- function(x,y)
                tryCatch(stats::var.test(x,y)$p.value, error=function(e) 1)
        } else if (test_custom=="chisq") {
            test_custom <- function(x,y)
                tryCatch(stats::chisq.test(x,y)$p.value, error=function(e) 1)
        }
    }

    fg_feat <- fg_get_feature_all(fg)
    for (type in c("node","edge")) {
        if (type=="node")
            features <- node_features
        if (type=="edge")
            features <- edge_features

        if (is.null(features)) features <- names(fg_feat[[type]])
        if (length(features)==0) next
        if (features[1]=="NONE") next
        features <- features[features%in%names(fg_feat[[type]])]
        if (length(features)==0) next

        for (feature in features) {
            if (!feature%in%names(fg_feat[[type]])) next
            for (classl in class_list) {
                start1 <- Sys.time()
                message("- claculating summary statistics (for class ",
                        class, ": ", paste0(classl, collapse=" & "), ") ",
                        test_name," for ", type, " feature ", feature)

                id1 <- classes_==classl[1]
                id2 <- classes_==classl[2]
                m <- fg_get_feature(fg, type=type, feature=feature)

                # calculate summary
                fg <- fg_add_summary(fg, type=type, summary_meta=list(
                    feature=feature, test_name=test_name,
                    class=class, label1=classl[1], label2=classl[2]),
                    overwrite=overwrite, p=NULL, summ_fun=fg_summary_,
                    m1=m[id1,,drop=FALSE], m2=m[id2,,drop=FALSE],
                    diminish=diminish, test_custom=test_custom,
                    p_thres=p_thres, p_rate=p_rate,
                    save_functions=save_functions, no_cores=no_cores)


                sm_name <- paste0(c(feature, class, classl), collapse="_")

                if (effect_size) {
                    if (is.null(fg_get_etc(fg)$effect_size))
                        fg@etc$effect_size <- list()
                    if (is.null(fg_get_etc(fg)$effect_size[[type]]))
                        fg@etc$effect_size[[type]] <- list()
                    if (is.null(fg_get_etc(fg)$effect_size[[type]][[sm_name]])) {
                        fg@etc$effect_size[[type]][[sm_name]] <-
                            fg_cohensd_(fg_get_feature(fg, type, feature), id1, id2)
                    }
                }

                if (adjust0 & grepl("SpecEnr", feature)) {
                    if (is.null(fg@etc$adjust0))
                        fg@etc$adjust0 <- list()
                    if (is.null(fg@etc$adjust0[[type]]))
                        fg@etc$adjust0[[type]] <- list()
                    if (is.null(fg@etc$adjust0[[type]][[sm_name]])) {
                        fg@etc$adjust0[[type]][[sm_name]] <-
                            fg_adjust0_(fg@feat[[type]][[feature]], id1, id2,
                                        adjust0_lim=adjust0_lim)
                    }
                }

                if (btwn & grepl("SpecEnr", feature)) {
                    # test functions
                    if (!is.function(btwn_test_custom)) {
                        if (btwn_test_custom=="t") {
                            btwn_test_custom <- function(x,y)
                                tryCatch(stats::t.test(x,y)$p.value,
                                         error=function(e) 1)
                        } else if (btwn_test_custom=="wilcox") {
                            btwn_test_custom <- function(x,y)
                                tryCatch(stats::wilcox.test(x,y)$p.value,
                                         error=function(e) 1)
                        } else if (btwn_test_custom=="ks") {
                            btwn_test_custom <- function(x,y)
                                tryCatch(stats::ks.test(x,y)$p.value,
                                         error=function(e) 1)
                        } else if (btwn_test_custom=="var") {
                            btwn_test_custom <- function(x,y)
                                tryCatch(stats::var.test(x,y)$p.value,
                                         error=function(e) 1)
                        } else if (btwn_test_custom=="chisq") {
                            btwn_test_custom <- function(x,y)
                                tryCatch(stats::chisq.test(x,y)$p.value,
                                         error=function(e) 1)
                        }
                    }

                    cf <- se_feats(feature)
                    if (is.null(fg@etc$actualVSexpect))
                        fg@etc$actualVSexpect <- list()
                    if (is.null(fg@etc$actualVSexpect[[type]]))
                        fg@etc$actualVSexpect[[type]] <- list()
                    if (is.null(fg@etc$actualVSexpect[[type]][[sm_name]])) {
                        mp <- fg_get_feature(fg, type=type, feature=cf[2])
                        mep <- fg_get_feature(fg, type=type, feature=cf[3])
                        fg@etc$actualVSexpect[[type]][[sm_name]] <-
                            fg_btwn_(id1, id2, mp, mep,
                                     btwn_test_custom=btwn_test_custom)
                    }
                }
                time_output(start1)
            }
        }
    }
    return(fg)
}


# calculates relation between expected and actual raw feature values
fg_btwn_ <- function(id1, id2, mp, mep, btwn_test_custom) {
    id1_ <- id1
    id2_ <- id2
    if (sum(id1_)>sum(id2_)) {
        id1 <- sample(which(id1_),sum(id2_))
    } else if (sum(id2_)>sum(id1_)) {
        id2 <- sample(which(id2_),sum(id1_))
    }

    # whether specenr should be significant based on how different actual
    # and expected raw values are; sigcands are those that are significantly different

    pp3v1 <- sapply(seq_len(ncol(mp)), function(ci)
        btwn_test_custom(mp[id1,ci], mep[id1,ci]) )
    pp3v2 <- sapply(seq_len(ncol(mp)), function(ci)
        btwn_test_custom(mp[id2,ci], mep[id2,ci]) )

    pp3v1[is.na(pp3v1)] <- pp3v2[is.na(pp3v2)] <- 1

    pp3c1 <- sapply(seq_len(ncol(mp)), function(ci)
        effsize::cohen.d(mp[id1,ci], mep[id1,ci])$estimate )
    pp3c2 <- sapply(seq_len(ncol(mp)), function(ci)
        effsize::cohen.d(mp[id2,ci], mep[id2,ci])$estimate )
    pp3c1[is.na(pp3c1)] <- pp3c2[is.na(pp3c2)] <- 0
    # |d|<0.2 "negligible", |d|<0.5 "small", |d|<0.8 "medium", otherwise "large"

    pp3b_raw <- lapply(seq_len(ncol(mp)), function(ci) {
        a <- mp[id1_,ci] - mp[id2_,ci]
        b <- mep[id1_,ci] - mep[id2_,ci]
        return(list(a=a,b=b))
    })
    pp3bt <- sapply(seq_len(ncol(mp)), function(ci)
        btwn_test_custom(
            pp3b_raw[[ci]]$a, pp3b_raw[[ci]]$b) )
    pp3bt[is.na(pp3bt)] <- 1
    pp3bc <- sapply(seq_len(ncol(mp)), function(ci)
        effsize::cohen.d( pp3b_raw[[ci]]$a, pp3b_raw[[ci]]$b)$estimate )
    pp3bc[is.na(pp3bc)] <- 0

    pp3b_raw_ <- lapply(seq_len(ncol(mp)), function(ci) {
        a <- log(mp[id1_,ci] / mp[id2_,ci])
        b <- log(mep[id1_,ci] / mep[id2_,ci])
        return(list(a=a,b=b))
    })
    pp3bt_ <- sapply(seq_len(ncol(mp)), function(ci)
        btwn_test_custom(
            pp3b_raw_[[ci]]$a, pp3b_raw_[[ci]]$b) )
    pp3bt_[is.na(pp3bt_)] <- 1
    pp3bc_ <- sapply(seq_len(ncol(mp)), function(ci)
        effsize::cohen.d( pp3b_raw_[[ci]]$a, pp3b_raw_[[ci]]$b)$estimate )
    pp3bc_[is.na(pp3bc_)] <- 0

    df_ps <- data.frame(
        phenotype=colnames(mp),
        tpv1=pp3v1, tpv2=pp3v2,
        cd1=pp3c1, cd2=pp3c2,
        btp=pp3bt, bcd=pp3bc,
        btp_=pp3bt_, bcd_=pp3bc_
    )
    return(df_ps)
}

# calculates min percentage of 0's between two groups of samples for each SpecEnr value
fg_adjust0_ <- function(m, id1, id2, adjust0_lim=c(-.1, .1)) {
    m <- as.matrix(m)
    xl <- min(sum(id1), length(id1))
    yl <- min(sum(id2), length(id2))
    purrr::map_dbl(seq_len(ncol(m)), function(ci)
        min(sum(m[id1,ci]<adjust0_lim[2] & m[id1,ci]>adjust0_lim[1])/xl,
            sum(m[id2,ci]<adjust0_lim[2] & m[id2,ci]>adjust0_lim[1])/yl))
}

fg_med0_ <- function(m, id1, id2) {
    m <- as.matrix(m)
    purrr::map_dbl(seq_len(ncol(m)), function(ci)
        max(abs(c(median(m[id1,ci]), stats::median(m[id2,ci])))))
}

# calculates effect size statistics cohens d
#' @importFrom purrr map map_dbl map_int
#' @importFrom effsize cohen.d
fg_cohensd_ <- function(m, id1, id2) {
    cd <- purrr::map(seq_len(ncol(m)), function(ci)
        effsize::cohen.d(m[id1,ci], m[id2,ci]))
    cde <- purrr::map_dbl(cd, function(x) x$estimate)
    cde[is.na(cde)] <- 0
    cdd <- purrr::map_int(cd, function(x) x$magnitude)
    cdd[is.na(cdd)] <- 1
    cdd <- factor(cdd)
    levels(cdd) <- levels(cd[[which(!is.na(cdd))[1]]]$magnitude)
    return(list(cohensd=cde, cohensd_size=cdd))
}

# calculates summary statistics
#' @importFrom future plan multiprocess
#' @importFrom furrr future_map
#' @importFrom purrr map
fg_summary_ <- function(
    fg, m1, m2, type, diminish,
    test_custom, p_thres, p_rate, save_functions,
    no_cores=1
) {
    if (no_cores>1) future::plan(future::multiprocess)

    # split feature matrix into sample classes
    mnames <- colnames(m1)
    ml <- length(mnames)

    fg_graph <- fg_get_graph(fg)

    # prepare layers
    if (type=="node") layers_ <- fg_graph$v$phenolayer
    if (type=="edge") {
        phen_to <- fg_graph$e$to
        phen_from <- fg_graph$e$from
        layers_ <- fg_graph$v$phenolayer[
            match(phen_to, fg_graph$v$phenotype)]
    }
    layers <- sort(unique(layers_))

    if (!diminish) {
        # calculate p values
        if (no_cores>1) {
            loop_ind <- flowGraph:::loop_ind_f(seq_len(ml), no_cores)
            p <- unlist(furrr::future_map(loop_ind, function(ii)
                purrr::map_dbl(ii, function(i) test_custom(m1[,i], m2[,i])) ))
        } else {
            p <- purrr::map_dbl(seq_len(ml), function(i)
                test_custom(m1[,i], m2[,i]))

        }

    } else {
        pchild <- fg_get_el(fg)$child

        # p value thresholds determining whether to continue down hierarchy
        p_thress <- rep(p_rate, length(layers))
        if (any(layers==0)) p_thress <- p_thress[-length(p_thress)]
        p_thress[1] <- p_thres
        p_thress <- cumprod(p_thress)
        p_thress <- sort(p_thress, decreasing=TRUE)

        # initialize p values vector
        p <- rep(NA, ml)
        root_ <- mnames==""
        if (sum(root_)>0) {
            p[root_] <- test_custom(m1[,root_], m2[,root_])
            if (is.nan(p[root_])) p[root_] <- 1
        }

        # calculate p values for layer 1
        testis <- testi <- layers_==min(layers_[!root_])
        p_thress_i <- 1
        while (sum(testi)>0) {
            if (no_cores>1) {
                loop_ind <- loop_ind_f(which(testi), no_cores)
                pl <- unlist(furrr::future_map(loop_ind, function(ii)
                    purrr::map_dbl(ii, function(i)
                        pt <- test_custom(m1[,i], m2[,i])
                    )
                ))
            } else {
                pl <- purrr::map_dbl(which(testi), function(i)
                        pt <- test_custom(m1[,i], m2[,i]))
            }
            p[testi] <- pl

            pl_cpops <- pl < p_thress[p_thress_i]
            pl_cpops_ <- mnames[which(testi)[pl_cpops]]
            p_thress_i <- p_thress_i + 1
            if (type=="node") {
                testi <- mnames %in% unlist(pchild[pl_cpops_])
            } else {
                good_edges <- mnames %in% pl_cpops_
                testi <- phen_from %in% phen_to[good_edges]
            }
            testis <- testis | testi
        }
    }

    names(p) <- mnames
    p_ <- list(values=p)
    if (save_functions)
        p_$test_fun <- test_custom

    return(p_)
}





