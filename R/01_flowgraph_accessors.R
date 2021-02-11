#' @title Saves flowGraph object to a specified path.
#' @description Saves flowGraph object to a specified path.
#' @param fg flowGraph object to save.
#' @param folder_path A string indicating the folder path to where the flowGraph
#'  object should save its elements; if this is the first time the object is
#'  being saved, this folder should be empty or if it is
#'  not yet created, the function will create it. If the object has previously
#'  been saved before and this parameter is set to \code{NULL}, the function
#'  will save the object into the save folder it was previously saved in.
#' @param save_plots A logical indicating whether or not to save plots.
#' @param paired A logical indicating whether the summary is paired; used in
#'  function \code{fg_plot_box}.
#' @param ... Other parameters for the \code{fg_save_plots} function.
#' @details See generated README.md file.
#' @return TRUE if flowGraph object successfully saved.
#' @examples
#'  no_cores <- 1
#'  data(fg_data_pos2)
#'  fg <- flowGraph(fg_data_pos2$count, class=fg_data_pos2$meta$class,
#'                  no_cores=no_cores)
#'
#'  fg_save(fg, "tmp")
#'
#' @seealso
#'  \code{length},\code{c("nrow", "nrow")},\code{NULL}
#'  \code{\link[purrr]{map}}
#' @rdname fg_save
#' @export
#' @importFrom purrr map
#' @importFrom data.table fwrite
#' @importFrom utils write.csv
fg_save <- function(fg, folder_path=NULL, save_plots=TRUE, paired=FALSE, ...) {
    options(encoding="UTF-8")

    # check folder path
    if (is.null(folder_path)) {
        folder_path <- fg_get_etc(fg)$save$path
        if (is.null(folder_path))
            stop("please provide valid folder path")
    }
    if (!dir.exists(folder_path))
        dir.create(folder_path, recursive=TRUE, showWarnings=FALSE)
    if (length(dir(folder_path, all.files=TRUE))==0) {
        stop_save <- TRUE

        etc_dir <- paste0(folder_path,"/etc/etc.rds")
        if (file.exists(etc_dir))
            if (readRDS(etc_dir)$save$id == fg_get_etc(fg)$save$id)
                stop_save <- FALSE

        if (stop_save)
            stop("the folder provided is non-empty")
    }

    # create readme file
    fg_graph <- fg_get_graph(fg)
    fc <- file(paste0(folder_path,"/README.md"))
    writeLines(c(
        "PLEASE DO NOT MAKE CHANGES TO THIS FILE.",
        paste0("Last saved: ",  Sys.time()),
        "",
        "This folder contains files in a flowGraph object created by the flowGraph package.",
        "",
        "For a given flow cytometry sample and its cell populations, flowGraph generates SpecEnr values for each cell population or immunophenotype of based their expected proportion.",
        "By doing so, flowGraph accounts for the relation between cell populations and flags truly enriched cell populations rather than those with induced changes.",
        "",
        "The Summary_statistics folder contains p-values obtained from features associated with each immunophenotype or/ relation between immunophenotypes. The plots folder contains plots to help interpret those p-values.",
        "",
        "The Features folder contains features for:",
        paste0("- ", nrow(fg_get_meta(fg)), " flow cytometry samples."),
        paste0("- ", length(fg_get_markers(fg)), " markers: ",
               paste0(fg_get_markers(fg), collapse=", ")),
        paste0("- ", nrow(fg_graph$v),
               " cell population nodes (nodes for short) and ",
               nrow(fg_graph$e), " edges on ",
               max(fg_get_graph(fg)$v$phenolayer), " layers."),
        "",
        "Cell hierarchy plots: By default, for SpecEnr features, we generate two cell_hierarchy plots. The original one is where the colours represent difference between mean SpecEnr values across sample classes. We also add one where the colours represent difference betwee mean original values across sample classes. Original here is usually proportion: the feature used to create SpecEnr. Note that if a node is coloured lightly on the second plot, then the difference is very small, meaning the SpecEnr value may become sporadic. Therefore, when analyzing the plots, for most cases we recommend looking at most important cell populations as the ones with large difference in both plots.",
        ""), fc)

    close(fc)


    # save sample meta data and markers
    meta <- fg_get_meta(fg)
    utils::write.csv(meta, file=paste0(folder_path, "/sample_meta.csv"),
                     row.names=FALSE)


    # save features
    fg_feat <- fg_get_feature_all(fg)
    fg_feat_desc <- fg_get_feature_desc(fg)
    fn_dir <- paste0(folder_path, "/features/nodes")
    dir.create(fn_dir, recursive=TRUE, showWarnings=FALSE)
    a <- purrr::map(names(fg_feat$node), function(fn) {
        m <- as.matrix(fg_feat$node[[fn]])
        utils::write.csv(m, file=paste0(fn_dir,"/", fn, ".csv"))
    })
    sm <- fg_feat_desc$node
    utils::write.csv(as.matrix(sm), file=paste0(fn_dir,".csv"), row.names=FALSE)

    if (!is.null(fg_feat$edge)) {
        fe_dir <- paste0(folder_path, "/features/edges")
        dir.create(fe_dir, recursive=TRUE, showWarnings=FALSE)
        a <- purrr::map(names(fg_feat$edge), function(fe) {
            m <- as.matrix(fg_feat$edge[[fe]])
            # if (fe=="prop") fn <- "Proportion"
            utils::write.csv(m, file=paste0(fe_dir,"/", fe, ".csv"))
        })
        sm <- fg_feat_desc$edge
        utils::write.csv(sm, file=paste0(fe_dir,".csv"), row.names=FALSE)
    }


    # save summaries
    fg_summary_desc <- fg_get_summary_desc(fg)
    fg_summary <- fg_get_summary_all(fg)
    if (!is.null(fg_summary$node)) {
        sn_dir <- paste0(folder_path, "/summary_statistics/nodes")
        dir.create(sn_dir, recursive=TRUE, showWarnings=FALSE)
        tb <- fg_get_summary_tables(fg)
        utils::write.csv(tb, file=paste0(sn_dir, "/pvalues.csv"))
        sm <- fg_summary_desc$node
        data.table::fwrite(sm, file=paste0(sn_dir,".csv"), sep=",",
                           row.names=FALSE, col.names=TRUE)
    }

    if (!is.null(fg_summary$edge)) {
        se_dir <- paste0(folder_path, "/summary_statistics/edges")
        dir.create(se_dir, recursive=TRUE, showWarnings=FALSE)
        tb <- fg_get_summary_tables(fg, type="edge")
        utils::write.csv(tb, file=paste0(se_dir, "/pvalues.csv"))
        sm <- as.matrix(fg_summary_desc$edge)
        data.table::fwrite(sm, file=paste0(sn_dir,".csv"), sep=",",
                           row.names=FALSE, col.names=TRUE)
    }

    # save edge lists, graph layout and others for reproducibility
    # add folder and save id to object
    elist <- fg_get_el(fg)
    markers <- fg_get_markers(fg)
    pl <- fg_get_plot_layout(fg)
    fg@etc$save$id <<- list(id=stringi::stri_rand_strings(1,5))
    fg@etc$save$path <<- folder_path
    etc <- fg_get_etc(fg)
    save(fg_graph, elist, fg_summary, markers, pl, etc,
         file=paste0(folder_path, "/other.Rdata"))

    if (save_plots)
        fg_save_plots(fg, plot_path=paste0(folder_path,"/plots"), paired=paired)

    return(TRUE)
}

#' @title Load a flowGraph object from a specified folder path.
#' @description Load a flowGraph object from a specified folder path.
#' @param folder_path A string indicating the folder path to where a flowGraph
#'  object was saved using the \code{fg_save} function.
#' @return flowGraph object
#' @details see function \code{fg_save}
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos2)
#'  fg <- flowGraph(fg_data_pos2$count, class=fg_data_pos2$meta$class,
#'                  no_cores=no_cores)
#'
#'  fg_save(fg, "tmp")
#'  fg <- fg_load("tmp")
#' @seealso
#'  \code{\link[flowGraph]{fg_save}}
#' @rdname fg_load
#' @export
#' @importFrom data.table fread
#' @importFrom purrr map
#' @importFrom methods new
fg_load <- function(folder_path) {
    options(encoding="UTF-8")
    options(stringsAsFactors=FALSE)

    # load sample meta data and markers
    meta <- data.table::fread(paste0(folder_path, "/sample_meta.csv"), data.table=FALSE)

    # load features
    fn_dir <- paste0(folder_path, "/features/nodes")
    fn_files <- gsub(".csv","",list.files(fn_dir))
    feats <- list(node=purrr::map(fn_files, function(fn) {
        a <- data.table::fread(file=paste0(fn_dir,"/", fn, ".csv"))
        rownames_a <- a[,1]
        a <- as.matrix(a[,-1,drop=FALSE])
        rownames(a) <- unlist(rownames_a)
        colnames(a)[1] <- ""
        a
    }))
    names(feats$node) <- fn_files
    feats_desc <- list(node=data.table::fread(paste0(fn_dir,".csv"), data.table=FALSE))

    fe_dir <- paste0(folder_path, "/features/edges")
    if (dir.exists(fe_dir)) {
        fe_files <- gsub(".csv","",list.files(fe_dir))
        feats$edge <- purrr::map(fe_files, function(fe) {
            a <- data.table::fread(file=paste0(fe_dir,"/", fe, ".csv"))
            rownames_a <- a[,1]
            a <- as.matrix(a[,-1,drop=FALSE])
            rownames(a) <- unlist(rownames_a)
            a
        })
        feats_desc$edge <- data.table::fread(paste0(fe_dir,".csv"), data.table=FALSE)
        names(feats$edge) <- fe_files
    }


    # load summary descriptions
    summary_desc <- list()
    sn_dir <- paste0(folder_path, "/summary_statistics/nodes")
    if (dir.exists(sn_dir))
        summary_desc$node <- data.table::fread(paste0(sn_dir,".csv"), data.table=FALSE)

    se_dir <- paste0(folder_path, "/summary_statistics/edges")
    if (dir.exists(se_dir)) {
        summary_desc$edge <- data.table::fread(paste0(se_dir,".csv"), data.table=FALSE)
    }

    # load edge lists, graph layout and others for reproducibility
    # originally: load(paste0(folder_path, "/other.Rdata"))
    # fg_graph, elist, fg_summary, markers, pl, etc
    other_path <- paste0(folder_path, "/other.Rdata")
    attach(other_path, warn.conflicts=FALSE)

    other_path_ <- paste0("file:", other_path)

    # organize into flowGraph object
    methods::new(
        "flowGraph",
        feat=feats,
        feat_desc=feats_desc,
        summary=get("fg_summary", pos=other_path_),
        summary_desc=summary_desc,
        markers=get("markers", pos=other_path_),
        graph=get("fg_graph", pos=other_path_),
        edge_list=get("elist", pos=other_path_),
        meta=meta,
        plot_layout=get("pl", pos=other_path_),
        etc=get("etc", pos=other_path_)
    )
}


#' @title Retrieves a feature matrix.
#' @description Retrieves a feature matrix from a given flowGraph object,
#'  the feature type, and feature name.
#' @param fg flowGraph object.
#' @param type A string indicating feature type 'node' or 'edge'.
#' @param feature A string indicating feature name;
#' @return A numeric matrix of the specified feature values.
#' @details Returns \code{NULL} if the requested feature does not exist.
#' @examples
#'
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=1)
#'
#'  feature_matrix <- fg_get_feature(fg, type='node', feature='count')
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_get_feature_desc}}
#'  \code{\link[flowGraph]{fg_add_feature}}
#'  \code{\link[flowGraph]{fg_rm_feature}}
#'  \code{\link[flowGraph]{fg_get_summary}}
#' @rdname fg_get_feature
#' @export
fg_get_feature <- function(fg, type="node", feature="count") {
    type <- match.arg(type, c("node", "edge"))
    as.matrix(fg_get_feature_all(fg)[[type]][[feature]])
}

fg_get_feature_all <- function(fg) fg@feat


#' @title Retrieves the index of the requested summary.
#' @description Retrieves the index of the requested summary
#'  from a given flowGraph object.
#' @param fg flowGraph object.
#' @param type A string indicating feature type the summary was created for
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
#'  \code{type} (feature type: \code{node} or \code{edge}),
#'  \code{feature} (feature name), \code{test_name} (summary statistic name),
#'  \code{class} (class), \code{lable1}, and \code{label2} (class labels compared).
#'  See \code{\link[flowGraph]{fg_get_summary_desc}} for details.
#' @return An integer analagous to \code{index}.
#'  If both \code{index} and \code{summary_meta} are \code{NULL}, returns 1.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  no_cores=no_cores)
#'
#'  # set features to NULL to apply summary statistic to all features.
#'  fg <- fg_summary(fg, no_cores=no_cores, class="class", label1="control",
#'                   overwrite=FALSE, test_name="t", diminish=FALSE,
#'                   node_features=NULL, edge_features=NULL)
#'  show(fg)
#'
#'  index <- flowGraph:::fg_get_summary_index(
#'   fg, type="node", summary_meta=list(
#'     feature="SpecEnr", test_name="t", class="class",
#'     label1="control", label2="exp"))
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_get_summary_desc}}
#'  \code{\link[flowGraph]{fg_add_summary}}
#'  \code{\link[flowGraph]{fg_rm_summary}}
#'  \code{\link[flowGraph]{fg_plot}}
#' @rdname fg_get_summary_index
#' @export
fg_get_summary_index <- function(fg, type="node", index=NULL, summary_meta=NULL) {
    type <- match.arg(type, c("node", "edge"))
    if (is.null(fg_get_summary_desc(fg))) stop("summary not found")
    if (is.null(fg_get_summary_desc(fg)[[type]])) stop()
    if (is.null(index) & is.null(summary_meta)) {
        warning("no summary statistic index/summary_meta provided, returning index=1")
        return(1)
    }
    if (is.null(index)) {
        index <- which(apply(fg_get_summary_desc(fg)[[type]],1, function(x)
            identical(as.character(unlist(x)),
                      as.character(unlist(summary_meta)))))
        if (length(index) == 0) stop("summary not found")
    }
    return(index)
}


#' @title Retrieves a summary statistic.
#' @description Retrieves a summary statistic from a given flowGraph object;
#'  while \code{fg} is required, the user can choose to input parameters
#'  \code{summary_meta}, \code{index}, or all of \code{type},
#'  \code{feat}, \code{test_name}, \code{class}, \code{label1}, and \code{label2}.
#'  See \code{\link[flowGraph]{fg_get_summary_desc}} for details.
#' @param fg flowGraph object.
#' @param type A string indicating feature type the summary was created for
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
#' @param adjust_custom A function or a string indicating the
#' test adjustment method to use.
#'  If a string is provided, it should be one of
#'  \code{c("holm", "hochberg", "hommel",
#'  "bonferroni", "BH", "BY", "fdr", "none")} (see \code{p.adjust.methods}).
#'  If a function is provided, it should take as input
#'  a numeric vector and output the
#'  same vector adjusted.
#' @param summary_fun A function that takes in a matrix and outputs a
#'  vector the same length as the number of columns this matrix has.
#'  Set to \code{NULL} to not calculate this summary
#'  (i.e. returned list will not contain \code{m1} and \code{m2}).
#'  See \code{\link[flowGraph]{fg_get_feature_means}}.
#' @param adjust0_lim A vector of two numeric values indicating a range around
#'  0, default set to -0.1 and 0.1.
#' @param filter_adjust0 A numeric variable indicating what percentage of
#'  SpecEnr values compared (minimum) should be not close to 0.
#'  Set to 1 to not conduct filtering. Original p-values stored in
#'  \code{values_original}.
#' @param filter_es A numeric variable between 0 and 1 indicating what the
#'  Cohen's D value of the nodes/edges in question must be greater or
#'  equal to, to be significant.
#' @param filter_btwn_tpthres A numeric variable between 0 and 1 indicating the
#'  unadjusted T-test p-value threshold used to test whether the actual
#'  and expected feature values used to calculate the specified SpecEnr
#'  feature are significantly different for each sample class. Note this only
#'  needs to be specified for SpecEnr features. Combined
#'  with \code{filter_btwn_es}, we conduct three tests to understand if
#'  there is an actual large difference between actual and expected features:
#'  (1,2) T-test of significance between the actual and expected raw feature value
#'  (e.g. proportion) for samples in each of the compared classes, (3) and the
#'  T-test of significance between the differences of actual and
#'  expected feature values of the two classes. If any two of the three tests
#'  come out as insignificant, we set the p-value for the associated node/edge
#'  to 1.
#' @param filter_btwn_es A numeric variable between 0 and 1 indicating what the
#'  Cohen's D value of the nodes/edges in question must be greater or
#'  equal to, to be significant -- see \code{filter_btwn_tpthres}.
#' @param default_p_thres A numeric variable indicating the p-value threshold
#'  user is using. Currently, all nodes/edges not passing the \code{filter}
#'  criterion will be defaulted to 1; if this parameter is set, then all
#'  of these nodes/edges will be set to a minimum of \code{default_p_thres}.
#' @return A list containing elements on feature summary retrieved by the user
#'  as in the \code{summary} slot of
#'  \code{\link[flowGraph]{flowGraph-class}}.
#'  If \code{summary_fun} is not \code{NULL}, this list also includes:
#'  \itemize{
#'   \item{\code{m1}: a numeric vector the same length as \code{values};
#'   this is a summary of the samples compared e.g. mean.}
#'   \item{\code{m2}: a numeric vector the same length as \code{values};
#'   this is a summary of the samples compared e.g. mean.}
#'   \item{\code{cohensd}: a numberic vector indicating cohen's d values
#'    considering effect size.}
#'   \item{\code{cohensd_size}: a factor vector interpreting cohen's d values.}
#'   \item{\code{adjust0}: a numeric vector indicating the percentage of
#'    samples that have a SpecEnr value in the range of \code{adjust0_lim}
#'    around 0; if there are two classes of samples being compared, we output
#'    the smaller percentage between the two classes.}
#'   \item{\code{btwn}: a data frame containing columns:
#'    \itemize{
#'     \item{\code{tpv1}: unadjusted p-value calculated
#'      between the actual and expected raw feature values of class 1.}
#'     \item{\code{tpv2}: unadjusted p-value calculated
#'      between the actual and expected raw feature values of class 2.}
#'     \item{\code{cd1}: Cohen's D between the actual and expected raw
#'      feature values of class 1.}
#'     \item{\code{cd2}: Cohen's D between the actual and expected raw
#'      feature values of class 2.}
#'     \item{\code{btp}: unadjusted p-value calculated between the
#'      difference between actual and expected raw feature of the two classes.}
#'     \item{\code{bcd}: Cohen's D calculated between the
#'      difference between actual and expected raw feature of the two classes.}
#'     \item{\code{btp_}: unadjusted p-value calculated between the
#'      log ratio between actual and expected raw feature of the two classes.}
#'     \item{\code{bcd_}: Cohen's D calculated between the
#'      log ratio between actual and expected raw feature of the two classes.}
#'    }
#'   }
#'  }
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  no_cores=no_cores)
#'
#'  # set features to NULL to apply summary statistic to all features.
#'  fg <- fg_summary(fg, no_cores=no_cores, class="class", label1="control",
#'                   overwrite=FALSE, test_name="t", diminish=FALSE,
#'                   node_features=NULL, edge_features=NULL)
#'  show(fg)
#'
#'  feat_summ <- fg_get_summary(fg, type="node", summary_meta=list(
#'      feature="SpecEnr", test_name="t", class="class",
#'      label1="control", label2="exp"))
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_get_feature_means}}
#'  \code{\link[flowGraph]{fg_get_summary_desc}}
#'  \code{\link[flowGraph]{fg_add_summary}}
#'  \code{\link[flowGraph]{fg_rm_summary}}
#'  \code{\link[flowGraph]{fg_get_feature}}
#' @rdname fg_get_summary
#' @export
#' @importFrom purrr map_lgl
#' @importFrom effsize cohen.d
fg_get_summary <- function(
    fg, type="node", index=NULL, summary_meta=NULL, adjust_custom="byLayer",
    summary_fun=colMeans,
    adjust0_lim=c(-.1,.1), filter_adjust0=1, filter_es=0,
    filter_btwn_tpthres=.05, filter_btwn_es=.5,
    default_p_thres=1
) {
    type <- match.arg(type, c("node", "edge"))
    index <- fg_get_summary_index(fg,type=type, index, summary_meta)
    summary_meta <- unlist(fg_get_summary_desc(fg)[[type]][index,])
    if (!grepl("SpecEnr",unlist(summary_meta["feat"])))
        filter_adjust0 <- 1

    sl <- fg_get_summary_all(fg)[[type]][[index]]
    sl <- append(sl, list(
        id1=fg_get_meta(fg)[,summary_meta["class"]]==summary_meta["label1"],
        id2=fg_get_meta(fg)[,summary_meta["class"]]==summary_meta["label2"]))

    # adjust function
    if (!is.function(adjust_custom))
        if (adjust_custom=="none") {
            adjust_custom <- function(x) x
        } else if (adjust_custom!="byLayer") {
            ac <- adjust_custom
            adjust_custom <- function(x) stats::p.adjust(x,method=ac)
        }

    pna <- is.na(sl$values)
    pnaw <- which(!pna)
    if (!is.function(adjust_custom)) {
        fg_graph <- fg_get_graph(fg)
        lyrs <- fg_graph$v$phenolayer[!pna]
        if (type=="edge")
            lyrs <- fg_graph$v$phenolayer[match(
                fg_graph$e$to[!pna], fg@graph$v$phenotype)]
        lyrsu <- unique(lyrs)
        for (lyr in lyrsu) {
            lyri <- lyrs==lyr
            sl$values[pnaw[lyri]] <- sl$values[pnaw[lyri]] * sum(lyri)
        }
        sl$values <- sl$values * length(lyrsu)
    } else {
        sl$values[pna] <- adjust_custom(sl$values[pna])
    }
    sl$values[pna] <- 1
    belowthres <- sl$values<default_p_thres

    # colmeans
    if (!is.null(summary_fun)) {
        feat <- unlist(fg_get_summary_desc(fg)[[type]]$feat[index])
        m1 <- fg_get_feature_means(
            fg, type, feat, id=sl$id1, summary_fun=summary_fun)
        m2 <- fg_get_feature_means(
            fg, type, feat, id=sl$id2, summary_fun=summary_fun)
        sl <- append(sl, list(m1=m1, m2=m2))
    }

    # cohen's d
    ename <- paste(summary_meta[c("feat", "class","label1","label2")], collapse="_")
    fg_etc <- fg_get_etc(fg)
    if (is.null(fg_etc$effect_size[[type]]))
        fg_etc$effect_size[[type]] <- list()
    if (is.null(fg_etc$effect_size[[type]][[ename]]))
        fg_etc$effect_size[[type]][[ename]] <-
            fg_cohensd_(fg@feat[[type]][[feat]], sl$id1, sl$id2)

    sl <- append(sl, fg_etc$effect_size[[type]][[ename]])

    if (grepl("SpecEnr", summary_meta["feat"])) {
        aname <- paste(summary_meta[c("feat","class","label1","label2")], collapse="_")
        if (is.null(fg_etc$adjust0[[type]]))
            fg_etc$adjust0[[type]] <- list()
        if (is.null(fg_etc$adjust0[[type]][[aname]]))
            fg_etc$adjust0[[type]][[aname]] <-
                fg_adjust0_(fg_get_feature_all(fg)[[type]][[feat]],
                            sl$id1,sl$id2, adjust0_lim)

        sl$adjust0 <- fg_etc$adjust0[[type]][[aname]]

        sm_name <- paste0(summary_meta[-2], collapse="_")
        if (is.null(fg_etc$actualVSexpect[[type]]))
            fg_etc$actualVSexpect[[type]] <- list()
        if (is.null(fg_etc$actualVSexpect[[type]][[sm_name]])) {
            feat_ <- se_feats(feat)
            mp <- fg_get_feature(fg, type=type, feature=feat_[2])
            mep <- fg_get_feature(fg, type=type, feature=feat_[3])
            fg_etc$actualVSexpect[[type]][[sm_name]] <- fg_btwn_(
                sl$id1, sl$id2, mp, mep,
                btwn_test_custom=function(x,y)
                    tryCatch(stats::t.test(x,y)$p.value, error=function(e) 1))
        }
        sl$btwn <- df_ps <- fg_etc$actualVSexpect[[type]][[sm_name]]

        sl0 <- sl$values
        if (filter_btwn_es>0 | filter_btwn_tpthres<1) {

            sig1 <- df_ps$tpv1<filter_btwn_tpthres
            sig2 <- df_ps$tpv2<filter_btwn_tpthres
            sig3 <- df_ps$btp<filter_btwn_tpthres
            sigcands <- (sig1 | sig2) & sig3

            sl$values_original <- sl0
            sl$values[!sigcands & belowthres] <- default_p_thres
        }

        # 0's
        if (filter_adjust0<1 & grepl("SpecEnr", unlist(summary_meta["feat"]))) {
            sl$values_original <- sl0
            sl$values[sl$adjust0>filter_adjust0 & belowthres] <- default_p_thres
        }
        # cohen's d
        if (filter_es>0) {
            sl$values_original <- sl0
            sl$values[abs(sl$cohensd)<filter_es & belowthres] <- default_p_thres
        }
    }
    fg@etc <<- fg_etc
    sl$values[is.na(sl$values) | sl$values>1] <- 1
    return(sl)
}

fg_get_summary_all <- function(fg) fg@summary


#' @title Retrieves a table containing all node or edge summary statistics.
#' @description Retrieves a table containing all node or edge summary statistics
#'  given a flowGraph object.
#' @param fg flowGraph object.
#' @param type A string indicating feature type the summaries the user
#'  wants to retrieve were created for, 'node' or 'edge'.
#' @return A list; this output is the same as that of function
#'  \code{fg_get_graph} with additional columns.
#'  These columns contain summary statistics from the \code{summary}
#'  slot of the flowGraph object. These columns are named:
#'  <feature type: node/edge>.<feature>.<summary name>.<class>.<class labels>.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  no_cores=no_cores)
#'
#'  fg <- fg_summary(fg, no_cores=no_cores, class="class", label1="control",
#'                   overwrite=FALSE, test_name="t", diminish=FALSE)
#'  show(fg)
#'
#'  feat_summ_table_node <- fg_get_summary_tables(fg, type="node")
#'  head(feat_summ_table_node)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_get_feature_means}}
#'  \code{\link[flowGraph]{fg_get_summary_desc}}
#'  \code{\link[flowGraph]{fg_add_summary}}
#'  \code{\link[flowGraph]{fg_rm_summary}}
#'  \code{\link[flowGraph]{fg_get_summary}}
#'  \code{\link[flowGraph]{fg_get_feature}}
#' @rdname fg_get_summary_tables
#' @export
fg_get_summary_tables <- function(fg, type="node") {
    type <- match.arg(type, c("node", "edge"))
    if (type=="node") {
        gr <- fg_get_graph(fg)$v
        gr <- gr[,!colnames(gr)%in%c("x", "y")]
    } else {
        gr <- fg_get_graph(fg)$e
        gr <- gr[,!colnames(gr)%in%c("from.x", "from.y", "to.x", "to.y")]
    }
    if (!type%in%names(fg_get_summary_all(fg))) stop("no summaries found")

    desc <- fg_get_summary_desc(fg)[[type]]
    descn <- apply(desc, 1, function(x) paste0(x, collapse="."))
    descn <- paste0(type, ".", descn)
    for (ft1 in seq_len(length(descn)))
        gr[[descn[ft1]]] <- fg_get_summary(fg, type, index=ft1)$values
    return(gr)
}


#' @title Retrieves and/or recalculates a feature description table.
#' @description Retrieves and/or recalculates a feature
#'  description table for a given flowGraph object.
#' @param fg flowGraph object.
#' @param re_calc A logical variable specifying whether or not a feature summary
#'  should be re-calculated or directly retrieved from \code{fg}.
#' @return A data frame where each row contains information on a feature
#' from the given flowGraph object; its columns is as in the \code{feat_desc}
#' slot of \code{\link[flowGraph]{flowGraph-class}}.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  no_cores=no_cores)
#'
#'  fg_get_feature_desc(fg, re_calc=TRUE)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_get_feature}}
#'  \code{\link[flowGraph]{fg_add_feature}}
#'  \code{\link[flowGraph]{fg_rm_feature}}
#'  \code{\link[flowGraph]{fg_get_summary_desc}}
#' @rdname fg_get_feature_desc
#' @export
fg_get_feature_desc <- function(fg, re_calc=FALSE) {
    if (re_calc) {
        fg_feat <- fg_get_feature_all(fg)
        result1 <- as.data.frame(sapply(names(fg_feat$node), function(x)
            summary_table(fg_feat$node[[x]], x)))

        result2 <- NULL
        if (length(fg_feat$edge) >
            0) {
            result2 <- as.data.frame(sapply(
                names(fg_feat$edge), function(x)
                    summary_table(fg_feat$edge[[x]], x)))
        }
        return(list(node=result1, edge=result2))
    }
    fg@feat_desc
}


#' @title Retrieves a feature summary description table.
#' @description Retrieves a feature summary description table for
#'  a given flowGraph object.
#' @param fg flowGraph object.
#' @return A data frame where each row contains information
#'  on a feature summary from \code{fg}:
#' \itemize{
#'   \item{\code{type}: feature type (i.e. 'node' or 'edge').}
#'   \item{\code{feat}: feature name.}
#'   \item{\code{test_name}: summary name.}
#'   \item{\code{class}: class or the column name of \code{fg_get_meta(fg)}
#'    whose values represent the class label of each sample on which
#'    the summary was created for.}
#'   \item{\code{label1}: A string from the \code{class} column of the
#'    \code{meta} slot indicating the label of samples compared.}
#'   \item{\code{label2}: A string from the \code{class} column of the
#'    \code{meta} slot indicating the label of samples compared.}
#' }
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  no_cores=no_cores)
#'
#'  fg_get_summary_desc(fg)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_get_summary}}
#'  \code{\link[flowGraph]{fg_add_summary}}
#'  \code{\link[flowGraph]{fg_rm_summary}}
#'  \code{\link[flowGraph]{fg_get_feature_desc}}
#' @rdname fg_get_summary_desc
#' @export
fg_get_summary_desc <- function(fg) fg@summary_desc



## get/replace fg_get_meta meta data

#' @title Retrieves sample meta.
#' @description Retrieves sample meta from a given flowGraph object.
#' @param fg flowGraph object.
#' @return A data frame containing sample meta data.
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_replace_meta}}
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=no_cores)
#'  head(fg_get_meta(fg))
#'
#' @export
#' @rdname fg_get_meta
fg_get_meta <- function(fg) fg@meta



## get graph

#' @title Retrieves a graph list from a given flowGraph object.
#' @description Retrieves a graph list from a given flowGraph object.
#' @param fg flowGraph object.
#' @return A list containing two data frames (\code{v} and ]code{e})
#'  from the \code{graph} slot of the given \code{flowGraph} object containing
#'  information on the cell populations phenotype nodes and edges representing
#'  relation between cell populations.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=no_cores)
#'  gr <- fg_get_graph(fg)
#'  head(gr$v)
#'  head(gr$e)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_plot}}
#'  \code{\link[flowGraph]{ggdf}}
#'  \code{\link[flowGraph]{plot_gr}}
#' @rdname fg_get_graph
#' @export
fg_get_graph <- function(fg) fg@graph

fg_get_el <- function(fg) fg@edge_list


## get markers
#' @title Retrieves the markers from a given flowGraph object.
#' @description Retrieves the markers from a given flowGraph object.
#' @param fg flowGraph object.
#' @return A character vector containing the markers used in a flowGraph object.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=no_cores)
#'  fg_get_markers(fg)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#' @rdname fg_get_markers
#' @export
fg_get_markers <- function(fg) fg@markers

#' @title Retrieves feature summaries.
#' @description Retrieves a feature summary (e.g. \code{colMeans}) for samples
#'  specified by sample id's \code{id} OR class label
#'  \code{label} for class \code{class} given a feature specified
#'  by \code{type} and \code{feat}.
#' @param fg flowGraph object.
#' @param type A string indicating feature type the summary was created for
#'  'node' or 'edge'.
#' @param feature A string indicating feature name the summary was created for;
#' @param class A string corresponding to a column name of the \code{meta} slot of
#'  \code{fg} whose values represent the class label of each sample
#'  on which the summary was created to compare or analyze;
#' @param label A string indicating a class label.
#' @param id A string vector containing the sample id's corresponding to the
#'  \code{id} column of the \code{meta} slot of \code{fg}.
#' @param summary_fun A function that takes in a matrix and outputs a
#'  vector the same length as the number of columns this matrix has.
#' @return A list containing two numeric vectors calculated using the
#'  \code{summary_fun} function on the subset of samples specified by
#'  sample id's \code{id} OR class label
#'  \code{label} for class \code{class} from a feature matrix
#'  specified by \code{type} and \code{feat}.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  no_cores=no_cores)
#'  fg <- fg_summary(fg, no_cores=no_cores, class="class", label1="control",
#'                   overwrite=FALSE, test_name="t", diminish=FALSE)
#'  show(fg)
#'  feat_mean <- fg_get_feature_means(fg, type="node", feature="count",
#'                                    class="class", label="control")
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_get_summary_desc}}
#'  \code{\link[flowGraph]{fg_add_summary}}
#'  \code{\link[flowGraph]{fg_rm_summary}}
#'  \code{\link[flowGraph]{fg_get_summary}}
#' @rdname fg_get_feature_means
#' @export
fg_get_feature_means <- function(
    fg, type=c("node","edge"), feature="count",
    class=NULL, label=NULL, id=NULL,
    summary_fun=colMeans
) {
    m <- fg_get_feature(fg, type, feature)
    if (!is.null(id)) {
        m1 <- m[id,,drop=FALSE]
    } else {
        m1 <- m[fg_get_meta(fg)[[class]]==label,,drop=FALSE]
    }
    summary_fun(as.matrix(m1))
}

fg_get_etc <- function(fg) fg@etc

fg_get_plot_layout <- function(fg) fg@plot_layout


# given number of markers, output number of nodes, edges
get_graph_info <- function(m) {
    npl <- epl <- rep(0,m)
    for (i in seq_len(m)) {
        npl[i] <- choose(m,i)*(2^i)
        epl[i] <- npl[i]*i
    }
    n <- 3^m
    e <- n*2*m/3
    return (list(n=n, npl=npl, e=e, epl=epl))
}


# given a "SpecEnr" feature, get the raw and expected features
se_feats <- function(feature) {
    if (feature=="SpecEnr")
        return(c(feature, "prop", "expect_prop"))
    if (grepl("SpecEnr", feature)) {
        fl <- gsub("SpecEnr_","",feature)
        return(c(feature, fl, paste0("expect_",fl)))
    }
    return(feature)
}
