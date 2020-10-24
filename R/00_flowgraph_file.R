
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
#'  data(fg_data_pos15)
#'  fg <- flowGraph(fg_data_pos15$count, class=fg_data_pos15$meta$class,
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
    if (base::is.null(folder_path)) {
        if (base::is.null(fg@etc$save$path))
            stop("please provide valid folder path")
        folder_path <- fg@etc$save$path
    }
    if (!dir.exists(folder_path))
        dir.create(folder_path, recursive=TRUE, showWarnings=FALSE)
    if (base::length(dir(folder_path, all.files=TRUE))==0) {
        stop_save <- TRUE

        etc_dir <- paste0(folder_path,"/etc/etc.rds")
        if (file.exists(etc_dir))
            if (readRDS(etc_dir)$save$id == fg@etc$save$id)
                stop_save <- FALSE

        if (stop_save)
            stop("the folder provided is non-empty")
    }

    # create readme file
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
        paste0("- ", nrow(fg@meta), " flow cytometry samples."),
        paste0("- ", length(fg@markers), " markers: ",
               paste0(fg@markers, collapse=", ")),
        paste0("- ", base::nrow(fg@graph$v),
               " cell population nodes (nodes for short) and ",
               base::nrow(fg@graph$e), " edges on ",
               max(fg@graph$v$phenolayer), " layers."),
        "",
        "Cell hierarchy plots: By default, for SpecEnr features, we generate two cell_hierarchy plots. The original one is where the colours represent difference between mean SpecEnr values across sample classes. We also add one where the colours represent difference betwee mean original values across sample classes. Original here is usually proportion: the feature used to create SpecEnr. Note that if a node is coloured lightly on the second plot, then the difference is very small, meaning the SpecEnr value may become sporadic. Therefore, when analyzing the plots, for most cases we recommend looking at most important cell populations as the ones with large difference in both plots.",
        ""), fc)

    close(fc)


    # save sample meta data and markers
    meta <- fg@meta
    utils::write.csv(meta, file=paste0(folder_path, "/sample_meta.csv"),
              row.names=FALSE)


    # save features
    fn_dir <- paste0(folder_path, "/features/nodes")
    dir.create(fn_dir, recursive=TRUE, showWarnings=FALSE)
    a <- purrr::map(names(fg@feat$node), function(fn) {
        m <- as.matrix(fg@feat$node[[fn]])
        # if (fn=="count") fn <- "Cell_count"
        # if (fn=="prop") fn <- "Proportion"
        # if (fn=="expect_prop") fn <- "Expected_proportion"
        utils::write.csv(m, file=paste0(fn_dir,"/", fn, ".csv"))
    })
    sm <- fg@feat_desc$node
    utils::write.csv(as.matrix(sm), file=paste0(fn_dir,".csv"), row.names=FALSE)

    if (!base::is.null(fg@feat$edge)) {
        fe_dir <- paste0(folder_path, "/features/edges")
        dir.create(fe_dir, recursive=TRUE, showWarnings=FALSE)
        a <- purrr::map(names(fg@feat$edge), function(fe) {
            m <- as.matrix(fg@feat$edge[[fe]])
            # if (fe=="prop") fn <- "Proportion"
            utils::write.csv(m, file=paste0(fe_dir,"/", fe, ".csv"))
        })
        sm <- fg@feat_desc$edge
        utils::write.csv(sm, file=paste0(fe_dir,".csv"), row.names=FALSE)
    }


    # save summaries
    if (!base::is.null(fg@summary$node)) {
        sn_dir <- paste0(folder_path, "/summary_statistics/nodes")
        dir.create(sn_dir, recursive=TRUE, showWarnings=FALSE)
        tb <- fg_get_summary_tables(fg)
        utils::write.csv(tb, file=paste0(sn_dir, "/pvalues.csv"))
        sm <- fg@summary_desc$node
        data.table::fwrite(sm, file=paste0(sn_dir,".csv"), sep=",",
                           row.names=FALSE, col.names=TRUE)
    }

    if (!base::is.null(fg@summary$edge)) {
        se_dir <- paste0(folder_path, "/summary_statistics/edges")
        dir.create(se_dir, recursive=TRUE, showWarnings=FALSE)
        tb <- fg_get_summary_tables(fg, type="edge")
        utils::write.csv(tb, file=paste0(se_dir, "/pvalues.csv"))
        sm <- as.matrix(fg@summary_desc$edge)
        data.table::fwrite(sm, file=paste0(sn_dir,".csv"), sep=",",
                           row.names=FALSE, col.names=TRUE)
    }

    # save edge lists, graph layout and others for reproducibility
    fo_dir <- paste0(folder_path,"/etc")
    dir.create(fo_dir, showWarnings=FALSE)

    gr <- fg@graph
    saveRDS(gr, file=paste0(fo_dir, "/graph.rds"))
    el <- fg@edge_list
    saveRDS(el, file=paste0(fo_dir, "/edge_list.rds"))

    ls <- fg@summary
    saveRDS(ls, file=paste0(fo_dir, "/summary.rds"))
    markers <- fg@markers
    saveRDS(markers, file=paste0(fo_dir, "/markers.rds"))
    pl <- fg@plot_layout
    saveRDS(pl, file=paste0(fo_dir, "/plot_layout.rds"))
    etc <- fg@etc
    saveRDS(etc, file=paste0(fo_dir, "/etc.rds"))

    if (save_plots)
        fg_save_plots(fg, plot_path=paste0(folder_path,"/plots"), paired=paired, ...)

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
#'  data(fg_data_pos15)
#'  fg <- flowGraph(fg_data_pos15$count, class=fg_data_pos15$meta$class,
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
    fo_dir <- paste0(folder_path,"/etc")
    gr <- readRDS(paste0(fo_dir, "/graph.rds"))
    el <- readRDS(paste0(fo_dir, "/edge_list.rds"))

    ls <- readRDS(paste0(fo_dir, "/summary.rds"))
    markers <- readRDS(paste0(fo_dir, "/markers.rds"))
    pl <- readRDS(paste0(fo_dir, "/plot_layout.rds"))
    etc <- readRDS(paste0(fo_dir, "/etc.rds"))

    # organize into flowGraph object
    fg <- methods::new(
        "flowGraph",
        feat=feats,
        feat_desc=feats_desc,
        summary=ls, summary_desc=summary_desc,
        markers=markers,
        graph=gr, edge_list=el,
        meta=meta, plot_layout=pl,
        etc=etc
    )
}


