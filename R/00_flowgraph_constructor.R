#' @name flowGraph
#' @title flowGraph object constructor.
#' @description Initializes a \code{flowGraph} object given the cell counts
#'  for one or more flow cytometry sample(s).
#'  The flowGraph object returned holds meta data
#'  for each sample, each cell population node,
#'  edges representing how each cell population node relate to one another,
#'  and features for these nodes and edges.
#' @param input_ Any of the following:
#' \itemize{
# #'   \item{a \code{Phenotypes} object from the \code{ftf***} package.}
#'   \item{a numeric matrix or vector of the cell counts;
#'   its column/names must be the phenotype names and its rownames
#'   must be sample ID's.}
# #'   \item{a list of \code{Phenotypes} object(s).}
# #'   \item{a vector of \code{Phenotypes} object full paths.}
#' }
#' All input samples should have the same \code{markers} and
#' \code{partitionsPerMarker}.
# #' If given \code{Phenotype} object(s),
# #' the function will look to generate phenotype names from the
# #' \code{Phenotypes@@PhenoCodes} or \code{rownames(Phenotypes@@MFIs)} slots.
#' @param meta A data frame with meta data for each \code{Phenotypes} or sample;
#'  One of its column names should be "id" whose values correspond to
#'  the name of each \code{Phenotypes} object. We also recommend for it to have
#'  a column named "class" where one of its unique values is "control".
#' @param class A string corresponding to the column name or index
#'  of \code{meta} whose values represent
#'  the class label of each sample; OR a vector the same length
#'  as the the number of samples in \code{input_} specifiying the class
#'  of each given sample --- this vector will be appended to \code{meta} under
#'  column name \code{class}.
#' @param no_cores An integer indicating how many cores to parallelize on.
#' @param markers A string vector of marker names used in \code{input_}.
#' @param layout_fun A string of a function from the \code{igraph} package that
#'  indicates what layout should be used if a cell hierarchy is to be ploted;
#'  all such functions have prefix \code{layout_}. This is defaulted to
#'  e.g. \code{layout_fun="layout.reingold.tilford"}.
#' @param max_layer And integer indicating the maximum layer in the cell
#'  hierarchy to analyze; set to `NULL` to analyze all layers.
#' @param cumsumpos A logical variable indicating whether
#'  or not to cumulate cell counts;
#'  this applies only when \code{partitionsPerMarker > 3} and will convert
#'  e.g. the count of A+ or A++ into the sum of the counts of
#'  A+, A++, A+++, ..., or A++, A+++, ... .
#' @param prop A logical variable indicating whether or not to
#'  calculate the proportion feature;
#'  this can be done later on with \code{flowGraph_prop}.
#' @param specenr logical variable: whether or not to calculate
#'  the SpecEnr feature, Default: T
#' @param path A string indicating the folder path to where the flowGraph
#'  object should save its elements, Default = NULL (don't save).
#' @param calculate_summary A logical variable indicating whether or not to
#'  calculate the summary statistics for SpecEnr based on default parameters
#'  using the \code{fg_summary} summary function on class
#'  specified in parameter \code{class}.
#' @param node_features A string vector indicating which node feature(s)
#'  to perform summary statistics on; set to \code{NULL} or \code{"NONE"}
#'  and the function will perform summary statistics on all or no node features.
#' @param edge_features A string vector indicating which edge feature(s)
#'  to perform summary statistics on; set to \code{NULL} or \code{"NONE"}
#'  and the function will perform summary statistics on all or no edge features.
#' @param test_name A string with the name of the test you are performing.
#' @param test_custom See \code{\link[flowGraph]{fg_summary}}.
#' @param diminish A logical variable; applicable if \code{calculate_summary} is
#'  \code{TRUE}; see \code{\link[flowGraph]{fg_summary}}.
#' @param label1 A string indicating a class label in
#'  \code{fg_get_meta(fg)[,class]}; set to \code{NULL} if you would
#'  like to compare all classes aganst all classes;
#'  applicable if \code{calculate_summary} is \code{TRUE}.
#' @param label2 A string indicating a class label in
#'  \code{fg_get_meta(fg)[,class]};
#'  applicable if \code{calculate_summary} is \code{TRUE}.
# #' @param normalize A logical variable indicating whether or not to
# #'  calculate normalized cell count;
# #'  this can be done later on with \code{fg_feat_node_norm}.
# #' @param norm_ind See \code{\link[flowGraph]{fg_feat_node_norm}}.
# #' @param norm_layer See \code{\link[flowGraph]{fg_feat_node_norm}}.
# #' @param norm_path See \code{\link[flowGraph]{fg_feat_node_norm}}.
#' @param save_plots A logical indicating whether or not to save plots.
# #' @param ... Other parameters used in the \code{fg_save_plots} function
# #'  of the \code{edgeR} package. See \code{\link[flowGraph]{tmm}}.
#' @return flowGraph object
#' @details \code{flowGraph} is the constructor for the \code{flowGraph} object.
#'  The user can choose to input as \code{input_} a vector, a \code{Phenotypes}
#'  object (meaning there is only one sample), a matrix, or a \code{Phenotypes}
#'  object list. If the user is also inputting a sample meta data frame, it
#'  must contain a \code{id} column corresponding to sample names.
#' @examples
#'
#' samplen <- 10
#' meta_file <- data.frame(
#'     id=1:samplen,
#'     class=append(rep("control", samplen/2), rep("exp", samplen/2)),
#'     stringsAsFactors=FALSE
#' )
#'
#'
#' ## using the constructor -----------------------
#'
#' data(fg_data_pos30)
#'
#' # input: vector of load-able Phenotypes paths
#' fg <- flowGraph(fg_data_pos30$count[1,], no_cores=no_cores)
#'
#' # input: matrix + vector of class corresponding to samples
#' fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                 no_cores=no_cores)
#' # - save to file directly
#' # fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#' #                no_cores=no_cores, path="path_to_folder)
#'
#' # input: matrix + meta data frame
#' # fg <- flowGraph(fg_data_pos30$count, meta=fg_data_pos30$meta,
#' #                 no_cores=no_cores)
#'
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_get_feature}}
#'  \code{\link[flowGraph]{fg_get_feature_desc}}
#'  \code{\link[flowGraph]{fg_get_summary}}
#'  \code{\link[flowGraph]{fg_get_summary_desc}}
#'  \code{\link[flowGraph]{fg_add_feature}}
#'  \code{\link[flowGraph]{fg_rm_feature}}
#'  \code{\link[flowGraph]{fg_add_summary}}
#'  \code{\link[flowGraph]{fg_rm_summary}}
#'  \code{\link[flowGraph]{fg_gsub_markers}}
#'  \code{\link[flowGraph]{fg_gsub_ids}}
#'  \code{\link[flowGraph]{fg_merge_samples}}
#'  \code{\link[flowGraph]{fg_extract_samples}}
#'  \code{\link[flowGraph]{fg_extract_phenotypes}}
#'  \code{\link[flowGraph]{fg_merge}}
#'  \code{\link[doParallel]{registerDoParallel}}
#'  \code{\link[Matrix]{Matrix}}
#' @rdname flowGraph
#' @export
#' @importFrom stringr str_split
#' @importFrom stringi stri_rand_strings
#' @importFrom future plan multiprocess
#' @importFrom furrr future_map
#' @importFrom purrr map_lgl map_chr map_dfr map compact map_int
#' @importFrom Matrix Matrix
#' @importFrom methods new
#' @importFrom igraph layout.reingold.tilford
flowGraph <- function(
    input_, meta=NULL, class="class", no_cores=1, markers=NULL,
    layout_fun="layout.reingold.tilford", # layout.circle

    max_layer=NULL,
    cumsumpos=FALSE, # whether to make positives +=+/++ (cumsum)
    prop=TRUE, specenr=TRUE,

    path=NULL,

    # summary parameters
    calculate_summary=TRUE, node_features="SpecEnr", edge_features="NONE",
    test_name="t_diminish",
    test_custom="t",
    diminish=TRUE, label1=NULL, label2=NULL,

    # plotting parameters
    save_plots=FALSE

) {
    options(stringsAsFactors=FALSE)

    if (!is.null(meta))
        if (!"id"%in%colnames(meta))
            stop("meta must have a column named \"id\" containing sample ID's")

    start <- start1 <- Sys.time()
    message("preparing input; ")

    # convert input_ into a cell count matrix ----------------------------------
    msg <- "make sure input type is correct"

    phenocode <- NULL
    if (any(class(input_)%in%c("integer","numeric")) |
        grepl("matrix",class(input_), ignore.case=TRUE)) {
        ## feature: count (sample x cell population)
        if (is.null(dim(input_))) {
            mc <- mc0 <- matrix(input_,nrow=1)
            colnames(mc) <- names(input_)
            rownames(mc) <- "s1"
        } else {
            mc <- mc0 <- input_
            if (is.null(rownames(input_)))
                rownames(input_) <-
                    paste0("s",seq_len(nrow(input_)))
            rownames(mc) <- rownames(input_)
        }

        ## make meta for samples
        if (is.null(meta)) sample_id <- rownames(mc)

        ## make meta for cell populations
        if (!any(grepl("[+-]",colnames(mc)))) stop(msg)

        phen <- colnames(mc)
        if (is.null(markers)) {
            markers <- unlist(stringr::str_split(phen,"[-+]+"))
            markers <- gsub("_","",markers) # underscore between marker conditions
            markers <- unique(markers[markers!=""])
        }

    } else {
        stop("input must be a numeric sample x cell population matrix.")
    }

    # prepare sample meta ----------------------------------
    if (!is.null(meta)) {
        if (nrow(meta)!= nrow(mc))
            stop("meta data and cell count matrix have different row counts")
        if (!"id"%in%colnames(meta))
            stop("meta must have an \'id\' column, id,
                 whose values corresponding with sample names")
        meta <- as.data.frame(meta)
        factor_col <- purrr::map_lgl(meta, is.factor)
        if (any(factor_col)) warning("converting factors in meta into strings")
        for (j in which(factor_col)) meta[[j]] <- as.character(meta[[j]])
    } else {
        meta <- data.frame(id=sample_id)
    }
    # insert class if there is one
    if (!identical(class, "class")) {
        if (length(class)==length(sample_id)) {
            if ("class"%in%colnames(meta)) {
                meta$class_ <- class
            } else {
                meta$class <- class
            }
        }
    }

    ## feature: count (sample x cell population) -----------
    message("preparing edge lists: mapping out cell population relations")

    keepinds <- apply(mc, 2, function(x) any(x>0))
    keepinds <- apply(mc, 2, function(x) any(x>0))
    mc <- mc[,keepinds,drop=FALSE]
    colnames(mc) <- phen[keepinds]
    rownames(mc) <- meta$id
    mc <- Matrix::Matrix(mc, sparse=TRUE) # rm all 0 cols
    if (!is.null(phenocode)) phenocode <- phenocode[keepinds]

    ## make meta for cell populations FINAL -----------------
    meta_cell <- get_phen_meta(phen[keepinds],phenocode)
    if (any(grepl("_",meta_cell$phenotype))) {
        meta_cell$phenotype_ <- meta_cell$phenotype
        meta_cell$phenotype <- gsub("_","", meta_cell$phenotype)
        colnames(mc) <- meta_cell$phenotype
    }
    # get rid of lower layers if requested
    if (!is.null(max_layer)) {
        phen_id <- meta_cell$phenolayer <= max_layer
        meta_cell <- meta_cell[phen_id,,drop=FALSE]
        mc <- mc[,phen_id,drop=FALSE]
    }


    # make parent list (for each cell popultion, list its parents)
    pccell <- get_phen_list(meta_cell=meta_cell, no_cores=no_cores)
    edf <- pccell$edf
    pchild <- pccell$pchild # not used, returned
    pparen <- pccell$pparen

    # extract cell populations that have parents to calculate
    # expected proportions for; just to be safe
    cells1 <- meta_cell$phenotype[meta_cell$phenolayer==1]
    cells1 <- append("",cells1[cells1%in%names(pparen)])
    cells <- meta_cell$phenotype[meta_cell$phenolayer>1]
    cells <- cells[cells%in%names(pparen)]
    cells_ <- append(cells1,cells)

    meta_cell <- meta_cell[match(cells_,meta_cell$phenotype),]
    rooti <- which(meta_cell$phenolayer==0)

    mc <- mc[,match(cells_, colnames(mc)),drop=FALSE] # final cpops

    # trim/order list of parents and children
    pparen_ <- pparen[names(pparen)%in%cells_]
    pparen_[cells] <- purrr::map(pparen_[cells], function(x) {
        a <- x[x%in%cells_]
        if (length(a)==0) return(NULL)
        a
    })
    pparen_ <- purrr::compact(pparen_)

    pchild_ <- pchild[names(pchild)%in%cells_]
    pchild_ <- purrr::map(pchild_, function(x) {
        a <- x
        for (i in seq_len(length(a)))
            a[[i]] <- x[[i]][x[[i]]%in%cells_]
        if (all(purrr::map_int(a, length)==0)) return(NULL)
        a
    })
    pchild_ <- purrr::compact(pchild_)

    edf <- edf[edf$to%in%cells_ & edf$from%in%cells_,,drop=FALSE]

    time_output(start1)
    start1 <- Sys.time()
    message("calculating features")


    ## list of outputs
    gr <- list(e=edf, v=meta_cell)
    gr <- set_layout_graph(gr, layout_fun) # layout cell hierarchy

    if (is.null(meta)) meta <- data.frame(id=sample_id)

    desc <- summary_table(mc,"count")
    fg <- methods::new(
        "flowGraph",
        feat=list(node=list(count=mc), edge=list()),
        feat_desc=list(node=desc),
        markers=markers,
        graph=gr, edge_list=list(child=pchild_,parent=pparen_),
        meta=meta, plot_layout=as.character(substitute(layout_fun)),
        etc=list(cumsumpos=FALSE, class_mean_normalized=FALSE)
    )
    time_output(start1)
    start1 <- Sys.time()

    ## features ----------------------------
    cumsumpos <- cumsumpos & any(grepl("3",meta_cell$phenocode))
    if (cumsumpos)
        fg <- fg_feat_cumsum(fg, no_cores=no_cores)

    if (prop) { # included in specenr
        fg <- fg_feat_node_prop(fg)
        fg <- fg_feat_edge_prop(fg, no_cores=no_cores)
    }

    if (specenr) {
        fg <- fg_feat_node_specenr(fg, no_cores=no_cores)
    }

    try({
        if (calculate_summary &
            ifelse(length(class)==1, class%in%colnames(meta), TRUE) &
            all(node_features%in%c(names(fg_get_feature_all(fg)$node), "NONE")) &
            all(edge_features%in%c(names(fg_get_feature_all(fg)$edge), "NONE"))) {

            time_output(start1)
            start1 <- Sys.time()
            message("calculating summary statistics")

            fg <- fg_summary(
                fg, no_cores=no_cores,
                class=ifelse(length(class)==1, class, "class"),
                label1=label1, label2=label2,
                node_features=node_features,
                edge_features=edge_features,
                overwrite=FALSE,
                test_custom=test_custom,
                test_name=test_name,
                diminish=diminish)
        }
    })
    time_output(start1)

    saved <- FALSE
    if (!is.null(path))
        tryCatch({
            saved <- fg_save(fg, path, save_plots=FALSE)
        }, error=function(e) {
            message("flowGraph object not saved; check path and try again with fg_save!")
        })

    if (calculate_summary & save_plots & saved)
        fg_save_plots(fg, plot_path=paste0(path, "/plots"))

    time_output(start, "total time used")
    return(fg)
}
