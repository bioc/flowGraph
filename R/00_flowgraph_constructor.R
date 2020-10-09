#' @name flowGraph
#' @title flowGraph object constructor.
#' @description Initializes a \code{flowGraph} object given the cell counts or
#'  Phenotypes objects, from the flowType package, for one or more
#'  flow cytometry sample(s). The flowGraph object returned holds meta data
#'  for each sample, each cell population node,
#'  edges representing how each cell population node relate to one another,
#'  and features for these nodes and edges.
#' @param input_ Any of the following:
#' \itemize{
#'   \item{a \code{Phenotypes} object from the \code{flowType} package.}
#'   \item{a numeric matrix or vector of the cell counts from a
#'   \code{Phenotypes} object and its column/names must be the phenotype names.}
#'   \item{a list of \code{Phenotypes} object(s).}
#'   \item{a vector of \code{Phenotypes} object full paths.}
#' }
#' All input samples should have the same \code{markers} and
#' \code{partitionsPerMarker}. If given \code{Phenotype} object(s),
#' the function will look to generate phenotype names from the
#' \code{Phenotypes@@PhenoCodes} or \code{rownames(Phenotypes@@MFIs)} slots.
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
#'  object should save its elements.
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
#' @param normalize A logical variable indicating whether or not to
#'  calculate normalized cell count;
#'  this can be done later on with \code{fg_feat_node_norm}.
#' @param norm_ind See \code{\link[flowGraph]{fg_feat_node_norm}}.
#' @param norm_layer See \code{\link[flowGraph]{fg_feat_node_norm}}.
#' @param norm_path See \code{\link[flowGraph]{fg_feat_node_norm}}.
#' @param save_plots A logical indicating whether or not to save plots.
#' @param ... Other parameters used in the \code{fg_save_plots} function
#'  of the \code{edgeR} package. See \code{\link[flowGraph]{tmm}}.
#' @return flowGraph object
#' @details \code{flowGraph} is the constructor for the \code{flowGraph} object.
#'  The user can choose to input as \code{input_} a vector, a \code{Phenotypes}
#'  object (meaning there is only one sample), a matrix, or a \code{Phenotypes}
#'  object list. If the user is also inputting a sample meta data frame, it
#'  must contain a \code{id} column corresponding to sample names.
#' @examples
#'
#' ## flowType package has been deprecated,
#' ## Phenotypes input format will be restored after flowTypeFilter is uploaded.
#' # library(flowType)
#' #
#' # # prepare parallel backend
#' # no_cores <- 1 #parallel::detectCores()-1
#' # # future::plan(future::multiprocess)
#'
#' # ## create Phenotypes data ----------------------
#'
#' # # define marker and total cell count
#' # celln <- 10000
#' # markern <- 3
#' # markers <- LETTERS[1:markern]
#'
#' # # define marker thresholds
#' # cvd <- rnorm(celln,2,1)
#' # p50 <- quantile(cvd, .5)
#' # thres <- lapply(markers, function(x) p50)
#' # names(thres) <- markers
#'
#' # # generate flowType Phenotypes list
#' # samplen <- 10
#' # # ftl <- furrr::future_map(1:samplen, function(i) {
#' # ftl <- lapply(1:samplen, function(i) {
#' #     # make flow frame
#' #     f <- new("flowFrame")
#' #     f@exprs <- matrix(rnorm(celln*markern,2,1), nrow=celln)
#' #     colnames(f@exprs) <- markers
#' #
#' #     # marker indices in flow frame
#' #     ci <- c(1:ncol(f@exprs))
#' #     names(ci) <- colnames(f@exprs)
#' #
#' #     # modify experiment samples such that ABC increases by 50%
#' #     if (i>(samplen/2)) {
#' #         ap <- f@exprs[,1]>thres[[1]]
#' #         bp <- f@exprs[,2]>thres[[2]]
#' #         cp <- f@exprs[,3]>thres[[3]]
#' #         tm <- sum(ap & bp & cp)/2
#' #         f@exprs <- rbind(f@exprs,
#' #                          f@exprs[sample(which(ap & bp & cp),tm),])
#' #     }
#' #
#' #     # make flowType Phenotypes
#' #     flowType(Frame=f, PropMarkers=ci, MarkerNames=colnames(f@exprs),
#' #              MaxMarkersPerPop=markern, PartitionsPerMarker=2,
#' #              Thresholds=thres,
#' #              Methods='Thresholds', verbose=FALSE, MemLimit=60)
#' # })
#'
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
#' ## flowType package is deprecated,
#' ## Phenotypes input format will be restored after flowTypeFilter is uploaded.
#' # # input: Phenotype object
#' # fg <- flowGraph(ftl[[1]], no_cores=no_cores)
#'
#' # # input: Phenotype list
#' # fg <- flowGraph(ftl, class=meta_file$class, no_cores=no_cores)
#' # # fg <- flowGraph(ftl, meta=meta_file, no_cores=no_cores)
#'
#' # input: vector of load-able Phenotypes paths
#' fg <- flowGraph(fg_data_pos30$count[1,], no_cores=no_cores)
#'
#' # input: matrix + vector of class corresponding to samples
#' fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                 no_cores=no_cores)
#'
#' # input: matrix + meta data frame
#' # fg <- flowGraph(fg_data_pos30$count, meta=fg_data_pos30$meta,
#' #                 no_cores=no_cores)
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
# #' @importFrom flowType decodePhenotype
#' @importFrom Matrix Matrix
#' @importFrom methods new
#' @importFrom igraph layout.reingold.tilford
flowGraph <- function(
    input_, meta=NULL, class="class", no_cores=1, markers=NULL,
    layout_fun="layout.reingold.tilford", # layout.circle

    cumsumpos=FALSE, # whether to make positives +=+/++ (cumsum)
    prop=TRUE, specenr=TRUE,

    path=NULL,

    # summary parameters
    calculate_summary=TRUE, node_features="SpecEnr", edge_features="NONE",
    test_name="t_diminish",
    test_custom="t",
    diminish=TRUE, label1=NULL, label2=NULL,

    # normalize parameters
    normalize=FALSE,
    norm_ind=0, norm_layer=3, norm_path=NULL,

    # plotting parameters
    save_plots=TRUE, ...

) {
    options(stringsAsFactors=FALSE)

    if (!base::is.null(meta))
        if (!"id"%in%base::colnames(meta))
            stop("meta must have a column named \"id\" containing sample ID's")

    start <- start1 <- Sys.time()
    message("preparing input; ")

    # convert input_ into a cell count matrix ----------------------------------
    msg <- "make sure input type is correct"

    phenocode <- NULL
    if (any(class(input_)%in%c("integer","numeric")) |
        grepl("matrix",class(input_), ignore.case=TRUE)) {
        ## feature: count (sample x cell population)
        # mc <- flowGraph_matrix(input_)
        if (base::is.null(base::dim(input_))) {
            mc <- mc0 <- base::matrix(input_,nrow=1)
            base::colnames(mc) <- base::names(input_)
            base::rownames(mc) <- "s1"
        } else {
            mc <- mc0 <- input_
            if (base::is.null(base::rownames(input_)))
                base::rownames(input_) <-
                    paste0("s",base::seq_len(nrow(input_)))
            base::rownames(mc) <- base::rownames(input_)
        }

        ## make meta for samples
        if (base::is.null(meta)) sample_id <- base::rownames(mc)

        ## make meta for cell populations
        if (!any(base::grepl("[+-]",base::colnames(mc)))) stop(msg)

        phen <- base::colnames(mc)
        if (base::is.null(markers)) {
            markers <-
                base::unique(base::unlist(stringr::str_split(phen,"[-+]+")))
            markers <- markers[markers!=""]
        }

    } else {
        stop("input must be a numeric sample x cell population matrix.")
        # if (is(input_, "Phenotypes"))
        #     input_ <- base::list(s1=input_)
        # if (is(input_,"list")) {
        #     testclass <- purrr::map_lgl(input_, is, "Phenotypes")
        #     if (!all(testclass)) stop(msg)
        #     ftl <- input_
        # } else if (is(input_, "character")) {
        #     # ftl <- flowGraph_load_ftl(input_)
        #     no_cores <- flowGraph:::ncores(no_cores)
        #     if (no_cores>1) future::plan(future::multiprocess)
        #
        #     ftl <- furrr::future_map(input_, function(x) base::get(load(x)))
        #     base::names(ftl) <- purrr::map_chr(stringr::str_split(input_,"/"),
        #                                        function(x) x[base::length(x)])
        # } else {
        #     stop(msg)
        # }
        #
        # # make sample_id
        # if (base::is.null(meta)) {
        #     sample_id <- base::names(ftl)
        #     if (base::is.null(sample_id))
        #         sample_id = paste0("s", seq_len(base::length(ftl)))
        # }
        #
        # ## make meta for cell populations
        # if (base::is.null(markers))
        #     markers <- ftl[[1]]@MarkerNames
        # phen <- NULL
        # try({
        #     phen <- purrr::map_chr(ftl[[1]]@PhenoCodes, function(x)
        #         flowType::decodePhenotype(x, markers,
        #                                   ftl[[1]]@PartitionsPerMarker))
        # }, silent=TRUE)
        # if (base::is.null(phen))
        #     try({ phen <- base::rownames(ftl[[1]]@MFIs) }, silent=TRUE)
        # if (base::is.null(phen))
        #     stop("no phenotype cell populations labels in Phenotype file.")
        #
        # try ({ phenocode <- ftl[[1]]@PhenoCodes }, silent=TRUE)
        #
        # ## feature: count (sample x cell population)
        # mc <- as.matrix(base::do.call(rbind,purrr::map(
        #     ftl, function(ft) ft@CellFreqs)))
        # if (is(mc[1],"character")){ # some versions of purrr might give char
        #     mc <- as.matrix(mc[,-1])
        #     class(mc) <- "numeric"
        # }
        # # if (base::is.null(dim(mc)[1])) mc <- matrix(mc, nrow=1)
        #
        # # remove Phenotype list to save memory
        # rm(ftl); gc()
    }

    time_output(start1)

    # prepare sample meta ----------------------------------
    if (!base::is.null(meta)) {
        if (base::nrow(meta)!= base::nrow(mc))
            stop("meta data and cell count matrix have different row counts")
        if (!"id"%in%base::colnames(meta))
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
        if (base::length(class)==base::length(sample_id)) {
            if ("class"%in%colnames(meta)) {
                meta$class_ <- class
            } else {
                meta$class <- class
            }
        }
    }

    ## feature: count (sample x cell population) -----------
    start1 <- Sys.time()
    message("preparing feature(s): mapping out cell population relations")

    keepinds <- base::apply(mc, 2, function(x) any(x>0))
    keepinds <- base::apply(mc, 2, function(x) any(x>0))
    mc <- mc[,keepinds,drop=FALSE]
    base::colnames(mc) <- phen[keepinds]
    base::rownames(mc) <- meta$id
    mc <- Matrix::Matrix(mc, sparse=TRUE) # rm all 0 cols
    if (!base::is.null(phenocode)) phenocode <- phenocode[keepinds]

    ## make meta for cell populations FINAL -----------------
    meta_cell <- get_phen_meta(phen[keepinds],phenocode)

    # make parent list (for each cell popultion, list its parents)
    pccell <- get_phen_list(meta_cell=meta_cell, no_cores=no_cores)
    edf <- pccell$edf
    pchild <- pccell$pchild # not used, returned
    pparen <- pccell$pparen

    # extract cell populations that have parents to calculate
    # expected proportions for; just to be safe
    cells1 <- meta_cell$phenotype[meta_cell$phenolayer==1]
    cells1 <- base::append("",cells1[cells1%in%base::names(pparen)])
    cells <- meta_cell$phenotype[meta_cell$phenolayer>1]
    cells <- cells[cells%in%base::names(pparen)]
    cells_ <- base::append(cells1,cells)

    meta_cell <- meta_cell[base::match(cells_,meta_cell$phenotype),]
    rooti <- base::which(meta_cell$phenolayer==0)

    mc <- mc[,base::match(cells_, base::colnames(mc)),drop=FALSE] # final cpops

    # trim/order list of parents and children
    pparen_ <- pparen[base::names(pparen)%in%cells_]
    pparen_[cells] <- purrr::map(pparen_[cells], function(x) {
        a <- x[x%in%cells_]
        if (base::length(a)==0) return(NULL)
        a
    })
    pparen_ <- purrr::compact(pparen_)

    pchild_ <- pchild[base::names(pchild)%in%cells_]
    pchild_ <- purrr::map(pchild_, function(x) {
        a <- x
        for (i in base::seq_len(base::length(a)))
            a[[i]] <- x[[i]][x[[i]]%in%cells_]
        if (all(purrr::map_int(a, base::length)==0)) return(NULL)
        a
    })
    pchild_ <- purrr::compact(pchild_)

    edf <- edf[edf$to%in%cells_ & edf$from%in%cells_,,drop=FALSE]

    time_output(start1, "cell population meta data created")


    ## list of outputs
    gr <- base::list(e=edf, v=meta_cell)
    gr <- set_layout_graph(gr, layout_fun) # layout cell hierarchy

    if (base::is.null(meta)) meta <- base::data.frame(id=sample_id)

    desc <- summary_table(mc,"count")

    fg <- methods::new(
        "flowGraph",
        feat=base::list(node=base::list(count=mc), edge=base::list()),
        feat_desc=base::list(node=desc),
        markers=markers,
        graph=gr, edge_list=base::list(child=pchild_,parent=pparen_),
        meta=meta, plot_layout=base::as.character(substitute(layout_fun)),
        etc=base::list(cumsumpos=FALSE, class_mean_normalized=FALSE)
    )

    ## features ----------------------------
    cumsumpos <- cumsumpos & any(base::grepl("3",meta_cell$phenocode))
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
            ifelse(base::length(class)==1, class%in%colnames(meta), TRUE) &
            all(node_features%in%c(names(fg@feat$node), "NONE")) &
            all(edge_features%in%c(names(fg@feat$edge), "NONE"))) {
            fg <- fg_summary(
                fg, no_cores=no_cores,
                class=ifelse(base::length(class)==1, class, "class"),
                label1=label1, label2=label2,
                node_features=node_features,
                edge_features=edge_features,
                overwrite=FALSE,
                test_custom=test_custom,
                test_name=test_name,
                diminish=diminish)
        }
    })

    fg@etc$save <- list(id=stringi::stri_rand_strings(1,5))
    saved <- FALSE
    if (!base::is.null(path))
        tryCatch({
            fg_save(fg, path, save_plots=FALSE)
            saved <- TRUE
            fg@etc$save$path <- path
        }, error=function(e) {
            message("flowGraph object not saved; check path and try again with fg_save!")
        })

    if (normalize) {
        if (is.null(norm_path)) norm_path <- paste0(path,"/Count_normalization")
        fg <- fg_feat_node_norm(fg, norm_ind=norm_ind, norm_layer=norm_layer,
                                norm_path=norm_path, no_cores=no_cores)

        if (calculate_summary & ("count_norm"%in%node_features |
                                 base::is.null(node_features)))
            fg <- fg_summary(
                fg, no_cores=no_cores, class=class,
                label1=label1, label2=label2,
                node_features="count_norm",
                edge_features="NONE",
                overwrite=FALSE,
                test_custom=test_custom,
                test_name=test_name,
                diminish=diminish)

        if (saved) {
            fn_dir <- paste0(path, "/Features/Nodes")
            fn <- "count_norm"
            m <- as.matrix(fg@feat$node[[fn]])
            write.csv(m, file=paste0(fn_dir,"/", fn, ".csv"))
            sm <- fg@feat_desc$node
            write.csv(sm, file=paste0(fn_dir,".csv"), row.names=FALSE)
        }
    }

    if (calculate_summary & save_plots & saved)
        fg_save_plots(fg, plot_path=paste0(path, "/plots"), ...)

    time_output(start, "total time used")
    return(fg)
}


#' ## input_ processing: matrix
#' #' @title Cleans input cell count matrix.
#' #' @description Cleans input cell count matrix.
#' #' @param input_ A numeric matrix.
#' #' @return A numeric matrix
#' #' @rdname flowGraph_matrix
#' flowGraph_matrix <- function(input_) {
#'     if (base::is.null(base::dim(input_))) {
#'         mc <- mc0 <- base::matrix(input_,nrow=1)
#'         base::colnames(mc) <- base::names(input_)
#'         base::rownames(mc) <- "s1"
#'     } else {
#'         mc <- mc0 <- input_
#'         if (base::is.null(base::rownames(input_)))
#'             base::rownames(input_) <-
#'                 paste0("s",base::seq_len(nrow(input_)))
#'         base::rownames(mc) <- base::rownames(input_)
#'     }
#'     return(mc)
#' }
#'
#' ## input_ processing: loading flowtypes from paths
#' #' @title Loads Phenotypes objects from vector of paths.
#' #' @description Loads Phenotypes objects from vector of paths.
#' #' @param input_ A string vector of \code{load()}-able Phenotypes
#' #' object paths.
#' #' @param no_cores An integer indicating number of cores to parallelize on.
#' #' @return A list of Phenotypes objects.
#' #' @rdname flowGraph_load_ftl
#' #' @importFrom future plan multiprocess
#' #' @importFrom furrr future_map
#' flowGraph_load_ftl <- function(input_, no_cores=1) {
#'     no_cores <- flowGraph:::ncores(no_cores)
#'     if (no_cores>1) future::plan(future::multiprocess)
#'
#'     ftl <- furrr::future_map(input_, function(x) base::get(load(x)))
#'
#'     return(ftl)
#' }
