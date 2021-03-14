#' @title Wrapper for map
#' @description Wrapper for \code{purrr::map} and \code{furrr::future_map}
#'  to handle parallel-ization
#' @param x Variable to recurse over; must be indices!
#' @param f Function to recurse over.
#' @param no_cores Number of cores to use; \code{future} must have
#'  been ran already.
#' @param prll If set to FALSE, forces use of purrr::map instead of
#' furrr::future_map, Default: TRUE
#' @param ... Other parameters used by \code{f}.
#' @return Unnested named list.
#' @details Wrapper for \code{purrr::map} and \code{furrr::future_map} to handle
#' parallel-ization easily; note that \code{future} must have been ran already
#' outside of the function and outputs will always be a list.
#' @seealso
#'  \code{\link[purrr]{map}}
#'  \code{\link[furrr]{future_map}}
#' @rdname fpurrr_map
# #' @export
#' @importFrom purrr map
#' @importFrom furrr future_map
fpurrr_map <- function(x, f, no_cores=1, prll=TRUE, ...) {
    if (length(x)==0) return(NULL)
    if (no_cores==1 | !prll) {
        res <- purrr::map(x, f)
    } else if (length(x)/no_cores < 2) {
        res <- furrr::future_map(x, f)
    } else {
        xxx <- loop_ind_f(x, no_cores)
        res <- furrr::future_map(xxx, function(xx) {
            a <- purrr::map(xx, f)
            if (is.list(a)) return(a)
            list(a)
        })
        res <- unlist(res, recursive=FALSE)
    }
    if (is.character(x))
        names(res) <- x
    return(res)
}


#' @title Extracts markers from cell population phenotypes
#' @description Extracts all unique markers from cell population phenotypes
#' @param phen A vector of cell population phenotypes.
#' @return A vector of unique markers
#' @seealso
#'  \code{\link[stringr]{str_split}}
#' @rdname extract_markers
# #' @export
#' @importFrom stringr str_split
extract_markers <- function(phen) {
    phen <- gsub("_not","_",phen)
    markers <- unique(unlist(stringr::str_split(phen,"[-+_]+")))
    markers[markers!=""]
}


#' @title Reformats phenotype
#' @description Reformats cell population phenotypes into flowGraph format
#' @param phen Vector of cell population phenotype names as character strings.
#' @param markers markers extracted from \code{phen}.
#' @return Vector with the same length as \code{phen} containing reformatted
#' and not necessarily changed cell population phenotype names.
#' @seealso
#'  \code{\link[stringr]{str_extract}},\code{\link[stringr]{str_split}}
#' @examples
#'
#' # fg_clean_phen(c("A+_B+","B+_notC","A-_C"))
#'
#' @rdname fg_clean_phen
# #' @export
#' @importFrom stringr str_extract str_split
fg_clean_phen <- function(phen, markers=NULL) {
    if (is.null(markers))
        markers <- extract_markers(phen)
    if (any(grepl("[_]not",phen))) {
        gate_markers <- gsub("not","",unique(stringr::str_extract(
            unlist(stringr::str_split(phen,"[-+_]+")),
            paste0("^not", paste0(markers),collapse="|") )))
        gate_markers <- gate_markers[!is.na(gate_markers)]
        for (gm in gate_markers) {
            phen <- gsub(paste0("^",gm,"$|_",gm,"$|_",gm,"_|^",gm), paste0("_",gm,"+"), phen)
            phen <- gsub(paste0("not",gm), paste0(gm,"-"), phen)
            phen <- gsub("^_", "", phen)
        }
    }
    if (any(grepl("[_]",phen)))
        phen <- gsub("_","",phen)
    phen
}

#' @title Converts input into a significance test function
#' @param test_custom a string \code{c("t", "wilcox","ks","var","chisq")} or a
#'  function.
#' @return a statistical significance test function.
#' @seealso
#'  \code{\link[stats]{t.test}},\code{\link[stats]{wilcox.test}},
#'  \code{\link[stats]{ks.test}},\code{\link[stats]{var.test}},
#'  \code{\link[stats]{chisq.test}}
#' @rdname test_c
# #' @export
#' @importFrom stats t.test wilcox.test ks.test var.test chisq.test
test_c <- function(test_custom) {
    if (!is.function(test_custom))
        if (test_custom=="t") {
            test_custom <- function(x,y)
                tryCatch(stats::t.test(x,y)$p.value,error=function(e)1)
        } else if (test_custom=="wilcox") {
            test_custom <- function(x,y)
                tryCatch(stats::wilcox.test(x,y)$p.value,error=function(e)1)
        } else if (test_custom=="ks") {
            test_custom <- function(x,y)
                tryCatch(stats::ks.test(x,y)$p.value,error=function(e)1)
        } else if (test_custom=="var") {
            test_custom <- function(x,y)
                tryCatch(stats::var.test(x,y)$p.value,error=function(e)1)
        } else if (test_custom=="chisq") {
            test_custom <- function(x,y)
                tryCatch(stats::chisq.test(x,y)$p.value,error=function(e)1)
        }
    test_custom
}


#' @title Calcuate SpecEnr from proportion and expected proportion
#' @description FUNCTION_DESCRIPTION
#' @param mp_ Numerical sample x cell population matrix w/ proportions.
#' @param me_ Numerical sample x cell population matrix w/ expected proportions.
#' @return Numerical sample x cell population matrix w/ SpecEnr.
# #' @rdname ms_create
# #' @export
ms_create <- function(mp_, me_) {
    mp_ <- as.matrix(mp_)
    me_ <- as.matrix(me_)
    ms_ <- as.matrix(mp_/me_)
    if (is.na(dim(ms_)[1]))
        ms_ <- matrix(ms_, ncol=ncol(mp_))
    ms_nan <- is.nan(ms_)
    ms_[is.infinite(ms_)] <- max(ms_[is.finite(ms_)])
    ms_[ms_nan] <- 0
    e0 <- as.matrix(me_)==0
    suppressWarnings({ ms_ <- as.matrix(log(ms_)) })
    ms_[is.nan(ms_)] <- 0
    ms_[e0] <- log(as.matrix(mp_)[e0])
    ms_[mp_==0] <- 0

    dimnames(ms_) <- dimnames(mp_)
    return(ms_)
}


#' @title Gets child populations of given cell populations
#' @description Gets the child populations of a vector of given cell populations
#'  \code{parens} and updates \code{pchild} the edge list if edge list doesn't contain
#'  the requested information.
#' @param parens Character vector of cell population phenotypes.
#' @param pchild Edge list where the name of the list is the cell population
#'  and the vector in each element contains the child cell populations of the
#'  named cell population.
#' @param pc_i A cell population x marker matrix where the values are 0/1/2/...
#'  correspondng to marker conditions /-/+/... for possible PARENT populations.
#' @param ac__ A list where the elements are marker index > "0"/"1"/"2"/... >
#'  a logical vector the same length as the number of cell population phenotypes
#'  indicating whether or not the marker condition exists in them; this is for
#'  the possible CHILD cell populations
#' @param meta_cell__ data frame with meta data for cell population phenotypes
#'  from the flowGraph object; this is for the possible CHILD cell populations.
#' @return A list containing child populations of \code{parens}; also globally
#'  updates \code{pchild}.
#' @seealso
#'  \code{\link[purrr]{map}},\code{\link[purrr]{keep}}
# #' @rdname get_child
# #' @export
#' @importFrom purrr map compact
get_child <- function(parens, pchild, pc_i, ac__, meta_cell__) {
    parens_ <- parens[!parens%in%names(pchild)]
    if (length(parens_)==0)
        return(pchild[match(parens,names(pchild))])

    res <- fpurrr_map(seq_len(nrow(pc_i)), function(j) {
        mcgrow <- pc_i[j,]
        colj1 <- which(mcgrow > 0)
        chi <- Reduce("&", purrr::map(colj1, function(coli)
            ac__[[coli]][[as.character(mcgrow[coli])]] ))
        meta_cell__$phenotype[chi]
    }, no_cores, prll=nrow(pc_i)>1000)#, pc=pc_i, ac=ac__, mcell=meta_cell__)
    names(res) <- rownames(pc_i)

    pchild <- append(pchild, purrr::compact(res))
    pchild <<- pchild
    return(pchild[match(parens,names(pchild))])
}

#' @title Gets parent populations of given cell populations
#' @description Gets the parent populations of a vector of given cell
#'  \code{childs} and updates \code{pparen} the edge list if edge list doesn't
#'  contain the requested information.
#' @param childs Character vector of cell population phenotypes.
#' @param pparen Edge list where the name of the list is the cell population
#'  and the vector in each element contains the parent cell populations of the
#'  named cell population.
#' @param pc__i A cell population x marker matrix where the values are 0/1/2/...
#'  correspondng to marker conditions /-/+/...;
#'  this is for the possible CHILD cell populations.
#' @param ac_ A list where the elements are marker index > "0"/"1"/"2"/... >
#'  a logical vector the same length as the number of cell population phenotypes
#'  indicating whether or not the marker condition exists in them;
#'  this is for the possible PARENT cell populations.
#' @param meta_cell_ data frame with meta data for cell population phenotypes
#'  from the flowGraph object; this is for the possible PARENT cell populations.
#' @return A list containing parent populations of \code{childs}; also globally
#'  updates \code{pparen}.
#' @rdname get_paren
# #' @export
#' @importFrom purrr map compact
get_paren <- function(childs, pparen, pc__i, ac_, meta_cell_) {
    childs_ <- childs[!childs%in%names(pparen)]
    if (length(childs_)==0)
        return(pparen[match(childs,names(pparen))])

    res <- fpurrr_map(seq_len(nrow(pc__i)), function(j) {
        mcgrow <- pc__i[j, ]
        colj1 <- which(mcgrow > 0)
        chidf <- do.call(cbind, purrr::map(colj1, function(coli)
            ac_[[coli]][[as.character(mcgrow[coli])]] ))
        chi <- apply(chidf, 1, function(x) sum(!x) == 1)
        meta_cell_$phenotype[chi]
    }, no_cores, prll=nrow(pc__i)>1000)
    names(res) <- rownames(pc__i)

    pparen <- append(pparen, purrr::compact(res))
    pparen <<- pparen
    return(pparen[match(childs,names(pparen))])
}

#' @title Gets edge proportions of a given edge matrix
#' @description Gets the edge proportions of the edges in edge matrix
#'  \code{edf_} and updates \code{ep} edge proportion matrix if it didn't
#'  contain the requested information.
#' @param edf_ edge x from&to data frame containing edges and their from and to
#'  cell population phenotypes.
#' @param ep sample x edge (parent_child) matrix with edge proportions.
#' @param mp_ sample x phenotype matrix with proportions.
#' @param no_cores Number of cores to use, Default: 1
#' @return \code{ep} with only the specific columns (edges) requested; also updates
#' \code{ep} globally.
#' @rdname get_eprop
# #' @export
get_eprop <- function(edf_, ep, mp_, no_cores=1) {
    clnm <- apply(edf_, 1, paste0, collapse="_")
    todoi <- which(!clnm%in%colnames(ep))
    if (length(todoi)==0)
        return(ep[,clnm,drop=FALSE])
    ep_ <- do.call(cbind, fpurrr_map(todoi, function(i) {
        if (edf_$from[i]=="") return(mp_[,edf_$to[i]])
        mp_[,edf_$to[i]]/mp_[,edf_$from[i]]
    }, no_cores, prll=length(todoi)>1000))
    colnames(ep_) <- clnm[todoi]
    if (is.null(ep))
        rownames(ep_) <- rownames(mp_)
    ep <- cbind(ep, ep_)
    ep <<- ep

    return(ep[,clnm,drop=FALSE])
}


#' @title Determines which phenotypes are statistically significant
#' @description Determines which phenotypes are statistically significant based
#'  on SpecEnr.
#' @param ms_ sample x phenotype SpecEnr matrix
#' @param summary_pars See \code{flowGraphSubset}.
#' @param summary_adjust See \code{flowGraphSubset}.
#' @param test_cust Final significance test function.
#' @param test_custom Raw significance test function.
#' @param lyrno An integer indicating total number of layers in the cell
#'  hierarchy including layer 0.
#' @param mp_ sample x phenotype proportion matrix.
#' @param me_ sample x phenotype expected proportion matrix.
#' @return A logical vector the same length as the number of columns in
#'  \code{ms_} indicating whether or not each phenotype is significant;
#'  used only for the fast version of flowGraph to determine whether or not
#'  to keep testing the phenotypes' children.
#' @rdname ms_psig
# #' @export
ms_psig <- function(ms_, summary_pars, summary_adjust,
                    test_cust, test_custom, lyrno, mp_, me_) {
    if (is.null(summary_adjust$filter_btwn_tpthres)) {
        warning("no p threshold specified in `summary_adjust$filter_btwn_tpthres`, aplying default threshold 0.05")
        summary_adjust$filter_btwn_tpthres <- .05
    }
    btwn_test_custom <- test_c(summary_adjust$btwn_test_custom)

    ## calculate p values
    p <- apply(ms_,2,test_cust,meta[[summary_pars$class]],test_custom)
    names(p) <- colnames(ms_)

    ## adjust p values
    dual <- FALSE
    if (!is.null(summary_pars$labels) & !is.null(summary_adjust))
        if (length(summary_pars$labels)==2)
            dual <- TRUE

    if (summary_adjust$adjust_custom=="byLayer")
        p <- sapply(p*ncol(ms_)*lyrno, min, 1)
    p[is.na(p)] <- 1
    ptf <- rep(FALSE, length(p))
    sigcands <- rep(TRUE, length(p))

    # cohen's d, number of 0s, pvalue in between
    if (dual & !is.null(btwn_test_custom) &
        summary_adjust$filter_btwn_es>0 |
        summary_adjust$filter_btwn_tpthres<1) {
        id1 <- summary_pars$labels[[1]]
        id2 <- summary_pars$labels[[2]]

        pp3v1 <- sapply(seq_len(ncol(mp_)), function(ci)
            btwn_test_custom(mp_[id1,ci], me_[id1,ci]) )
        pp3v2 <- sapply(seq_len(ncol(mp_)), function(ci)
            btwn_test_custom(mp_[id2,ci], me_[id2,ci]) )

        pp3b_raw <- lapply(seq_len(ncol(mp_)), function(ci) {
            a <- mp_[id1,ci] - mp_[id2,ci]
            b <- me_[id1,ci] - me_[id2,ci]
            return(list(a=a,b=b))
        })
        pp3bt <- sapply(seq_len(ncol(mp_)), function(ci)
            btwn_test_custom( pp3b_raw[[ci]]$a, pp3b_raw[[ci]]$b) )
        pp3bt[is.na(pp3bt)] <- 1

        sig1 <- pp3v1<summary_adjust$filter_btwn_tpthres
        sig2 <- pp3v2<summary_adjust$filter_btwn_tpthres
        sig3 <- pp3bt<summary_adjust$filter_btwn_tpthres
        sigcands <- (sig1 | sig2) & sig3
    }
    ptf[sigcands & p<summary_adjust$filter_btwn_tpthres] <- TRUE
    return(ptf)
}


#' @name flowGraphSubset_summary_pars
#' @title Default for flowGraphSubset's summary_pars
#' @description Default input for flowGraphSubset's \code{summary_pars} parameter.
#' @return Default list parameter flowGraphSubset's \code{summary_pars} parameter.
#' @rdname flowGraphSubset_summary_pars
#' @examples
#'
#'  flowGraphSubset_summary_pars()
#'
#' @export
flowGraphSubset_summary_pars <- function()
    list(
        node_feature="SpecEnr",
        edge_feature="NONE",
        test_name="t_diminish",
        test_custom="t",
        diminish=TRUE,
        class="class",
        labels=NULL)


#' @name flowGraphSubset_summary_adjust
#' @title Default for flowGraphSubset's summary_adjust
#' @description Default input for flowGraphSubset's \code{summary_adjust}
#'  parameter. ONLY USE THIS OVER flowGraph IF: 1) your data set has more than
#'  10,000 cell populations and you want to speed up your calculation time AND
#'  2) you only have one set of classes you want to test on the
#'  SAME SET OF SAMPLES (e.g. control vs experiment). As flowGraphSubset doesn't
#'  calculate the SpecEnr for all cell populations, so if you want to test other
#'  sets of classes on the same sample, you will not be able to test all
#'  possible cell populations on the new set of classes.
#' @return Default list parameter flowGraphSubset's \code{summary_adjust} parameter.
#' @rdname flowGraphSubset_summary_adjust
#' @examples
#'
#'  flowGraphSubset_summary_adjust()
#'
#' @export
flowGraphSubset_summary_adjust <- function()
    list(
        adjust_custom="byLayer", btwn_test_custom="t",
        adjust0_lim=c(-.1,.1), filter_adjust0=1, filter_es=0,
        filter_btwn_tpthres=.05, filter_btwn_es=.5)


#' @name flowGraphSubset
#' @title flowGraph object constructor.
#' @description Initializes a \code{flowGraph} object given the cell counts
#'  for one or more flow cytometry sample(s).
#'  The flowGraph object returned holds meta data
#'  for each sample, each cell population node,
#'  edges representing how each cell population node relate to one another,
#'  and features for these nodes and edges.
#' @param input_ a numeric matrix of the cell counts;
#'   its column/names must be the phenotype names and its rownames
#'   must be sample ID's.
#' @param meta A data frame with meta data for each \code{Phenotypes} or sample;
#'  One of its column names should be "id" whose values correspond to
#'  the name of each \code{Phenotypes} object. We also recommend for it to have
#'  a column named "class" where one of its unique values is "control".
#' @param class A string corresponding to the column name or index
#'  of \code{meta} whose values represent
#'  the class label of each sample, Default: 'class'
#' @param no_cores An integer indicating how many cores to parallelize on,
#'  Default: 1
#' @param markers A string vector of marker names used in \code{input_},
#'  Default: NULL
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
#'  A+, A++, A+++, ..., or A++, A+++, ... , Default: FALSE
#' @param path A string indicating the folder path to where the flowGraph
#'  object should save its elements, Default = NULL (don't save).
#' @param summary_pars A list containing parameters for calculating the
#'  statistical significance summary significance that will determine whether
#'  to trim out phenotypes for this fast version of flowGraph. The lists'
#'  elements are:
#' \itemize{
#'   \item{\code{node_feature}: "SpecEnr"; this is the feature we will be
#'    testing, don't change this.}
#'   \item{\code{edge_feature}: "NONE"; this unneeded for now.}
#'   \item{\code{test_name}: "t_diminish"; this unneeded for now.}
#'   \item{\code{test_custom}: "t"; a string or a function indicating the
#'    statistical test desires. These tests can be
#'    \code{c("t", "wilcox","ks","var","chisq")} corresponding to functions
#'    \code{\link[stats]{t.test}},\code{\link[stats]{wilcox.test}},
#'    \code{\link[stats]{ks.test}},\code{\link[stats]{var.test}},
#'    \code{\link[stats]{chisq.test}}}
#'   \item{\code{diminish}: TRUE; whether or not to continue testing
#'    phenotypes whos parent phenotypes are all insignificant.}
#'   \item{\code{class}: "class"; the column name in \code{meta} that contains
#'    class labels you want to test.}
#'   \item{\code{labels}: c("aml", "control") for the flowcap data set;
#'   SET THIS!! to the class labels you want to test using \code{test_custom}.}
#' }
#' @param summary_adjust A list of parameters on how to adjust the p-values;
#'  this also affects which phenotypes are tested. The elements in the list are:
#' \itemize{
#'   \item{\code{adjust_custom}: "byLayer"; this is a string (corresponding to
#'    an option in \code{\link[stats]{p.adjust}}) or a function used to
#'    adjust p-values.}
#'   \item{\code{btwn_test_custom}: "t"; see \code{test_custom} in
#'    \code{summary_pars}; this statistical significance test is used in
#'    the filters.}
#'   \item{\code{adjust0_lim}: see \code{fg_get_summary}.}
#'   \item{\code{filter_adjust0}: see \code{fg_get_summary}.}
#'   \item{\code{filter_es}: see \code{fg_get_summary}.}
#'   \item{\code{filter_btwn_tpthres}: see \code{fg_get_summary}.}
#'   \item{\code{filter_btwn_es}: see \code{fg_get_summary}.}
#' }
#' @param save_plots A logical indicating whether or not to save plots.
#' @return flowGraph object
#' @details All node and edge features are trimmed such that only the
#'  significant phenotypes are left; the original input is stored in the slot
#'  \code{etc$original_count} of the returned flowGraph object.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   data(fg_data_pos2)
#'   fg <- flowGraph(fg_data_pos2$count, meta=fg_data_pos2$meta, no_cores=1)
#'  }
#' }
#' @rdname flowGraphSubset
#' @export
#' @importFrom Matrix Matrix
#' @importFrom purrr map_lgl map map_dfr compact
#' @importFrom methods new
#' @importFrom future plan multiprocess
#' @importFrom data.table as.data.table setattr ".I"
#' @importFrom stringr str_split
#' @importFrom matrixStats rowMins rowMaxs
#' @importFrom stringi stri_rand_strings
flowGraphSubset <- function(
    input_,
    meta=NULL,
    class="class",
    no_cores=1,
    markers=NULL, # including gate names
    layout_fun="layout.reingold.tilford", # layout.circle

    max_layer=NULL,

    # feature parameters
    cumsumpos=FALSE, # whether to make positives +=+/++ (cumsum)

    # save flowGraph object as folder
    path=NULL,

    # summary parameters
    summary_pars=flowGraphSubset_summary_pars(),
    summary_adjust=flowGraphSubset_summary_adjust(),

    # plotting parameters
    save_plots=TRUE
) {

    warning("The fast version of flowGraph is in beta <(@u@<)")

    ## check parameters/input ####
    if (!is.null(meta))
        if (!"id"%in%colnames(meta))
            stop("'meta' must have a column named \"id\" with sample ID's.")

    if (is.na(dim(input_)[1]))
        stop("input must be a numeric sample x cell population matrix.")

    start <- start1 <- Sys.time()
    message("preparing input; ")

    options(stringsAsFactors=FALSE)


    ## feature: count (sample x cell population)
    if (is.null(dim(input_)[1]))
        input_ <- matrix(input_, nrow=1, dimnames=list("s1",names(input_)))
    if (is.null(rownames(input_)))
        rownames(input_) <- paste0("s",seq_len(nrow(input_)))

    # remove columns with all 0's (these won't have children so it's ok)
    keepinds <- apply(input_, 2, function(x) any(x>0))
    input_ <- Matrix::Matrix(input_[,keepinds,drop=FALSE], sparse=TRUE)


    ## meta (samples)
    if (is.null(meta))
        meta <- data.frame(id=rownames(input_))
    if (nrow(meta)!= nrow(input_))
        stop("meta data and cell count matrix have different row counts")
    if (!"id"%in%colnames(meta))
        stop("meta must have an \'id\' column, id,
             whose values corresponding with sample names")
    if (!summary_pars$class%in%colnames(meta))
        stop("significance test class not found (summary_pars$class)")

    meta <- as.data.frame(meta)
    factor_col <- purrr::map_lgl(meta, is.factor)
    if (any(factor_col))
        warning("converting factors in meta into strings")
    for (j in which(factor_col))
        meta[[j]] <- as.character(meta[[j]])
    rownames(input_) <- meta$id


    ## meta (cell populations) ####
    # remove "_" between marker conditions,
    # replace "not" with "-",
    # if no "not", add "+"
    if (!any(grepl("[+-]", colnames(input_))))
        stop("make sure input type is correct")
    if (is.null(markers))
        markers <- extract_markers(colnames(input_))
    colnames(input_) <- fg_clean_phen(colnames(input_), markers)

    # meta (cell populations)
    meta_cell <- get_phen_meta(colnames(input_))
    input_ <- input_[,match(meta_cell$phenotype,colnames(input_)),drop=FALSE]

    # get rid of lower layers if requested
    if (!is.null(max_layer)) {
        phen_id <- meta_cell$phenolayer <= max_layer
        meta_cell <- meta_cell[phen_id,,drop=FALSE]
        input_ <- input_[,phen_id,drop=FALSE]
    }
    time_output(start1)
    start1 <- Sys.time()
    message("preparing to calculate features")

    ## flowGraph ####
    fg <- methods::new(
        "flowGraph",
        feat=list(node=list(count=input_), edge=list()),
        feat_desc=list(node=summary_table(input_,"count")),
        meta=meta,
        markers=markers,

        graph=list(v=meta_cell), edge_list=list(child=NULL,parent=NULL),
        plot_layout=as.character(substitute(layout_fun)),
        etc=list(cumsumpos=FALSE, class_mean_normalized=FALSE)
    )

    # feature: make + cpops' count a cumulated sum "++"="++"+"+++"
    cumsumpos <- cumsumpos & any(grepl("3",meta_cell$phenocode))
    if (cumsumpos)
        fg <- fg_feat_cumsum(fg, no_cores=no_cores)
    mc <- fg_get_feature(fg, "node", "count")

    # feature: proportion
    mp <- mc/mc[,colnames(mc)==""]
    dimnames(mp) <- dimnames(mc)


    ## START ####
    if (no_cores>1) future::plan(future::multiprocess)


    ## initialize significance test ####
    labels <- NULL
    if (is.null(summary_pars$class) & length(class)==1)
        summary_pars$class <- class
    test_custom <- test_c(summary_pars$test_custom)
    btwn_test_custom <- NULL
    if (!is.null(summary_adjust$btwn_test_custom))
        btwn_test_custom <- test_c(summary_adjust$btwn_test_custom)
    us <- unique(meta[,summary_pars$class])
    if (length(summary_pars$labels)==2 & !is.list(summary_pars$labels)) {
        labels <- summary_pars$labels
        summary_pars$labels <- lapply(summary_pars$labels, function(x)
            meta[,summary_pars$class]==x)
    }
    if (length(us)==2 &
        is.null(summary_pars$labels)) {
        summary_pars$labels <- lapply(us, function(x)
            meta[,summary_pars$class]==x)
    }
    # if it's a 2 way significance test AND labels is specified,
    # x and y are the two groups we're comparing
    # ELSE x are the values, y are the classes
    test_cust <- function(x,y, test_custom) {
        if (is.null(summary_pars$labels))
            return(test_custom(x,y))
        test_custom(x[summary_pars$labels[[1]]],x[summary_pars$labels[[2]]])
    }


    ## initialize cell population meta material for making edge list ####
    # indices for each layer
    dt <- data.table::as.data.table(meta_cell$phenolayer)[
        , list(list(.I)), by=meta_cell$phenolayer]
    lyril <- data.table::setattr(dt$V1, 'names', dt$meta_cell)
    lyrno <- length(lyril)
    lyrs <- as.numeric(names(lyril))
    lyrstf <- sapply(lyrs, function(x) (x-1)%in%lyrs & (x-2)%in%lyrs)

    # make phenocode matrix
    pc <- Matrix::Matrix(Reduce(
        rbind, lapply(stringr::str_split(meta_cell$phenocode,""),
                      as.numeric)), sparse=TRUE)
    rownames(pc) <- meta_cell$phenotype
    allcolu <- apply(pc, 2,unique)
    if (!is.list(allcolu))
        allcolu <- split(allcolu, seq(ncol(allcolu)))

    # make allcol > markeri > 0/1/2 > T/F per phenotype
    allcol <- fpurrr_map(seq_len(length(allcolu)), function(ci) {
        a <- purrr::map(allcolu[[ci]], function(ui) pc[, ci]==ui)
        names(a) <- allcolu[[ci]]; return(a)
    }, no_cores=no_cores, prll=length(allcolu)>1000)

    # split everything above by layer
    pcs <- acs <- meta_cells <- list()
    for (lyri in names(lyril)) {
        li <- lyril[[lyri]]
        pcs[[lyri]] <- pc[li,,drop=FALSE]
        acs[[lyri]] <- purrr::map(allcol, purrr::map, function(y) y[li])
        meta_cells[[lyri]] <- meta_cell[li,, drop=FALSE]
    }
    time_output(start1)
    start1 <- Sys.time()
    message("calculating features + edge lists by + p-values for cell populations with significant parents only (this is the fast version of flowGraph);")


    ## initialize features (SpecEnr, expect_prop, edge prop) for lyr0/1 ####
    # expected prop/SpecEnr for node(m), edge(e)
    me <- ep <- ms <- p <- lyrp <- sig_phens <- p1 <- NULL
    pchild <- pparen <- list()

    # root proportion is just 1, but for completeness we add it in
    if ("0"%in%names(lyril)) {
        me <- mp[,colnames(mp)=="",drop=FALSE]
        ms <- matrix(0, ncol=1, nrow=nrow(me), dimnames=dimnames(me))
        p <- 1; names(p) <- ""
        p1 <- sig_phens <- ""
        lyrp <- 0
    }
    if ("1"%in%names(lyril)) {
        p1 <- meta_cells[["1"]]$phenotype
        me <- cbind(me, matrix(.5, nrow=nrow(mp), ncol=length(p1),
                               dimnames=list(rownames(mp),p1)))
        ep <- mp[,p1,drop=FALSE]
        colnames(ep) <- paste0("_",p1)

        ## initiate child and parent edge lists for layers 0/1
        pchild <- list(p1)
        names(pchild) <- ""
        pparen <- purrr::map(seq_len(length(pchild[[1]])), function(x) "")
        names(pparen) <- pchild[[1]]

        ms_ <- ms_create(mp[,p1,drop=FALSE], me[,p1,drop=FALSE])
        ms <- cbind(ms, ms_)

        # FAST trim
        p1 <- p1[ms_psig(ms_, summary_pars, summary_adjust,
                         test_cust, test_custom, lyrno,
                         mp[,colnames(ms_),drop=FALSE],
                         me[,colnames(ms_),drop=FALSE])]
        sig_phens <- append(sig_phens, p1)
        lyrp <- 1
    }


    ## calculate features for each layer ####
    for (lyr in sort(lyrs[lyrstf])) {
        if (length(p1)==0) break

        start2 <- Sys.time()
        message("- ", length(p1), "/", length(lyril[[as.character(lyr)]]),
                " pops @ layer ", lyr)

        lyrc_ <- as.character(lyr-1)
        lyrc__ <- as.character(lyr)

        # get cell pops
        p2 <- meta_cells[[lyrc__]]$phenotype
        if (!is.null(lyrp)) #  prev layer exists?
            if (lyrp+1==lyr)
                p2 <- unique(unlist(get_child(
                    p1, pchild, pc_i=pcs[[lyrc_]], ac__=acs[[lyrc__]],
                    meta_cell__=meta_cells[[lyrc__]])))
        if (length(p2)==0) break
        p1 <- p2

        # get cell pops' parents
        pparen_ofchd <- get_paren(
            p2, pparen, pc__i=pcs[[lyrc__]], ac_=acs[[lyrc_]],
            meta_cell_=meta_cells[[lyrc_]])

        # get cell pops' grandparents
        grprnt_ofchd <- get_paren(
            unique(unlist(pparen_ofchd)), pparen,
            pc__i=pcs[[lyrc_]], ac_=acs[[as.character(lyr-2)]],
            meta_cell_=meta_cells[[as.character(lyr-2)]])

        # get min edge proportions between parent and grandparents
        edf_ <- purrr::map_dfr(names(grprnt_ofchd), function(paren)
            data.frame(from=grprnt_ofchd[[paren]], to=paren))
        ep_ <- get_eprop(
            edf_, ep,
            mp[,colnames(mp)%in%unique(append(edf_$from,edf_$to)),drop=FALSE],
            no_cores)
        ep_min <- do.call(cbind, fpurrr_map(names(grprnt_ofchd), function(x) {
            ep_col <- grepl(paste0("_",x), colnames(ep_))
            matrixStats::rowMins(ep_[,ep_col,drop=FALSE])
        }, no_cores, prll=length(grprnt_ofchd)>1000))

        # calculate expected proportion
        me_ <- do.call(cbind, fpurrr_map(p2, function(phen) {
            matrixStats::rowMins(ep_min[,pparen[[phen]],drop=FALSE]) *
                matrixStats::rowMaxs(mp[,pparen[[phen]],drop=FALSE])
        }, no_cores, prll=length(p2)>1000))
        if (is.na(dim(me_)[1]))
            me_ <- matrix(me_, ncol=length(pparen_ofchd))
        me_[is.nan(me_)] <- 0
        me_[as.matrix(me_)<0] <- 0
        colnames(me_) <- p2
        me <- cbind(me, me_)

        # calculte SpecEnr
        ms_ <- ms_create(mp[,p2,drop=FALSE], me[,p2,drop=FALSE])

        # FAST trim
        p1 <- p2[ms_psig(ms_, summary_pars, summary_adjust,
                         test_cust, test_custom, lyrno,
                         mp[,colnames(ms_),drop=FALSE],
                         me[,colnames(ms_),drop=FALSE])]
        sig_phens <- append(sig_phens, p1)
        ms <- cbind(ms, ms_)

        lyrp <- lyr
        time_output(start2)
    }
    time_output(start1, "total time used to calcuate features + edge lists")
    start1 <- Sys.time()
    if (length(sig_phens[!sig_phens%in%""])==0) {
        warning("no significant cell population phenotypes found. Try again with another set of class labels or use `flowGraph` the ful constructor instead.")
        return(NULL)
    }

    ## trim everything so there are only significant cell populations ####
    ms <- ms[,match(sig_phens, colnames(ms)),drop=FALSE]
    me <- me[,match(sig_phens, colnames(me)),drop=FALSE]
    mp <- mp[,match(sig_phens, colnames(mp)),drop=FALSE]
    rownames(ms) <- rownames(me) <- rownames(mp) <- rownames(input_)
    mc <- input_[,match(sig_phens, colnames(input_)),drop=FALSE]
    pparen <- pparen[names(pparen)%in%sig_phens]
    pparen <- purrr::map(pparen, function(x) x[x%in%sig_phens])
    pparen <- purrr::compact(pparen)
    pchild <- pchild[names(pchild)%in%sig_phens]
    pchild <- purrr::map(pchild, function(x) x[x%in%sig_phens])
    pchild <- purrr::compact(pchild)
    eps <- stringr::str_split(colnames(ep),"_")
    ep <- ep[,sapply(eps, function(x) all(x%in%sig_phens) ),drop=FALSE]
    edf <- as.data.frame(do.call(rbind, eps))
    names(edf) <- c("from","to")
    meta_cell <- meta_cell[meta_cell$phenotype%in%sig_phens,,drop=FALSE]
    markers <- extract_markers(sig_phens)

    # graph
    gr <- list(e=edf, v=meta_cell)
    gr <- set_layout_graph(gr, layout_fun) # layout cell hierarchy


    ## flowGraph ####
    desc <- do.call(rbind, list(
        summary_table(mc,"count"),
        summary_table(mp,"prop"),
        summary_table(me,"expect_prop"),
        summary_table(ms,"SpecEnr_prop")))
    desc_ <- summary_table(ep,"prop")
    fg <- methods::new(
        "flowGraph",
        feat=list(node=list(
            count=mc, prop=mp, expect_prop=me, SpecEnr_prop=ms),
            edge=list(prop=ep)),
        feat_desc=list(node=desc, edge=desc_),
        meta=meta,
        markers=markers,

        graph=gr, edge_list=list(child=pchild,parent=pparen),
        plot_layout=as.character(substitute(layout_fun)),
        etc=list(cumsumpos=fg@etc$cumsumpos, class_mean_normalized=FALSE,
                 original_count=input_)
    )


    ## calculate summary
    time_output(start1)
    start1 <- Sys.time()
    message("calculating summary statistics for SpecEnr")
    try({
        fg <- fg_summary(
            fg, no_cores=no_cores,
            class=summary_pars$class,
            label1=labels[1], label2=labels[2],
            node_features="SpecEnr",
            edge_features="NONE",
            overwrite=TRUE,
            test_custom=summary_pars$test_custom,
            test_name=summary_pars$test_name,
            diminish=summary_pars$diminish)
    })
    time_output(start1)


    ## save flowGraph object
    fg@etc$save <- list(id=stringi::stri_rand_strings(1,5))
    saved <- FALSE
    if (!is.null(path))
        tryCatch({
            fg_save(fg, path, save_plots=FALSE)
            saved <- TRUE
            fg@etc$save$path <- path
        }, error=function(e) {
            message("flowGraph object not saved; check path and try again with fg_save!")
        })

    if (save_plots & saved)
        fg_save_plots(fg, plot_path=paste0(path, "/plots"))

    time_output(start, "total time used")
    return(fg)
}
