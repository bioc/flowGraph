## GRAPH -----------------------------------------------------------------------

#' @title Determines the layer on which a phenotype resides.
#' @description Determines the layer on which the given phenotypes reside.
#' @param phen A string vector of phenotype or cell population name labels.
#' @return A numeric vector with the same length as \code{phen} indicating which
#'    layer each phenotype resides on.
#' @details Given a vector of phenotypes, returns an equal length vector of
#'    the number of markers in each phenotype.
#' @examples
#'
#'    phen <- c('A+B+C-D++', 'A+B-', '', 'B++D-E+')
#'    cell_type_layers(phen)
#'
#' @seealso
#'    \code{\link[flowGraph]{get_phen_list}}
#'    \code{\link[flowGraph]{get_phen_meta}}
#' @rdname cell_type_layers
#' @export
#' @importFrom purrr map_int
#' @importFrom stringr str_split
cell_type_layers <- function(phen)
    purrr::map_int(stringr::str_split(phen, "[+-]+"), function(x) sum(x!=""))


#' @title Genrates phenotype meta data.
#' @description Generates phenotype meta data given a vector of
#'    phenotypes and optionally phenocodes.
#' @param phen A string vector of phenotype or cell population name labels.
#' @param phenocode A string vector of phenocodes corresponding
#'    to the phenotypes in \code{phen}.
#' @return A data frame with columns containing meta data on
#'    cell poulation nodes with columns:
#'    \itemize{
#'        \item{\code{phenotype}: cell population node label e.g. "A+B+".}
#'        \item{\code{phenocode}: a string penocode containing a numeric
#'         corresponding to the phenotype column e.g. "2200".}
#'        \item{\code{phenolayer}: a numeric layer on which a cell population
#'         resides in e.g. 2.}
#'    }
#' @examples
#'
#'    phen <- c('A+B+C-D++', 'A+B-', '', 'B++D-E+')
#'    phenc <- c('22130','21000','00000','03012')
#'    get_phen_meta(phen, phenc)
#'
#' @seealso
#'    \code{\link[flowGraph]{get_phen_list}}
#'    \code{\link[flowGraph]{cell_type_layers}}
#' @rdname get_phen_meta
#' @export
#' @importFrom stringr str_split str_extract_all str_count
#' @importFrom purrr map_chr map_int map
get_phen_meta <- function(phen, phenocode=NULL) {
    pm <- data.frame(phenotype=phen)
    markers <- unlist(stringr::str_split(phen, "[-+]"))
    markers <- gsub("_","",markers) # underscore between marker conditions
    markers <- unique(markers[markers != ""])
    if (is.null(phenocode)) {
        phenocode <- purrr::map_chr(phen, function(x) {
            if (x == "") return(paste0(rep(0, length(markers)), collapse=""))
            b <- stringr::str_split(x, "[+-]+[_]*")[[1]]
            b <- b[-length(b)]
            bo <- match(markers, b)
            phec <- purrr::map_int(
                unlist(stringr::str_extract_all(x, "[+-]+")), function(y)
                    stringr::str_count(y, "[+]")) + 1
            c <- phec[bo]
            c[is.na(c)] <- 0
            paste0(c, collapse="")
        })
    }
    pm$phenocode <- phenocode
    pm$phenolayer <- cell_type_layers(phen)
    pm$phenotype <- as.character(pm$phenotype)
    pm$phenogroup <- gsub("[+-]+", "_", pm$phenotype)

    return(pm)
}


#' @title Creates edge lists.
#' @description Creates edge lists indicating relationships
#' between cell populations
#'    given meta data on these cell populations
#'    produced by the \code{get_phen_meta} function.
#' @param meta_cell A data frame containing meta data on cell populations
#'    as produced by the \code{get_phen_meta} function.
#' @param phen A string vector of phenotype or cell population name labels.
#'    Cannot be set to \code{NULL} if \code{meta_cell} is set to \code{NULL}.
#' @param no_cores An integer indicating how many cores to parallelize on.
#' @return A list containing 'pchild', an edge list
#'    indicating where edges point to,
#'    'pparen', an edge list indicating where edges point from, and
#'    'edf', a data frame where each row contains
#'    the nodes an edge points 'from' and 'to'.
#' @examples
#'
#'    phen <- c('A+B-C+', 'A+B-', 'A+')
#'    get_phen_list(phen=phen)
#'
#' @seealso
#'    \code{\link[flowGraph]{get_phen_meta}}
#'    \code{\link[flowGraph]{cell_type_layers}}
#' @rdname get_phen_list
#' @export
#' @importFrom future plan multisession sequential
#' @importFrom stringr str_split str_extract_all
#' @importFrom Matrix Matrix
#' @importFrom purrr map compact map_chr
#' @importFrom furrr future_map
#' @importFrom data.table as.data.table setattr
get_phen_list <- function(meta_cell=NULL, phen=NULL, no_cores=1) {
    if (no_cores>1) future::plan(future::multisession)

    if (is.null(phen) & is.null(meta_cell))
        stop("give me something!")
    if (is.null(meta_cell))
        meta_cell <- get_phen_meta(phen)


    ## initialize cell population meta material for making edge list ####
    # indices for each layer
    dt <- data.table::as.data.table(meta_cell$phenolayer)[
        , list(list(.I)), by=meta_cell$phenolayer]
    lyril <- data.table::setattr(dt$V1, 'names', dt$meta_cell)
    lyrs <- length(lyril)

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
    }, no_cores=no_cores)

    # split everything above by layer
    pcs <- acs <- meta_cells <- list()
    for (lyri in names(lyril)) {
        li <- lyril[[lyri]]
        pcs[[lyri]] <- pc[li,,drop=FALSE]
        acs[[lyri]] <- purrr::map(allcol, purrr::map, function(y) y[li])
        meta_cells[[lyri]] <- meta_cell[li,, drop=FALSE]
    }

    ## initialize
    pchild <- pparen <- list()

    # root proportion is just 1, but for completeness we add it in
    if (all(c("0", "1")%in%names(lyril))) {
        ## initiate child and parent edge lists for layers 0/1
        pchild <- list(meta_cells[["1"]]$phenotype)
        names(pchild) <- ""
        pparen <- purrr::map(seq_len(length(pchild[[1]])), function(x) "")
        names(pparen) <- pchild[[1]]
    }


    ## calculate features for each layer ####
    lyrs <- sort(unique(meta_cell$phenolayer))
    lyrstf <- sapply(lyrs, function(x) (x-1)%in%lyrs)

    for (lyr in lyrs[lyrstf]) {
        start2 <- Sys.time()
        message("- ", length(lyril[[as.character(lyr)]]),
                        " pops @ layer ", lyr)

        lyrc_ <- as.character(lyr-1)
        lyrc__ <- as.character(lyr)
        meta_cell_grid_ <- pcs[[lyrc_]]
        meta_cell_grid__ <- pcs[[lyrc__]]
        allcol_ <- acs[[lyrc_]]
        allcol__ <- acs[[lyrc__]]
        meta_cell_ <- meta_cells[[lyrc_]]
        meta_cell__ <- meta_cells[[lyrc__]]

        # child
        pchildl <- fpurrr_map(seq_len(nrow(meta_cell_)), function(j) {
            mcgrow <- meta_cell_grid_[j, ]
            colj1 <- which(mcgrow > 0)
            chi <- Reduce("&", purrr::map(colj1, function(coli)
                allcol__[[coli]][[as.character(mcgrow[coli])]]))
            meta_cell__$phenotype[chi]
        }, no_cores=no_cores, prll=nrow(meta_cell__) > 5000)
        names(pchildl) <- meta_cell_$phenotype
        pchildl <- purrr::compact(pchildl)
        pchild <- append(pchild, pchildl)

        # paren
        pparenl <- fpurrr_map(seq_len(nrow(meta_cell__)), function(j) {
            mcgrow <- meta_cell_grid__[j, ]
            colj1 <- which(mcgrow > 0)
            chidf <- do.call(cbind, purrr::map(colj1, function(coli)
                allcol_[[coli]][[as.character(mcgrow[coli])]] ))
            chi <- apply(chidf, 1, function(x) sum(!x) == 1)
            meta_cell_$phenotype[chi]
        }, no_cores=no_cores, prll=nrow(meta_cell__) > 5000)
        names(pparenl) <- meta_cell__$phenotype
        pchildl <- purrr::compact(pparenl)
        pparen <- append(pparen, pparenl)

        time_output(start2)
    }
    edf <- do.call(rbind,purrr::map(seq_len(length(pchild)), function(x)
        data.frame(from=names(pchild)[x], to=pchild[[x]])) )

    temp_se <- function(x) stringr::str_extract_all(x, "[^_^+^-]+[+-]+")
    from_ <- temp_se(edf$from)
    to_ <- temp_se(edf$to)
    edf$marker <- unlist(fpurrr_map(seq_len(length(from_)), function(x)
        setdiff(to_[[x]], from_[[x]]), no_cores=no_cores, prll=nrow(edf)>5000))

    future::plan(future::sequential)
    return(list(pchild=pchild, pparen=pparen, edf=edf))
}



#' @title Determines cell hierarchy layout.
#' @description Determines cell hierarchy layout and returns the X, Y coordinate
#'    of each cell population.
#' @param gr A list containing data frames \code{e} and \code{v}.
#' @param layout_fun A string of a function from the \code{igraph} package that
#'    indicates what layout should be used if a cell hierarchy is to be ploted;
#'    all such functions have prefix \code{layout_}
#'    e.g. \code{layout_fun="layout.reingold.tilford"}.
#' @return A list containing data frames \code{e} and \code{v}; each data frame
#'    contains an X, Y column or coordinate for each node and edge.
#' @examples
#'
#'    no_cores <- 1
#'    data(fg_data_pos30)
#'    fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                                    prop=FALSE, specenr=FALSE,
#'                                    no_cores=no_cores)
#'
#'    head(set_layout_graph(fg_get_graph(fg)))
#'
#' @rdname set_layout_graph
#' @export
#' @import igraph
#' @importFrom igraph graph_from_data_frame layout.reingold.tilford V
#' @importFrom purrr map_int
set_layout_graph <- function(gr, layout_fun="layout.reingold.tilford") {
    # layout.circle, layout.reingold.tilford FUN is a layout
    # function from the igraph package assume
    # graph is connected, used internally

    gr_e <- gr$e
    gr_v <- gr$v

    # edit layout
    gr0 <- igraph::graph_from_data_frame(gr_e[,c("from","to")])
    if (layout_fun!="layout.reingold.tilford") {
        layout_fun <- get(layout_fun)
        gr_vxy_ <- layout_fun(gr0)
        gr_vxy <- as.data.frame(gr_vxy_)
    } else {
        gr_vxy_ <- igraph::layout.reingold.tilford(gr0)
        gr_vxy <- as.data.frame(gr_vxy_)

        ## edit layout manually
        gr_vxy <- space_hierarchy(gr_vxy)
        # turn plot sideways
        gr_vxy <- gr_vxy[, 2:1]
        xmax <- max(gr_vxy[, 1])
        if (gr_vxy[1, 1] == xmax)
            gr_vxy[, 1] <- xmax - gr_vxy[, 1]
    }

    # get node
    colnames(gr_vxy) <- c("x", "y")
    gr_v_xy <- gr_vxy[match(gr_v$phenotype, names(igraph::V(gr0)[[]])), ]
    gr_v$x <- gr_v_xy$x
    gr_v$y <- gr_v_xy$y

    # get edge
    e_match <- match(gr_e$from, gr_v$phenotype)
    gr_e$from.x <- gr_v$x[e_match]
    gr_e$from.y <- gr_v$y[e_match]
    e_match2 <- match(gr_e$to, gr_v$phenotype)
    gr_e$to.x <- gr_v$x[e_match2]
    gr_e$to.y <- gr_v$y[e_match2]

    return(list(v=gr_v, e=gr_e))
}

space_hierarchy <- function(gr_vxy) {
    gys <- sort(unique(gr_vxy[, 2]))
    gxns <- purrr::map_int(gys, function(y) sum(gr_vxy[, 2] == y))
    maxlayern <- max(gxns)
    maxlayer <- gys[which.max(gxns)]
    maxlayertf <- gr_vxy[, 2] == maxlayer
    gr_vxy[maxlayertf, 1] <- rank(gr_vxy[maxlayertf, 1]) - 1
    for (gy in gys) {
        if (gy == maxlayer)
            next()
        layertf <- gr_vxy[, 2] == gy
        gr_vxy[layertf, 1]=rank(gr_vxy[layertf, 1]) - 1 +
            floor((maxlayern - sum(layertf))/2)
    }
    return(gr_vxy)
}


#' @title Determines cell hierarchy layout.
#' @description Determines cell hierarchy layout and returns the X, Y coordinate
#'    of each cell population. This function is a wrapper for
#'    \code{\link[flowGraph]{set_layout_graph}}.
#' @param fg flowGraph object.
#' @param layout_fun A string version of a function name from
#'    the \code{igraph} package that
#'    indicates what layout should be used if a cell hierarchy is to be ploted;
#'    all such functions have prefix \code{layout_}
#'    e.g. \code{layout_fun="layout.reingold.tilford"}.
#' @return flowGraph object with coordinate meta data on cell populations
#'    and edges for plotting use.
#' @details Given a flowGraph object, modifies the \code{graph} slot such that
#'    it contains X, Y axes for each node in accordance to a user specified
#'    layout.
#' @examples
#'
#'    no_cores <- 1
#'    data(fg_data_pos30)
#'    fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                                    prop=FALSE, specenr=FALSE,
#'                                    no_cores=no_cores)
#'
#'    fg <- fg_set_layout(fg)
#'    head(fg_get_graph(fg)$v)
#'
#' @rdname fg_set_layout
#' @export
fg_set_layout <- function(fg, layout_fun="layout.reingold.tilford") {
    fg@plot_layout <- layout_fun
    fg@graph <- set_layout_graph(fg_get_graph(fg), layout_fun)
    return(fg)
}
