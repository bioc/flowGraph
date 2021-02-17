#' @title Generates the proportion node feature.
#' @description Generates the proportion node feature and returns it
#'  inside the returned flowGraph object.
#' @param fg flowGraph object.
#' @param overwrite A logical variable indicating whether to
#'  overwrite the existing proportion node feature if it exists.
#' @return flowGraph object containing the proportion node feature.
#' @details Given a flowGraph object, \code{fg_feat_node_prop} returns the same
#'  flowGraph object, inside of which is an additional proportions \code{prop}
#'  \code{node} feature
#'  and its meta data. The proportions feature is made using the node count
#'  feature and is the cell count of each cell population over the total
#'  cell count.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=no_cores)
#'
#'  fg <- fg_feat_node_prop(fg)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_feat_node_specenr}}
#'  \code{\link[flowGraph]{fg_add_feature}}
#'  \code{\link[flowGraph]{fg_get_feature}}
#'  \code{\link[flowGraph]{fg_rm_feature}}
#'  \code{\link[flowGraph]{fg_get_feature_desc}}
#' @rdname fg_feat_node_prop
#' @export
fg_feat_node_prop <- function(fg, overwrite=FALSE) {
    start1 <- Sys.time()
    message("preparing feature(s): node proportion ")

    fg <- fg_add_feature(fg, type="node", feature="prop", overwrite=overwrite,
                         m=NULL, feat_fun=fg_feat_node_prop_)

    time_output(start1)
    return(fg)
}

fg_feat_node_prop_ <- function(fg) {
    mc <- fg_get_feature(fg, "node", "count")
    mp <- mc/mc[,colnames(mc)==""]
    dimnames(mp) <- dimnames(mc)
    return(mp)
}


#' @title Generates the proportion edge feature.
#' @description Generates the proportion edge feature and returns it
#'  inside the flowGraph object.
#' @param fg flowGraph object.
#' @param no_cores An integer indicating how many cores to parallelize on.
#' @param overwrite A logical variable indicating whether to
#'  overwrite the existing proportion edge feature if it exists.
#' @return flowGraph object containing the proportion edge feature.
#' @details Given a flowGraph object, \code{fg_feat_edge_prop} returns the same
#'  flowGraph object with an additional proportions \code{prop} \code{edge}
#'  feature and its meta data. The proportions feature is
#'  made using the node count
#'  feature and is the cell count of each cell population (e.g. A+B+)
#'  over the cell count of its parent (e.g. A+);
#'  each edge then corresponds with such a relationship.
#'  The edge feature matrix has column names <from>_<to> e.g. A+_A+B+.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=no_cores)
#'
#'  fg <- fg_feat_edge_prop(fg)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_feat_node_prop}}
#'  \code{\link[flowGraph]{fg_feat_node_specenr}}
#'  \code{\link[flowGraph]{fg_add_feature}}
#'  \code{\link[flowGraph]{fg_get_feature}}
#'  \code{\link[flowGraph]{fg_rm_feature}}
#'  \code{\link[flowGraph]{fg_get_feature_desc}}
#' @rdname fg_feat_prop_edge
#' @export
fg_feat_edge_prop <- function(fg, no_cores=1, overwrite=FALSE) {
    start1 <- Sys.time()
    message("preparing feature(s): proportion on edges ")

    fg <- fg_add_feature(fg, type="edge", feature="prop", overwrite=overwrite,
                         m=NULL, feat_fun=fg_feat_edge_prop_, no_cores=no_cores)

    time_output(start1)
    return(fg)
}

#' @importFrom future plan multiprocess
#' @importFrom furrr future_map
#' @importFrom purrr map
fg_feat_edge_prop_ <- function(fg, no_cores=1) {
    if (no_cores>1) future::plan(future::multiprocess)

    edf <- fg_get_graph(fg)$e
    mc <- fg_get_feature(fg, "node", "count")

    rooti <- which(colnames(mc)=="")

    childprop_ <- do.call(cbind, fpurrr_map(seq_len(nrow(edf)), function(i) {
        pname <- ifelse(edf$from[i]=="", rooti, edf$from[i])
        mc[,edf$to[i]]/mc[,pname]
    }, no_cores=no_cores, prll=nrow(edf)>1000))
    colnames(childprop_) <- paste0(edf$from,"__",edf$to)
    childprop_[is.nan(as.matrix(childprop_))] <- 0

    return(childprop_)
}


#' @title Generates the SpecEnr edge feature.
#' @description Generates the SpecEnr edge feature and returns it
#'  inside the flowGraph object.
#' @param fg flowGraph object.
#' @param no_cores An integer indicating how many cores to parallelize on.
#' @param overwrite A logical variable indicating whether to
#'  overwrite the existing proportion edge feature if it exists.
#' @return flowGraph object containing the proportion edge feature.
#' @details Given a flowGraph object, \code{fg_feat_edge_SpecEnr} returns the same
#'  flowGraph object with an additional SpecEnr and expected proportions
#'  \code{expect_prop} \code{edge}
#'  feature and its meta data. The expected proportions edge feature is
#'  calculated by taking the ratio of the child nodes' (e.g. A+B+)
#'  expected proportion value over
#'  its parent nodes' (e.g. A+) actual proportion value.
#'  The SpecEnr feature is the actual over expected proportion ratio, logged.
#'  The edge feature matrix has column names <from>_<to> e.g. A+_A+B+.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=no_cores)
#'
#'  fg <- fg_feat_edge_specenr(fg)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_feat_node_prop}}
#'  \code{\link[flowGraph]{fg_feat_node_specenr}}
#'  \code{\link[flowGraph]{fg_add_feature}}
#'  \code{\link[flowGraph]{fg_get_feature}}
#'  \code{\link[flowGraph]{fg_rm_feature}}
#'  \code{\link[flowGraph]{fg_get_feature_desc}}
#' @rdname fg_feat_edge_specenr
#' @export
fg_feat_edge_specenr <- function(fg, no_cores=1, overwrite=FALSE) {
    start1 <- Sys.time()
    message("preparing feature(s): expected proportion on edges ")

    # mp: sample x cell population, cell count matrix
    if (is.null(fg_get_feature_all(fg)$node$prop) | overwrite)
        fg <- fg_feat_node_prop(fg)

    # ep: sample x edge, proportion matrix i.e. if column name is A+B+C+_A+B+,
    # then the value is the count of A+B+C+ over A+B+.
    if (is.null(fg_get_feature_all(fg)$edge$prop) | overwrite)
        fg <- fg_feat_edge_prop(fg, no_cores=no_cores)

    # create expected features
    if (is.null(fg_get_feature_all(fg)$node$expect_prop) | overwrite)
        fg <- fg_add_feature(fg, type="node", feature="expect_prop",
                             overwrite=overwrite,
                             m=NULL, feat_fun=fg_feat_node_exprop_,
                             no_cores=no_cores)

    if (is.null(fg_get_feature_all(fg)$edge$expect_prop) | overwrite)
    fg <- fg_add_feature(fg, type="edge", feature="expect_prop",
                         overwrite=overwrite,
                         m=NULL, feat_fun=fg_feat_edge_exprop_, no_cores=no_cores)

    mp <- fg_get_feature(fg, "edge", "prop")
    mep <- fg_get_feature(fg, "edge", "expect_prop")

    a <- mp/mep
    aa <- as.matrix(a)
    a[is.infinite(aa)] <- max(a[is.finite(aa)])
    a[is.nan(aa)] <- 0
    e0 <- as.matrix(mep)==0
    suppressWarnings({ specenr1 <- log(a) })
    specenr1[is.nan(as.matrix(specenr1))] <- 0
    specenr1[e0] <- log(as.matrix(mp)[e0])
    specenr1[mp==0] <- 0

    fg <- fg_add_feature(
        fg, type="edge", feature="SpecEnr",
        m=specenr1, no_cores=no_cores, overwrite=overwrite)

    time_output(start1)
    return(fg)
}

#' @importFrom future plan multiprocess
#' @importFrom furrr future_map
#' @importFrom purrr map
fg_feat_edge_exprop_ <- function(fg, no_cores=1) {
    if (no_cores>1) future::plan(future::multiprocess)

    edf <- fg_get_graph(fg)$e
    mp <- fg_get_feature(fg, "node", "prop")
    mep <- fg_get_feature(fg, "node", "expect_prop")

    rooti <- which(colnames(mp)=="")

    childprop_ <- do.call(cbind, fpurrr_map(seq_len(nrow(edf)), function(i) {
        pname <- ifelse(edf$from[i]=="", rooti, edf$from[i])
        mep[,edf$to[i]]/mp[,pname]
    }, no_cores=no_cores, prll=nrow(edf)>1000))
    colnames(childprop_) <- paste0(edf$from,"__",edf$to)
    childprop_[is.nan(as.matrix(childprop_))] <- 0

    return(childprop_)
}


#' @title Generates the SpecEnr node feature.
#' @description Generates the SpecEnr node feature and returns it
#'  inside the returned flowGraph object.
#' @param fg flowGraph object
#' @param no_cores An integer indicating how many cores to parallelize on.
#' @param feature A string indicating feature name; this is the feature
#'  SpecEnr will be calculated on.
#' @param overwrite A logical variable indicating whether to
#'  overwrite the existing SpecEnr node feature if it exists.
#' @return flowGraph object containing the SpecEnr node feature.
#' @details Given a flowGraph object, \code{fg_feat_node_specenr}
#'  returns the same
#'  flowGraph object with an additional \code{SpecEnr} and \code{expect_prop}
#'  \code{node} feature and its meta data.
#'  The expected proportions feature is made using the \code{prop} \code{node}
#'  and \code{edge} features; therefore, the returned flowGraph will also
#'  contain these two features. For details on how these feature is calculated.
#'
#' @references
#'  \insertRef{yue2019identifying}{flowGraph}
#'
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=no_cores)
#'
#'  # SpecEnr is by default calculated based on proportions
#'  fg <- fg_feat_node_specenr(fg, no_cores=no_cores)
#'
#'  # SpecEnr can be calculated for other feature values too
#'  fg <- fg_feat_node_specenr(fg, feature="count")
#'
#'  show(fg)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_feat_node_prop}}
#'  \code{\link[flowGraph]{fg_add_feature}}
#'  \code{\link[flowGraph]{fg_get_feature}}
#'  \code{\link[flowGraph]{fg_rm_feature}}
#'  \code{\link[flowGraph]{fg_get_feature_desc}}
#' @rdname fg_feat_node_specenr
#' @export
#' @importFrom Rdpack reprompt
fg_feat_node_specenr <- function(fg,no_cores=1,feature="prop",overwrite=FALSE) {

    start1 <- Sys.time()
    message("preparing feature(s): SpecEnr ")

    fg_feat <- fg_get_feature_all(fg)

    # mp: sample x cell population, cell count matrix
    if (is.null(fg_feat$node$prop) | overwrite)
        fg <- fg_feat_node_prop(fg)

    # ep: sample x edge, proportion matrix i.e. if column name is A+B+C+_A+B+,
    # then the value is the count of A+B+C+ over A+B+.
    if (is.null(fg_feat$edge$prop) | overwrite)
        fg <- fg_feat_edge_prop(fg, no_cores=no_cores)

    # create expeced features
    if (is.null(fg_feat$node$expect_prop) | overwrite)
        fg <- fg_add_feature(fg, type="node", feature="expect_prop",
                             overwrite=overwrite,
                             m=NULL, feat_fun=fg_feat_node_exprop_,
                             no_cores=no_cores)

    if (is.null(fg_feat$node[[paste0("expect_",feature)]]) | overwrite)
        fg <- fg_add_feature(
            fg, type="node", feature=paste0("expect_",feature),
            m=fg_get_feature(fg, type="node", feature=feature)[,1] *
                fg_get_feature(fg, type="node", feature="expect_prop"))

    mp <- fg_get_feature(fg, "node", feature)
    exp1 <- fg_get_feature(fg, "node", paste0("expect_",feature))

    a <- mp/exp1 ### SPECENR ---------------------------
    aa <- as.matrix(a)
    a[is.infinite(aa)] <- max(a[is.finite(aa)])
    a[is.nan(aa)] <- 0
    e0 <- as.matrix(exp1)==0
    suppressWarnings({ specenr1 <- log(a) })
    specenr1[is.nan(as.matrix(specenr1))] <- 0
    specenr1[e0] <- log(as.matrix(mp)[e0])
    specenr1[mp==0] <- 0

    fg <- fg_add_feature(
        fg, type="node", feature=paste0( "SpecEnr", ifelse(
            feature=="prop", "", paste0("_",feature)) ),
        m=specenr1, no_cores=no_cores, overwrite=overwrite)

    time_output(start1)
    return(fg)
}

# expected proportion: this is the version currently used.
# this version is the same as the new version, but the new
# version is easier to experiment on and calculates everything at once.
#' @importFrom future plan multiprocess
#' @importFrom stringr str_extract_all str_extract str_count
#' @importFrom furrr future_map
#' @importFrom purrr map
#' @importFrom stats median
fg_feat_node_exprop_ <- function(fg, no_cores=1) {
    # prepare parallel backend
    if (no_cores>1) future::plan(future::multiprocess)

    fg_graph <- fg_get_graph(fg)

    # meta data for samples
    meta_cell <- fg_graph$v

    # pparen: edge list i.e. the name of elements in the list are cell pops;
    # the vector in each element lists corresponding parent cell populations.
    # pchild: edge list i.e. the name of elements in the list are cell pops;
    # the vector in each element lists corresponding child cell populations.
    pparen <- fg_get_el(fg)$parent
    # child: pchild <- fg_get_el(fg)$child

    # mp: sample x cell population, cell count matrix
    mp <- fg_get_feature(fg, "node", "prop")

    # ep: sample x edge, proportion matrix i.e. if column name is A+B+C+_A+B+,
    # then the value is the count of A+B+C+ over A+B+.
    ep <- fg_get_feature(fg, "edge", "prop")
    phens <- fg_graph$v$phenotype[
        fg_graph$v$phenolayer != max(fg_graph$v$phenolayer)]
    phens <- phens[phens!=""]
    # max edges

    ep_ <- do.call(cbind, fpurrr_map(phens, function(phen)
        matrixStats::rowMins(ep[,paste0(pparen[[phen]], "__", phen),drop=FALSE])
        , no_cores=no_cores, prll=length(phens)>1000))
    colnames(ep_) <- phens

    rooti <- which(colnames(mp)=="")

    ## start calculating expected proportion
    ## first we prepare the list of cell populations;
    ## note we only calc expected proportion for cell populations in layers 2+.
    # vector of cell populations in the zeroth and first layer.
    cells1 <- append("",meta_cell$phenotype[meta_cell$phenolayer==1])
    expe1 <- matrix(.5,nrow=nrow(mp), ncol=length(cells1),
                          dimnames=list(rownames(mp),cells1))
    expe1[,1] <- 1

    # vector of cell populations in layers 2+.
    cells <- meta_cell$phenotype[meta_cell$phenolayer>1]
    # vector of the layer in which the cell populations in layers 2+ reside on.
    cellsn <- meta_cell$phenolayer[meta_cell$phenolayer>1]
    # logical: if there is any e.g. A++, in this data set
    # (as opposed to just A+, A-)
    multiplus <- any(nchar(unlist(
        stringr::str_extract_all(cells,"[+]+")))>1)
    cells1_p <- stringr::str_extract(cells1,"[+]+")
    maxp <- max(nchar(cells1_p[!is.na(cells1_p)]))

    ## calc expected props for cell pops w/ positive marker conds only ---------
    # this reduces computation time e.g. A+B+C+, not A+B+C-
    # cpind <- seq_len(length(cells))
    cpind <- which(!grepl("[-]",cells))
    # all: loop_ind <- loop_ind_f(1:length(cells),no_cores)
    expecp <- do.call(cbind, fpurrr_map(cpind, function(ic) {
        # cell population name e.g. A+B+C+
        cpop <- cells[ic]

        # node proportion values extracted from mp for cpop's parents
        # e.g. A+B+, A+C+, B+C+; and grandparents e.g. A+, B+, C+.
        pnames <- pparen[[cpop]]
        parent <- mp[,pnames,drop=FALSE] # proportion values of parents.

        # edge proportion values extracted from ep for edges between
        # cpop's parents and grandparents.
        pedges <- ep_[,pnames,drop=FALSE]

        # using the above, get expected proportion; see formula in
        # https://www.biorxiv.org/content/10.1101/837765v2
        expect1 <- apply(pedges,1,min) * apply(parent,1,max)
        return(expect1)
    }, no_cores=no_cores, prll=length(cbind)>1000))
    expecp[is.nan(expecp)] <- 0
    colnames(expecp) <- cells[cpind]


    ## infer the expected prop of all other cell pops --------------------------
    # this isn't necessary, we could've just done it all above,
    # but this saves time; see same paper as above link
    expec0 <- expec <- as.matrix(cbind(expe1,expecp))
    cpopneg <- setdiff(cells,colnames(expec))
    cpopnegl <- cell_type_layers(cpopneg)

    # placeholder;
    # "count" is a variable given by gsubfn to any function it is given
    # count <- 0
    csp <- fg_get_etc(fg)$cumsumpos | !multiplus

    for (lev in sort(unique(cpopnegl))) {
        sibsl <- cells[cellsn==lev]
        cpopl <- cpopneg[cpopnegl==lev]
        # number of negative marker conditions
        cpopnegno <- stringr::str_count(cpopl,"[-]")
        cpopnegnos <- purrr::map(sort(unique(cpopnegno)),
                                 function(x) cpopl[cpopnegno==x])
        for (cpops in cpopnegnos) {
            if (csp) {
                expecn <- do.call(cbind, fpurrr_map(cpops, function(cpop) {
                    sib <- gsub('(.*?)-(.*)', '\\1+\\2', cpop)
                    if (!sib%in%colnames(expec)) {
                        pari <- which(purrr::map_lgl(
                            pparen[[cpop]],
                            ~grepl(.x,sib)))
                        if (length(pari)==0) return(rep(0,nrow(mp)))
                        return(mp[,pparen[[cpop]][pari[1]]])
                    }
                    pname <- intersect(pparen[[cpop]],
                                       pparen[[sib]])
                    return(mp[,pname] - expec[,sib])
                }, no_cores=no_cores, prll=length(cpops)>1000))
            } else {
                expecn <- do.call(cbind, fpurrr_map(cpops, function(cpop) {
                    cpopgsub <- gsub("[+]","[+]", cpop)
                    cpopgsub <- gsub('(.*?)-(.*)', '\\1[+]+\\2',
                                     cpopgsub)
                    cpopgsub <- gsub("[-]","[-]", cpopgsub)
                    sibs <- sibsl[grepl(cpopgsub,sibsl)]
                    sibs_ <- sibs%in%colnames(expec)
                    if (sum(sibs_)==0) {
                        sibs__ <- paste0(sibs, collapse="")
                        pari <- which(purrr::map_lgl(
                            pparen[[cpop]],
                            ~grepl(.x,sibs__)))
                        if (length(pari)==0)
                            return(rep(0,nrow(mp)))
                        return(mp[,pparen[[cpop]][pari[1]]])
                    }
                    sibs <- sibs[sibs_]
                    pname <- intersect(pparen[[cpop]],
                                       pparen[[sibs[1]]])
                    return(mp[,pname] - rowSums(expec[,sibs,drop=FALSE]))
                }, no_cores=no_cores, prll=length(cpops)>1000))
            }
            colnames(expecn) <- cpops
            expec <- cbind(expec, expecn)
        }
    }

    exp1 <- cbind(expe1,expec[,match(cells,colnames(expec)),
                                    drop=FALSE])
    # some: exp1 <- cbind(expe1,expecp[,match(cells,colnames(expecp)),drop=FALSE])

    exp1[as.matrix(exp1)<0] <- 0
    # some: exp1 <- cbind(expe1,expecp[,match(cells,colnames(expecp)),drop=FALSE])
    return(exp1) ### EXPECTED PROPORTION
    # mp/exp1 ### SPECENR
}


fg_feat_node_exprop_new <- function(fg, no_cores=1) {
    # prepare parallel backend
    if (no_cores>1) future::plan(future::multiprocess)

    fg_graph <- fg_get_graph(fg)

    # meta data for samples
    meta_cell <- fg_graph$v

    # pparen: edge list i.e. the name of elements in the list are cell pops;
    # the vector in each element lists corresponding parent cell populations.
    # pchild: edge list i.e. the name of elements in the list are cell pops;
    # the vector in each element lists corresponding child cell populations.
    pparen <- fg_get_el(fg)$parent
    # child: pchild <- fg@edge_list$child

    # mp: sample x cell population, cell count matrix
    mp <- fg_get_feature(fg, "node", "prop")

    # ep: sample x edge, proportion matrix i.e. if column name is A+B+C+_A+B+,
    # then the value is the count of A+B+C+ over A+B+.
    ep <- fg_get_feature(fg, "edge", "prop")
    phens <- fg_graph$v$phenotype[
        fg_graph$v$phenolayer != max(fg_graph$v$phenolayer)]
    phens <- phens[phens!=""]

    rooti <- which(colnames(mp)=="")

    ep_ <- do.call(cbind, fpurrr_map(phens, function(phen)
        matrixStats::rowMins(ep[,paste0(pparen[[phen]], "__", phen),drop=FALSE])
        , no_cores=no_cores, prll=length(phens)>1000))
    colnames(ep_) <- phens


    ## start calculating expected proportion
    ## first we prepare the list of cell populations;
    ## note we only calc expected proportion for cell populations in layers 2+.
    # vector of cell populations in the zeroth and first layer.
    cells1 <- append("",meta_cell$phenotype[meta_cell$phenolayer==1])
    expe1 <- matrix(.5,nrow=nrow(mp), ncol=length(cells1),
                          dimnames=list(rownames(mp),cells1))
    expe1[,1] <- 1

    # vector of cell populations in layers 2+.
    cells <- meta_cell$phenotype[meta_cell$phenolayer>1]
    # vector of the layer in which the cell populations in layers 2+ reside on.
    cellsn <- meta_cell$phenolayer[meta_cell$phenolayer>1]
    # logical: if there is any e.g. A++, in this data set
    # (as opposed to just A+, A-)
    multiplus <- any(nchar(unlist(
        stringr::str_extract_all(cells,"[+]+")))>1)
    cells1_p <- stringr::str_extract(cells1,"[+]+")
    maxp <- max(nchar(cells1_p[!is.na(cells1_p)]))

    ## calc expected props for cell pops w/ positive marker conds only ---------
    # this reduces computation time e.g. A+B+C+, not A+B+C-
    cpind <- seq_len(length(cells))
    # all: cpind <- which(!grepl("[-]",cells))

    # get all pairs of marker condition indices
    markers <- fg_get_markers(fg)
    ml <- max(fg_graph$v$phenolayer)
    mlcomb <- purrr::map(2:ml, function(mly)
        Reduce(cbind,purrr::map(seq_len(mly-1), function(x) {
        xl <- (x+1):mly
        rbind(xl,rep(x,mly-x))
    })))


    cells_ <- stringr::str_extract_all(cells, "[^_^+^-]+[+-]+")
    # some: loop_ind <- loop_ind_f(1:length(cells),no_cores)
    expecp <- do.call(cbind, fpurrr_map(cpind, function(ic) {
        # cell population name e.g. A+B+C+
        cpop <- cells[ic]
        cmarkers <- cells_[[ic]]

        # node proportion values extracted from mp for cpop's parents
        # e.g. A+B+, A+C+, B+C+; and grandparents e.g. A+, B+, C+.
        pnames <- pparen[[cpop]]
        parent <- mp[,pnames,drop=FALSE] # proportion values of parents.
        if (length(cmarkers)==2) return(apply(parent,1,prod))


        # try all combinations, because...
        mpi <- mp[,cpop]
        mlcombi <- mlcomb[[length(cmarkers)-1]]
        allcombexp <- purrr::map(seq_len(ncol(mlcombi)),function(yi) {
            y <- mlcombi[,yi]
            par1 <- gsub(cmarkers[y[1]],"",cpop,fixed=TRUE)
            mpe1p <- mp[,gsub(cmarkers[y[1]],"",cpop,fixed=TRUE)]
            mpe1e <- ep[,paste0(gsub(cmarkers[y[2]],"",par1,fixed=TRUE),
                                "_",gsub(cmarkers[y[2]],"",cpop,fixed=TRUE))]
            mpe1 <- mpe1p*mpe1e

            par1 <- gsub(cmarkers[y[2]],"",cpop,fixed=TRUE)
            mpe2p <- mp[,gsub(cmarkers[y[2]],"",cpop,fixed=TRUE)]
            mpe2e <- ep[,paste0(gsub(cmarkers[y[1]],"",par1,fixed=TRUE),
                                "_",gsub(cmarkers[y[1]],"",cpop,fixed=TRUE))]
            mpe2 <- mpe2p*mpe2e # mpe1 and mpe are the same
            return(cbind(mpe1p,mpe1e,mpe2p,mpe2e,mpe1))
        })
        acediff <- sapply(allcombexp, function(x)
            stats::median(abs(x[,5]-mp[,cpop])))
        return(allcombexp[[which.min(acediff)]][,5])

        allcombexp_par <- sapply(allcombexp, function(x)
            apply(x[,c(1,3)],1,min) )
        allcombexp_edg <- sapply(allcombexp,function(x)
            apply(x[,c(2,4)],1,max) )
        allcombexp_which <- apply(abs(allcombexp_par-allcombexp_edg),1,which.max)
        expect1 <- purrr::map_dbl(
            seq_len(length(allcombexp_which)), function(x) {
                allcombexp[[allcombexp_which[x]]][x,5]
            })

        # # edge proportion values extracted from ep for edges between
        # # cpop's parents and grandparents.
        # some: pes <- ep[,paste0(pparen[[cpop]], "_", cpop),drop=FALSE]
        pedges <- ep_[,pnames,drop=FALSE]

        # # using the above, get expected proportion; see formula in
        # # https://www.biorxiv.org/content/10.1101/837765v2
        # some: expect1 <- apply(pedges,1,min) * apply(parent,1,max)
        return(expect1)
    }, no_cores=no_cores, prll=length(cpind)>1000))
    expecp[is.nan(expecp)] <- 0
    colnames(expecp) <- cells[cpind]


    ## infer the expected prop of all other cell pops --------------------------
    # this isn't necessary, we could've just done it all above,
    # but this saves time; see same paper as above link
    expec0 <- as.matrix(cbind(expe1,expecp))

    exp1 <- expec0
    exp1[as.matrix(exp1)<0] <- 0
    # exp1 <- cbind(expe1,expecp[,match(cells,colnames(expecp)),drop=FALSE])
    return(exp1) ### EXPECTED PROPORTION
}



#' @title Converts cell counts into cumulated cell counts.
#' @description Converts the cell counts in a flowGraph object into
#'  cumulated cell counts; this is optional and can be done only for there is
#'  more than one threshold for one or more markers.
#'  This should also only be ran when initializing a flowGraph object as
#'  converting back and forth is computationally expensive.
#'  If the user is interested in seeing non- and cumulated counts, we recommend
#'  keeping two flowGraph objects, one for each version.
#'  This function simply converts e.g. the count of A+ or A++ into
#'  the sum of count of A+, A++, and A+++ or A++, and A+++.
#' @param fg flowGraph object.
#' @param no_cores An integer indicating how many cores to parallelize on.
#' @return flowGraph object with cumulated counts.
#' @details \code{fg_feat_cumsum} returns the given flowGraph object with an
#'  adjusted count feature. As in our example,
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=no_cores)
#'
#'  fg <- flowGraph:::fg_feat_cumsum(fg, no_cores=no_cores)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[Matrix]{Matrix}}
#' @rdname fg_feat_cumsum
#' @importFrom future plan multiprocess
#' @importFrom stringr str_split
#' @importFrom Matrix Matrix
#' @importFrom furrr future_map
#' @importFrom purrr map
fg_feat_cumsum <- function(fg, no_cores) {
    if (no_cores>1) future::plan(future::multiprocess)

    # check if already cumsum
    if (fg_get_etc(fg)$cumsumpos) return(fg)

    fg_graph <- fg_get_graph(fg)

    # check if do-able (there exists multple ++)
    if (!any(grepl("3",fg_graph$v$phenocode))) return(fg)

    mc <- as.matrix(fg_get_feature(fg, "node", "count"))
    meta_cell <- fg_graph$v
    meta_cell_grid <-
        do.call(rbind, stringr::str_split(meta_cell$phenocode,""))
    meta_cell_grid <- apply(meta_cell_grid, 2, as.numeric)
    meta_cell_grid_TF1 <- apply(meta_cell_grid,2, function(x) {
        xmax <- max(x)
        if (xmax==2) return(rep(FALSE,length(x)))
        x>1 & x<xmax
    })
    cany1 <- which(apply(meta_cell_grid_TF1, 2, any))
    for (marker in cany1) {
        coldo <- which(meta_cell_grid_TF1[,marker])
        coldo <- coldo[order(meta_cell_grid[coldo,marker],
                             decreasing=TRUE)]
        mc[,coldo] <- do.call(
            cbind, fpurrr_map(coldo, function(j) {
                mcgi <- mcgis <- mcgip <- meta_cell_grid[j,]
                mcgis[marker] <- mcgis[marker]+1
                jsib <- meta_cell$phenocode==paste0(mcgis,collapse="")
                mc[,j] + mc[,jsib]
            }, no_cores=no_cores, prll=length(coldo)>1000))

    }
    fg@feat$node$count_original <- fg_get_feature(fg, "node", "count")
    fg@feat$node$count <- mc
    fg@etc$cumsumpos <- TRUE
    if (length(fg_get_feature_all(fg)$node)>2 |
        length(fg_get_feature_all(fg)$edge)>0)
        warning("IMPORTANT: see function fg_rm_features to remove
                node features (other than count), and edge features.")
    return(fg)
}


#' @title Normalizes all features for class.
#' @description For each class label in column \code{class} of \code{meta},
#'  \code{fg_feat_mean_class}
#'  takes the column mean of the rows in the given feature matrices
#'  (as specified in \code{node_features} and \code{edge_features})
#'  associated with that class;
#'  it then takes the difference point by point between these means and
#'  the original rows for that class.
#' @return A numeric matrix whose dimensions equate to that of the input
#'  and whose values are normalized per class.
#' @title Normalizes matrix values in a flowGraph object by class.
#' @description FUNCTION_DESCRIPTION
#' @param fg PARAM_DESCRIPTION
#' @param class a column name in \code{fg_get_meta(fg)} indicating the
#'  meta data that should be used as the class label of each sample
#'  while conudcting normalization.
#' @param no_cores An integer indicating how many cores to parallelize on.
#' @param node_features A string vector indicating the node features to
#'  perform normalization on; set as \code{NULL} to normalize all.
#' @param edge_features A string vector indicating the edge features to
#'  perform normalization on; set as \code{NULL} to normalize all.
#' @return flowGraph object with normalized features.
#' @details For all features in the given \code{flowGraph} object and for
#'  each class label in column \code{class} of \code{meta},
#'  \code{fg_feat_mean_class}. It takes the column mean of the rows in the given
#'  feature matrices
#'  (as specified in \code{node_features} and \code{edge_features})
#'  associated with that class;
#'  it then takes the difference point by point between these means and
#'  the original rows for that class.
#'  \code{fg_feat_mean_class}
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=no_cores)
#'
#'  fg <- fg_feat_mean_class(fg, class="class", node_features="count",
#'                         no_cores=no_cores)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#' @rdname fg_feat_mean_class
#' @export
#' @importFrom future plan multiprocess
#' @importFrom furrr future_map
fg_feat_mean_class <- function(
    fg, class, no_cores=1, node_features=NULL, edge_features=NULL
) {
    if (no_cores>1) future::plan(future::multiprocess)

    fg_meta <- fg_get_meta(fg)

    if (!class%in%colnames(fg_meta))
        stop("invalid class name, choose one from fg_get_meta(fg)\n")

    if (length(unique(fg_meta[,class]))<2) return(fg)
    label <- paste0("_MEAN_",class)

    fg_feat <- fg_get_feature_all(fg)

    if (is.null(node_features)) node_features <- names(fg_feat$node)
    node_features <- node_features[node_features%in%names(fg_feat$node)]
    node_features <- node_features[!grepl("MEAN", node_features)]
    if (length(node_features)>0) {
        fg_nodefs0 <- fg_feat$node[node_features]
        fg_nodefs <- fpurrr_map(fg_nodefs0, function(m0)
            mean_diff(as.matrix(m0), fg_meta[,class]),
            no_cores=no_cores, prll=length(fg_nodefs0)>1000)
        names(fg_nodefs) <- paste0(names(fg_nodefs),label)
        fg@feat$node <- append(fg_feat$node, fg_nodefs)
    }

    if (is.null(edge_features)) edge_features <- names(fg_feat$edge)
    edge_features <- edge_features[edge_features%in%names(fg_feat$edge)]
    edge_features <- edge_features[!grepl("MEAN", edge_features)]
    if (length(edge_features)>0) {
        fg_edgefs0 <- fg_feat$edge[edge_features]
        fg_edgefs <- fpurrr_map(fg_edgefs0, function(m0)
            mean_diff(as.matrix(m0), fg_meta[,class]),
            no_cores=no_cores, prll=length(fg_nodefs0)>1000)
        names(fg_nodefs) <- paste0(names(fg_edgefs),label)
        fg@feat$edge <- append(fg_feat$edge, fg_edgefs)
    }

    fg@feat_desc <- fg_get_feature_desc(fg, re_calc=TRUE)
    return(fg)
}



