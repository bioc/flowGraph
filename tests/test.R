library(flowGraph)
# library(ftf***)
library(testthat)

# prepare parallel backend
no_cores <- 1#parallel::detectCores()-1
# future::plan(future::multiprocess)

## create Phenotypes data ----------------------

data(fg_data_pos30)

# # define marker and cell count
# celln <- 10000
# markern <- 3
# markers <- LETTERS[1:markern]
#
# # define marker thresholds
# cvd <- rnorm(celln,2,1)
# p50 <- quantile(cvd, .5)
# thres <- furrr::future_map(markers, function(x) p50)
# names(thres) <- markers
#
# # generate Phenotypes list
samplen <- 10
# ftl <- furrr::future_map(1:samplen, function(i) {
#     # make flow frame
#     f <- new("flowFrame")
#     f@exprs <- matrix(rnorm(celln*markern,2,1), nrow=celln)
#     colnames(f@exprs) <- markers
#
#     # marker indices in flow frame
#     ci <- c(1:ncol(f@exprs))
#     names(ci) <- colnames(f@exprs)
#
#     # modify experiment samples such that ABC increases by 50%
#     if (i>(samplen/2)) {
#         ap <- f@exprs[,1]>thres[[1]]
#         bp <- f@exprs[,2]>thres[[2]]
#         cp <- f@exprs[,3]>thres[[3]]
#         tm <- sum(ap & bp & cp)/2
#         f@exprs <- rbind(f@exprs, f@exprs[sample(which(ap & bp & cp),tm),])
#     }
#
#     # make Phenotypes
#     ftf***(Frame=f, PropMarkers=ci, MarkerNames=colnames(f@exprs),
#              MaxMarkersPerPop=markern, PartitionsPerMarker=2, Thresholds=thres,
#              Methods='Thresholds', verbose=FALSE, MemLimit=60)#@CellFreqs
# })

meta_file <- data.frame(
    id=1:samplen,
    class=append(rep("control", samplen/2), rep("exp", samplen/2)),
    stringsAsFactors=FALSE
)


## 00 helpers ----------------------------------

# context("00_helpers")

# test_that("ncores", {
#     expect_equal(ncores(1), 1)
#     expect_equal(ncores(Inf), parallel::detectCores())
# })

# test_that("loop_ind_f", {
#     expect_equal(length(loop_ind_f(1,10)), 1)
#     expect_equal(length(loop_ind_f(rep(1,2),10)), 2)
#     expect_equal(length(loop_ind_f(rep(1,10),10)), 10)
#     expect_equal(length(loop_ind_f(rep(1,15),10)), 10)
# })


## 00 flowgraph constructor --------------------

context("00_flowGraph_constructor")

test_that("errors", {
    expect_error(flowGraph(NA))#, "make sure input type is correct")
    expect_error(flowGraph(fg_data_pos30$count[1:2,],
                           meta=fg_data_pos30$meta[1:5,]))
    expect_error(flowGraph(fg_data_pos30$count[1:2,],
                           meta=fg_data_pos30$meta[1:2,-1]))
})

contexts <- c("vector", "matrix") #"phenotype", "phenotype list",
for (context_no in contexts) {
    if (context_no == "phenotype")
        fg <- flowGraph(ftl[[1]], no_cores=no_cores)
    if (context_no == "phenotype list")
        fg <- flowGraph(ftl[1:2], no_cores=no_cores)
    if (context_no == "vector")
        fg <- flowGraph(fg_data_pos30$count[1,], no_cores=no_cores)
    if (context_no == "matrix")
        fg <- flowGraph(fg_data_pos30$count[1:2,], no_cores=no_cores)

    test_that(context_no, {
        expect_equal(nrow(fg@meta), nrow(fg@feat$node$count))
        expect_equal(nrow(fg@graph$v), ncol(fg@feat$node$count))
        expect_equal(nrow(fg@graph$e), ncol(fg@feat$edge$prop))

        expect_equal(nrow(fg@feat_desc$node), length(fg@feat$node))
        expect_equal(nrow(fg@feat_desc$edge), length(fg@feat$edge))

        a = lapply(fg@feat$node, function(x)
            expect_equal(nrow(x), nrow(fg@feat$node$count)) )
        a = lapply(fg@feat$edge, function(x)
            expect_equal(nrow(x), nrow(fg@feat$node$count)) )
    })
}

test_that("graph", {
    expect_setequal(fg@edge_list$child[["A++B++"]],
                 c("A++B++C-", "A++B++C+", "A++B++D-", "A++B++D+"))
    expect_setequal(fg@edge_list$parent[["A++B++"]],
                 c("A++", "B++"))
    expect_setequal(fg@graph$e$from[fg@graph$e$to=="A++B++"],
                 c("A++", "B++"))
})


context("02_flowGraph_features")

test_that("prop", {
    expect_equal(all(fg@feat$node$prop[,1]==1), TRUE)
})

# fg <- flowGraph(fg_data_pos30$count[1:2,], no_cores=no_cores)
fg_ <- flowGraph(fg_data_pos30$count[1:2,], no_cores=no_cores, cumsumpos=TRUE)

test_that("cumsumpos", {
    expect_equal(fg_@feat$node$count[1,"A+"],
                 sum(fg@feat$node$count[1,c("A+","A++")]))
    expect_equal(fg_@feat$node$count[1,"B++"],
                 sum(fg@feat$node$count[1,c("B++","B+++")]))
    expect_equal(fg_@feat$node$count[1,"B++C+"],
                 sum(fg@feat$node$count[1,c("B++C+","B+++C+")]))
    expect_equal(fg_@feat$node$count[1,"A++B++"],
                 sum(fg@feat$node$count[1,c("A++B++","A++B+++")]))
})


context("03_flowGraph_summary")

# fg <- flowGraph(ftl, meta=meta_file, no_cores=no_cores)

fg <- fg_summary(fg, no_cores=no_cores, class="class", label1="control",
                overwrite=FALSE, test_name="wilcox", diminish=FALSE)

fg <- fg_summary(fg, no_cores=no_cores, class="class", label1="control",
                overwrite=FALSE, test_name="wilcox_diminish", diminish=TRUE)

test_that("summary", {
    expect_equal(length(fg@summary$node), nrow(fg@summary_desc$node))
})

fg <- fg_clear_summary(fg)

test_that("summary clear", {
    expect_equal(fg@summary, list())
    expect_equal(fg@summary_desc, list())
})



