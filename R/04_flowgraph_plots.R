## FG PLOT ---------------------------------------------------------------------

#' @title Creates a cell hierarchy plot.
#' @description Creates a cell hierarchy plot given a flowGraph object. If a path is not provided for \code{fg_plot} to save the plot, please use \code{plot_gr} to view plot given the output of \code{fg_plot}.
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
#'  \code{feature} (feature name), \code{test_name} (summary statistic name),
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
#' @param show_nodes_edges A logical vector indicating which nodes/edges (type)
#'  to show in the plot; if this is not specified, only nodes/edges with
#'  significant summary statistics will be shown.
#' @param label_max An integer specifying the maximum number of nodes to label.
#' @param p_thres A double indicating a summary statistic threshold
#'  e.g. if we are plotting a T test summary statistic, we can set the threshold
#'  to .05; nodes with a p-value greater than .05 will not be plotted.
#' @param filter_adjust0 A numeric variable indicating what percentage of
#'  SpecEnr values compared (minimum) should be not close to 0.
#'  Set to 1 to not conduct filtering.
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
#' @param node_labels A string vector indicating which node feature(s)
#'  should be used to label a node.
#'  We recommend keeping the length of this vector to below 2.
#'  Set to "NONE"/NULL if no p-value labels are needed.
#' @param summary_fun A function that takes in a matrix and outputs a
#'  vector the same length as the number of columns this matrix has;
#'  see \code{\link[flowGraph]{fg_summary}}.
#' @param adjust_custom A function or a string indicating the
#' test adjustment method to use.
#'  If a string is provided, it should be one of
#'  \code{c("holm", "hochberg", "hommel",
#'  "bonferroni", "BH", "BY", "fdr", "none")} (see \code{p.adjust.methods}).
#'  If a function is provided, it should take as input
#'  a numeric vector and output the
#'  same vector adjusted.
#' @param layout_fun A string representing a function from the
#'  \code{igraph} package that
#'  indicates what layout should be used if a cell hierarchy is to be ploted;
#'  all such functions have prefix \code{layout_}. Only specify if different
#'  from the default one already calculated in the \code{fg} flowGraph object
#'  given.
#' @param show_bgedges A logical variable indicating whether or not
#'  edges not specified for plotting should be plotted as
#'  light grey in the background.
# #' @param mod_layout A logical variable indicating whether or not to re-space
# #'  nodes if the layout is a cell hierarchy and \code{show_bgedges=FALSE}.
#' @param main A string or the title of the plot; if left as \code{NULL},
#'  a default title will be applied.
#' @param interactive A logical variable indicating whether the plot should be
#'  an interactive plot; see package \code{ggiraph}.
#' @param visNet_plot A logical variable indicating if an interactive plot is
#'  chosen, if function should output a visNetwork plot; if set to \code{FALSE},
#'  ggplot's girafe will be used instead.
#' @param path A string indicating the path to where the function should save
#'  the plot; leave as \code{NULL} to not save the plot. Static plots are saved
#'  as PNG, interactive plots are saved as HTML.
#' @param width A numeric variable specifying, in inches,
#'  what the plot width should be.
#' @param height A numeric variable specifying, in inches,
#'  what the plot height should be.
#' @return A list of nodes and edges for plotting with the \code{plot_gr}
#'  function. Other elements in this list include \code{show_bgedges},
#'  which has the same value as parameter \code{show_bgedges}, and \code{main},
#'  the title of the plot.
#' @details \code{fg_plot} takes a flowGraph object as input and returns the
#'  \code{graph} slot of the given object with additional columns to serve as
#'  input into \code{\link[flowGraph]{plot_gr}} for plotting using functions
#'  in the \code{ggplot2} package. Users can choose to save a PNG version of
#'  the plot by filling out the \code{path} parameter with a full path to the
#'  PNG plot. In addition to specifying columns added from
#'  \code{\link[flowGraph]{ggdf}}, \code{fg_plot} also adds label column(s)
#'  whose values serve as labels in the interactive version of the plot.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos15)
#'  fg <- flowGraph(fg_data_pos15$count, class=fg_data_pos15$meta$class,
#'                  no_cores=no_cores)
#'
#'  gr <- fg_plot(fg, type="node", index=1, label_max=30,
#'    show_nodes_edges=NULL, p_thres=.01, node_labels=c("prop", "expect_prop"),
#'    path=NULL) # set path to a full path to save plot as a PNG
#'  # plot_gr(gr)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{get_phen_meta}}
#'  \code{\link[flowGraph]{ggdf}}
#'  \code{\link[flowGraph]{plot_gr}}
#'  \code{\link[flowGraph]{fg_get_feature}}
#'  \code{\link[flowGraph]{fg_get_summary}}
#' @rdname fg_plot
#' @export
#' @importFrom htmlwidgets saveWidget
#' @importFrom ggplot2 ggsave
#' @importFrom visNetwork visSave
fg_plot <- function(
    fg, type="node", index=1, summary_meta=NULL, adjust_custom="byLayer",
    show_nodes_edges=NULL,
    label_max=30, p_thres=.05, filter_adjust0=1, filter_es=0,
    filter_btwn_tpthres=1, filter_btwn_es=0,
    # show_nodes_edges: a vector of true false, length equals to number of nodes
    node_labels=c("prop", "expect_prop"),
    summary_fun=base::colMeans, layout_fun=NULL,
    show_bgedges=TRUE,
    # mod_layout=FALSE,
    main=NULL, interactive=FALSE, visNet_plot=TRUE,
    path=NULL, width=9, height=9
) {
    # width in inches feat must be list with names node and edge
    # cellhierarchy plots p values

    type <- base::match.arg(type, c("node", "edge"))

    gr <- ggdf(fg@graph)

    # change layout if needed
    if (!base::is.null(layout_fun)) {
        if (layout_fun!=fg@plot_layout)
            gr <- set_layout_graph(gr, layout_fun) # layout cell hierarchy
    }

    # get summary statistics
    summary_meta <- base::unlist(fg@summary_desc[[type]][
        flowGraph:::fg_get_summary_index(fg, type, index, summary_meta),])
    if (!base::grepl("SpecEnr",base::unlist(summary_meta["feat"])))
        filter_adjust0 <- 1
    pms <- fg_get_summary(
        fg, type, index, summary_meta, adjust_custom=adjust_custom,
        summary_fun=summary_fun, default_p_thres=p_thres,
        filter_adjust0=filter_adjust0, filter_es=filter_es,
        filter_btwn_tpthres=filter_btwn_tpthres, filter_btwn_es=filter_btwn_es
    )
    id1 <- pms$id1
    id2 <- pms$id2
    m1 <- pms$m1
    m2 <- pms$m2
    p <- pms$values

    # show_nodes_edges: which nodes/edges to show
    sn_ <- p < p_thres
    if (!base::is.null(show_nodes_edges))
        if (base::length(show_nodes_edges)==base::length(m1) &
            base::all(show_nodes_edges%in%c(TRUE,FALSE)) |
            base::length(show_nodes_edges)>1)
            sn_ <- show_nodes_edges
    if (base::sum(sn_)==0) {
        warning("no significant nodes found, maybe increase p_thres?")
        return(NULL)
    }

    if (type=="node") {
        gr$v$v_ind <- gr$v$label_ind <- sn_
        gr$e$e_ind <- gr$e$from%in%gr$v$phenotype[sn_] &
            gr$e$to%in%gr$v$phenotype[sn_]
        if (base::sum(sn_)>label_max) {
            gr$v$label_ind <- FALSE
            gr$v$label_ind[utils::head(base::order(p), label_max)] <- TRUE
        }
    } else {
        gr$e$e_ind <- sn_
        gr$v$v_ind <- gr$v$label_ind <-
            gr$v$phenotype%in%base::unlist(gr$e[sn_,c("to","from")])
        if (base::sum(gr$v$v_ind)>label_max) {
            gr$v$label_ind <- FALSE
            po <- base::order(p); po <- po[po %in% base::which(sn_)]
            while (base::sum(gr$v$label_ind) < label_max |
                   base::length(po)==0) {
                po1 <- gr$v$phenotype%in%gr$e[po[1],c("to","from")]
                gr$v$label_ind[po1] <- TRUE
                po <- po[-1]
            }
        }
    }

    # plot title
    if (base::is.null(main))
        main <- base::paste0(
            "cell hierarchy plot.\n",
            "- ", type, " feature: ", summary_meta["feat"], ".\n",
            "- p-value: ", summary_meta["test_name"],
            "\n- comparing ",
            summary_meta["class"], ", labels: ", summary_meta["label1"], " & ", summary_meta["label2"], ".")

    # create node labels
    if (type=="node") {
        gr$v$label <- gr$v$label_long <- gr$v$phenotype
        if (node_labels[1]!="NONE")
            gr$v$label <- gr$v$label_long <-
                base::paste0(gr$v$label, " (",base::signif(p,3),")")
        if (!"NONE"%in%node_labels) {
            node_labels <-
                node_labels[node_labels%in%base::names(fg@feat[[type]])]
            if (base::length(node_labels)==0) node_labels <- NULL
            if (type=="node") node_labels <-
                    base::append(summary_meta["feat"], node_labels)
            if (!base::is.null(node_labels))
                node_labels <- node_labels[!base::is.null(node_labels)]
            if (!base::is.null(node_labels))
                node_labels <- node_labels[!base::duplicated(node_labels)]
            if (!base::is.null(node_labels))
                for (label in node_labels) {
                    vls <- base::list(
                        m1=fg_get_feature_means(fg, "node", label, id=id1),
                        m2=fg_get_feature_means(fg, "node", label, id=id2))
                    vl <- base::paste0(base::round(vls$m1,3), "/",
                                       base::round(vls$m2,3))
                    gr$v$label_long <-
                        base::paste0(gr$v$label_long, "\n", label, ": ", vl)
                    gr$v$label <- base::paste0(gr$v$label, " ", vl)
                }
        }
    } else {
        gr$v$label <- gr$v$label_long <- base::paste0(gr$v$phenotype)
    }

    # node colour, size
    if (type=="node") {
        gr$v$colour <- m2-m1
        gr$v$size <- -log(p)
        gr$v$size[!base::is.finite(gr$v$size)] <-
            base::max(gr$v$size[base::is.finite(gr$v$size)])
        gr$v$size[base::is.nan(gr$v$size)] <- 0
        ew1 <- fg_get_feature_means(fg, type="edge", feature="prop", id=pms$id1)
        ew2 <- fg_get_feature_means(fg, type="edge", feature="prop", id=pms$id2)
        gr$e$colour <- ew1 - ew2
        gr$e$colour1 <- ew1
        gr$e$colour2 <- ew2
    } else {
        gr$e$colour <- m2-m1
        gr$e$size <- -log(p)
        gr$e$size[!base::is.finite(gr$e$size)] <-
            base::max(gr$e$size[base::is.finite(gr$e$size)])
        gr$e$size[base::is.nan(gr$e$size)] <- 0
    }

    # # space out nodes again if no background edges
    # if (!show_bgedges & fg@plot_layout=="layout.reingold.tilford" & mod_layout) {
    #     gr_vxy <- space_hierarchy(gr$v[,c("y","x")])
    #     gr$v[,c("x","y")] <- gr_vxy[,c(2,1)]
    # }

    # plot and save
    gp <- plot_gr(gr, main=main, show_bgedges=show_bgedges,
                  interactive=interactive,
                  visNet_plot=visNet_plot)
    suppressMessages({
        if (!base::is.null(path) & !interactive)
            ggplot2::ggsave(base::ifelse(
                base::grepl("[.]png$",path, ignore.case=TRUE),
                path, base::paste0(path, ".png")),
                            plot=gp, scale=1, width=width, height=height,
                            units="in", dpi=600, limitsize=TRUE)
        if (!base::is.null(path) & interactive & !visNet_plot)
            htmlwidgets::saveWidget(
                gp, base::ifelse(base::grepl("[.]html$",path, ignore.case=TRUE),
                       path, base::paste0(path, ".html")))
        if (!base::is.null(path) & interactive & visNet_plot)
            visNetwork::visSave(
                gp, file=ifelse(base::grepl("[.]html$",path, ignore.case=TRUE),
                                path, base::paste0(path, ".html")),
                background="white")
    })
    if (base::is.null(path)) message("use function plot_gr to plot fg_plot output")
    gr$main <- main
    gr$show_bgedges <- show_bgedges
    gr$visNet_plot <- visNet_plot
    gp
    return(gr)
}




## graph plot functions

#' @title Prepares a given node and edge graph list for plotting.
#' @description Prepares a given node and edge graph list
#'  for plotting by function plot_gr;
#'  do not use this function on its own.
#' @param gr0 A list containing data frames \code{e} and \code{v}.
#' @return A list containing data frames \code{e} and \code{v}, each
#'  with additional meta data column.
#' @details code{ggdf} adds to the data frames \code{v} and \code{e} in slot
#'  \code{graph} from a \code{flowGraph} object specifying plotting options as
#'  required by \code{\link[flowGraph]{plot_gr}}:
#'  \itemize{
#'    \item{\code{v}}
#'    \itemize{
#'      \item{\code{size}: a numeric indicating node size.}
#'      \item{\code{colour}: a numeric or string indicating node colour.}
#'      \item{\code{label}: a string indicating the label of a node.}
#'      \item{\code{label_long}: a string indicating teh long label of a node;
#'       used in interactive plots in \code{\link[flowGraph]{plot_gr}}}.
#'      \item{\code{label_ind}: a vector of logical variables indicating which
#'       nodes to add a label to in a static plot.}
#'      \item{\code{v_ind}: a vector of logical variables indicating which
#'       nodes to plot.}
#'    }
#'    \item{\code{e}}
#'    \itemize{
#'      \item{\code{colour}: a numeric or string indicating edge colour.}
# #'      \item{\code{size}: a numeric or string indicating edge size.}
#'      \item{\code{e_ind}: a vector of logical variables indicating which
#'       edges to plot.}
#'    }
#'  }
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos30)
#'  fg <- flowGraph(fg_data_pos30$count, class=fg_data_pos30$meta$class,
#'                  prop=FALSE, specenr=FALSE,
#'                  no_cores=no_cores)
#'
#'  gr_ <- ggdf(fg_get_graph(fg))
#'  head(gr_$v)
#'  head(gr_$e)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{get_phen_meta}}
#'  \code{\link[flowGraph]{plot_gr}}
#' @rdname ggdf
#' @export
ggdf <- function(gr0) {
    base::list(e=base::data.frame(gr0$e, colour=0, size=1, e_ind=FALSE),
               v=base::data.frame(gr0$v,
                                  size=1, colour=0,
                                  #sizeb=1, colourb="", fill="",
                                  label=gr0$v$phenotype,
                                  label_ind=FALSE, v_ind=FALSE))
}

#' @title Plots a cell hierarchy.
#' @description Plots a cell hierarchy given the output from \code{fg_plot}, a list of nodes and edges.
#' @param gr A list containing data frames \code{e} and \code{v}.
#' @param main A string containing the plot title. If this is set to NULL, the
#'  function will look for a plot title in the \code{main} slot of \code{gr};
#'  otherwise, this defaults to "".
#' @param show_bgedges A logical variable indicating whether or not
#'  edges not specified for plotting should be plotted as light grey
#'  in the background. If this is \code{NULL}, the function will look for a
#'  \code{show_bgedges} in the \code{show_bgedges} slot of \code{gr};
#'  otherwise, this defaults to \code{TRUE}.
#' @param colour_palette A colour palette e.g. the default palette if the user
#'  sets this to \code{NULL} is \code{rev(RColorBrewer::brewer.pal(10,"RdBu"))}.
#' @param label_coloured A logical indicating whether to colour the node
#'  labels using the same colours as the nodes in the non-interactive plot.
#' @param shiny_plot A logical indicating whether this plot is made for shiny;
#'  users don't need to change this.
#' @param interactive A logical variable indicating whether the plot should be
#'  an interactive plot; see package \code{ggiraph}.
#' @param visNet_plot A logical variable indicating if an interactive plot is
#'  chosen, if function should output a visNetwork plot; if set to \code{FALSE},
#'  ggplot's girafe will be used instead.
#' @param colour_edges A logical variable indicating whether to colour edges if
#'  plotting a node feature summary.
#' @param ... Other parameters for \code{ggplot} if \code{interactive}
#'  is set to \code{FALSE}; other parameters for \code{plot_ly}
#'  if \code{interactive} is set to \code{TRUE}.
#' @return  A \code{ggplot} object if \code{interactive} is set to \code{FALSE};
#'  a \code{ggiraph} object if \code{interactive} is set to \code{TRUE}.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos15)
#'  fg <- flowGraph(fg_data_pos15$count, class=fg_data_pos15$meta$class,
#'                  no_cores=no_cores)
#'
#'  # fg <- fg_summary(fg, no_cores=no_cores, class="class", control="control",
#'  #                  overwrite=FALSE, test_name="t_byLayer", diminish=FALSE)
#'
#'  gr_summary <- fg_plot(
#'    fg, type="node", p_thres=.05, show_bgedges=TRUE,
#'    path=NULL) # set path to a full path to save plot as a PNG
#'
#'  plot_gr(gr_summary, main=gr_summary$main, show_bgedges=TRUE)
#'
#'  plot_gr(gr_summary, main=gr_summary$main, show_bgedges=TRUE, interactive=TRUE)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_plot}}
#'  \code{\link[flowGraph]{get_phen_meta}}
#'  \code{\link[flowGraph]{ggdf}}
#'  \code{\link[flowGraph]{fg_get_feature}}
#'  \code{\link[flowGraph]{fg_get_summary}}
# #'  \code{\link[plotly]{plot_ly}}
#' @rdname plot_gr
#' @export
#' @importFrom ggplot2 ggplot aes scale_x_continuous scale_y_continuous
#'  theme element_blank ggtitle scale_fill_brewer geom_segment geom_point
#'  scale_colour_gradientn labs
#' @importFrom ggiraph girafe geom_point_interactive geom_segment_interactive
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggrepel geom_label_repel
#' @importFrom visNetwork visNetwork visNodes visEdges visOptions visInteraction visHierarchicalLayout visEvents
#' @importFrom purrr map
plot_gr <- function(
    gr, main=NULL, show_bgedges=TRUE, colour_palette=NULL, label_coloured=TRUE,
    shiny_plot=FALSE, interactive=FALSE, visNet_plot=TRUE, colour_edges=FALSE,
    ...
) {
    # gr_v: name x y label size colour sizeb colourb
    # gr_e: from to from.x from.y to.x to.y colour

    # place holder; all variables specified in plotting functions
    # refers to columns names in the given data frame
    size <- 0

    gr_v <- gr$v
    gr_e <- gr$e

    if (base::is.null(main))
        if (base::is.null(gr$main)) {
            main <- ""
        } else {
            main <- base::paste0("(",base::sum(gr_v$v_ind),"/",
                                 base::nrow(gr_v),") ", gr$main)
        }
    if (!base::is.null(gr$show_bgedges)) show_bgedges <- gr$show_bgedges
    if (!base::is.null(gr$interactive)) interactive <- gr$interactive
    if (!base::is.null(gr$visNet_plot)) show_bgedges <- gr$show_bgedges

    if (shiny_plot) {
        interactive <- TRUE
        if (!visNet_plot) {
            max_x <- base::max(gr$v$x)
            gr_v$x <- gr$v$y
            gr_v$y <- max_x-gr$v$x
            gr_e$from.x <- gr$e$from.y
            gr_e$from.y <- max_x-gr$e$from.x
            gr_e$to.x <- gr$e$to.y
            gr_e$to.y <- max_x-gr$e$to.x
        }
    }

    if(base::is.null(colour_palette))
        colour_palette <- c('blue','cyan','yellow','red')
    # colour_palette <- rev(RColorBrewer::brewer.pal(7,"RdBu"))

    # prepare base plot
    gp <- ggplot2::ggplot(gr_v, ggplot2::aes(x=x, y=y, ...)) +
        ggplot2::scale_x_continuous(expand=c(0,1)) +  # expand x limits
        ggplot2::scale_y_continuous(expand=c(0,1)) + # expand y limits
        # theme_bw()+  # use the ggplot black and white theme
        ggplot2::theme(
            axis.text.x=ggplot2::element_blank(),  # rm x-axis text
            axis.text.y=ggplot2::element_blank(), # rm y-axis text
            axis.ticks=ggplot2::element_blank(),  # rm axis ticks
            axis.title.x=ggplot2::element_blank(), # rm x-axis labels
            axis.title.y=ggplot2::element_blank(), # rm y-axis labels
            panel.background=ggplot2::element_blank(),
            panel.border=ggplot2::element_blank(),
            panel.grid.major=ggplot2::element_blank(),  #rm grid labels
            panel.grid.minor=ggplot2::element_blank(),  #rm grid labels
            plot.background=ggplot2::element_blank()) +
        ggplot2::ggtitle(main) +
        ggplot2::scale_colour_gradientn(colours=colour_palette) +
        ggplot2::labs(size="-ln(p-value)",
                      col="mean feat diff")

    if (show_bgedges)  # keep greyed out edges on
        gp <- gp +
        ggplot2::geom_segment(
            data=gr_e[!gr_e$e_ind,], color="grey",
            ggplot2::aes(x=from.x,xend=to.x, y=from.y,yend=to.y))

    if (!interactive) {
        if (!base::all(gr_e$colour==gr_e$colour[1]) &
            !base::all(gr_e$size==gr_e$size[1])) {
            gp <- gp +
                ggplot2::geom_segment(
                    data=gr_e[gr_e$e_ind,],
                    ggplot2::aes(x=from.x, xend=to.x, y=from.y, yend=to.y,
                                 color=colour, size=size)) +
                ggplot2::geom_point(
                    data=gr_v[gr_v$v_ind,], colour="grey50",
                    ggplot2::aes(x=x,y=y))
        } else {
            # gr_e$colour <- apply(gr_e[,c("from","to")], 1, function(x)
            #     exp(-gr_v$size[gr_v$phenotype==x[2]]) -
            #         exp(-gr_v$size[gr_v$phenotype==x[1]]))
            if (colour_edges) {
                gp <- gp +
                    ggplot2::geom_segment(
                        data=gr_e[gr_e$e_ind,], # colour="grey50",
                        ggplot2::aes(x=from.x, xend=to.x, y=from.y, yend=to.y,
                                     colour=colour))
            } else {
                gp <- gp +
                    ggplot2::geom_segment(
                        data=gr_e[gr_e$e_ind,], colour="grey50",
                        ggplot2::aes(x=from.x, xend=to.x, y=from.y, yend=to.y))
            }
            gp <- gp +
                ggplot2::geom_point(
                    data=gr_v[gr_v$v_ind,],
                    ggplot2::aes(x=x,y=y, color=colour, size=size))
            if (label_coloured) {
                gp <- gp +
                    ggrepel::geom_label_repel(
                        data=gr_v[gr_v$label_ind,],
                        ggplot2::aes(x=x, y=y,label=label, color=colour),
                        nudge_x=-.1, direction="y", hjust=1,
                        segment.size=0.2,force=1.5)

            } else {
                gp <- gp +
                    ggrepel::geom_label_repel(
                        data=gr_v[gr_v$label_ind,],
                        ggplot2::aes(x=x, y=y,label=label),
                        nudge_x=-.1, direction="y", hjust=1,
                        segment.size=0.2,force=1.5)
            }
        }
    }
    else if (!visNet_plot) {
        if (!base::all(gr_e$colour==gr_e$colour[1]) &
            !base::all(gr_e$size==gr_e$size[1])) {
            gp <- gp +
                ggplot2::geom_segment(
                    data=gr_e[gr_e$e_ind,],
                    ggplot2::aes(x=from.x, xend=to.x, y=from.y, yend=to.y,
                                 color=colour, size=size,
                                 tooltip=base::paste0(from, " > ", to))) +
                ggplot2::geom_point(
                    data=gr_v[gr_v$v_ind,], colour="grey50",
                    ggplot2::aes(x=x,y=y))
        } else {
            # gr_e$colour <- apply(gr_e[,c("from","to")], 1, function(x)
            #     exp(-gr_v$size[gr_v$phenotype==x[2]]) -
            #         exp(-gr_v$size[gr_v$phenotype==x[1]]))
            if (colour_edges) {
                gp <- gp +
                    ggplot2::geom_segment(
                        data=gr_e[gr_e$e_ind,], # colour="grey50",
                        ggplot2::aes(x=from.x, xend=to.x, y=from.y, yend=to.y,
                                     colour=colour))
            } else {
                gp <- gp +
                    ggplot2::geom_segment(
                        data=gr_e[gr_e$e_ind,], colour="grey50",
                        ggplot2::aes(x=from.x, xend=to.x, y=from.y, yend=to.y))
            }

            if (shiny_plot) {
                gp <- gp + ggiraph::geom_point_interactive(
                    data=gr_v[gr_v$v_ind,],
                    ggplot2::aes(x=x,y=y, color=colour, size=size,
                                 tooltip=label_long, data_id=phenotype))
            } else {
                gp <- gp + ggiraph::geom_point_interactive(
                    data=gr_v[gr_v$v_ind,],
                    ggplot2::aes(x=x,y=y, color=colour, size=size,
                                 tooltip=label_long, data_id=phenogroup))
            }
        }

        if (!shiny_plot & interactive) gp <- ggiraph::girafe(ggobj=gp)



        # warning("interactive plotting is currently in beta")
        # vx <- gr_v$x[gr_v$v_ind]
        # vy <- gr_v$y[gr_v$v_ind]
        # gp <- plotly::plot_ly(x=~vx, y=~vy, mode="markers", opacity=1,
        #                       color=gr_v$colour[gr_v$v_ind],
        #                       colorscale="YlGnBu",
        #                       size=gr_v$size[gr_v$v_ind],
        #                       text=gr_v$label_long[gr_v$v_ind],
        #                       hoverinfo="text", ...)
        #
        # # gr_e$colour_new <- gr_e$colour
        # # gr_e$colour_new[gr_e$colour<0] <- "turquoise"
        # # gr_e$colour_new[gr_e$colour>0] <- "cherry"
        # # gr_e$colour_new[gr_e$colour==0] <- "grey50"
        # if (!show_bgedges) {
        #     el <- purrr::map(which(gr_e$e_ind), function(j) {
        #         list(type="line", line=list(color=gr_e$colour[j], width=.5),
        #              x0=gr_e$from.x[j], x1=gr_e$to.x[j],
        #              y0=gr_e$from.y[j], y1=gr_e$to.y[j])
        #     })
        # } else {
        #     el <- purrr::map(seq_len(nrow(gr_e)), function(j) {
        #         list(type="line", line=list(color=gr_e$colour[j], width=.5),
        #              x0=gr_e$from.x[j], x1=gr_e$to.x[j],
        #              y0=gr_e$from.y[j], y1=gr_e$to.y[j])
        #     })
        # }
        #
        # axis <- list(title="", showgrid=FALSE, showticklabels=FALSE,
        #              zeroline=FALSE)
        #
        # gp <- plotly::layout(gp, title=main, shapes=el, xaxis=axis, yaxis=axis)
    }
    else { # visNet
        v_ind <- gr$v$v_ind
        v_label <- purrr::map_chr(
            stringr::str_split(gr$v$label[v_ind]," "),
            function(x) base::paste0(x[1:2], collapse=" "))
        nodes <- base::data.frame(
            id=gr$v$phenotype[v_ind],
            # get only phenotype and p
            label=v_label,
            title=gr$v$label_long[v_ind],
            color=flowGraph:::noTOcol(gr$v$colour[v_ind],
                                      colourp=colour_palette),
            size=gr$v$size[v_ind],
            group=gr$v$phenogroup[v_ind],
            level=gr$v$phenolayer[v_ind],
            stringsAsFactors=FALSE)

        e_ind <- gr$e$from%in%gr$v$phenotype[v_ind] &
            gr$e$to%in%gr$v$phenotype[v_ind]
        if (base::sum(e_ind)==0) {
            edges <- base::data.frame(from=gr$e$from[1], to=gr$e$to[1],
                                      hidden=TRUE)
        } else {
            edges <- base::data.frame(
                from=gr$e$from[e_ind],
                to=gr$e$to[e_ind],
                width=gr$e$size[e_ind],
                label=gr$e$marker[e_ind],
                stringsAsFactors=FALSE)
            if (!is.null(gr$e$size1)) {
                edges$title=base::paste0(
                    base::signif(gr$e$size1,3), " vs ",
                    base::signif(gr$e$size2,3))[e_ind] # tooltip
            }
        }

        gp <- visNetwork::visNetwork(nodes, edges) %>%
            visNetwork::visNodes(shape="dot") %>%
            visNetwork::visEdges(arrows="to") %>%
            visNetwork::visOptions(highlightNearest=TRUE, nodesIdSelection=TRUE) %>%
            visNetwork::visInteraction(hover=TRUE, multiselect=!shiny_plot) %>%
            # visIgraphLayout(layout="layout_as_tree") %>% # "layout_in_circle"
            visNetwork::visHierarchicalLayout(
                parentCentralization=TRUE,
                treeSpacing=50,
                nodeSpacing=300,
                levelSeparation=200
            ) %>%
            visNetwork::visEvents(select="function(nodes) {
                                  Shiny.onInputChange('current_node_selection', nodes.id);
                                  ;}")

}
    # plot(gp)
    return(gp)
}

# colour palette; cols = a colour for each node
noTOcol <- function(
    values, colourp=c('blue','cyan','yellow','red'), col_names=NULL,
    colour_n=1000) {
    colorFunc <- grDevices::colorRampPalette(colourp)
    z <- base::pretty(c(min(values), max(values))) # colour breakpoints
    colinds <- base::unlist(
        base::lapply(seq_len(base::length(values)), function(i)
            base::max(1, base::ceiling(
                (values-min(z))*1000/(max(z)-min(z)))[i]) ))
    cols=colorFunc(colour_n)[colinds]
    if (base::is.null(col_names) & !base::is.null(base::names(values)))
        col_names <- base::names(values)
    base::names(cols) <- col_names
    return(cols)
}


#' @title Creates a QQ plot of a summary statistic.
#' @description Creates a QQ plot of a summary statistic.
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
#'  \code{feature} (feature name), \code{test_name} (summary statistic name),
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
#' @param logged A logical indicating whether or not to log the  summary
#'  statistic p value.
#' @param p_thres A double indicating a summary statistic threshold
#'  e.g. if we are plotting a T test summary statistic, we can set the threshold
#'  to .05; nodes with a p-value greater than .05 will not be plotted.
#' @param filter_adjust0 A numeric variable indicating what percentage of
#'  SpecEnr values compared (minimum) should be not close to 0.
#'  Set to 1 to not conduct filtering.
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
#' @param shiny_plot A logical indicating whether this plot is made for shiny;
#'  users don't need to change this.
#' @param main A string or the title of the plot; if left as \code{NULL},
#'  a default title will be applied.
#' @param interactive A logical indicating whether or not plot should be an
#'  interactive ggiraph plot as opposed to a static plot.
#' @param path A string indicating the path to where the function should save
#'  the plot; leave as \code{NULL} to not save the plot. Static plots are saved
#'  as PNG, interactive plots are saved as HTML.
#' @return A static or interactive qq plot.
#' @details The interactive plot is made using the \code{ggiraph} package.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos15)
#'  fg <- flowGraph(fg_data_pos15$count, class=fg_data_pos15$meta$class,
#'                  no_cores=no_cores)
#'
#'  fg_plot_qq(fg, type="node", summary_meta=NULL, adjust_custom="byLayer", index=1,
#'          interactive=TRUE, logged=FALSE)
#'
#'  fg_plot_qq(fg, type="node", summary_meta=NULL, adjust_custom="byLayer", index=1,
#'          interactive=FALSE, logged=FALSE)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_plot}}
#'  \code{\link[flowGraph]{plot_gr}}
#'  \code{\link[flowGraph]{fg_get_feature}}
#'  \code{\link[flowGraph]{fg_get_summary}}
#' @rdname fg_plot_qq
#' @export
#' @importFrom htmlwidgets saveWidget
#' @importFrom ggiraph girafe geom_point_interactive
#' @importFrom ggplot2 ggplot ggtitle geom_point geom_abline labs ggsave
fg_plot_qq <- function(
    fg, type="node", index=1, summary_meta=NULL, adjust_custom="byLayer",
    logged=TRUE, p_thres=.05,
    filter_adjust0=1, filter_es=0,
    filter_btwn_tpthres=1, filter_btwn_es=0,
    shiny_plot=FALSE,
    main=NULL, interactive=FALSE, path=NULL
) {
    type <- base::match.arg(type, c("node", "edge"))

    if (shiny_plot) interactive <- TRUE

    index <- flowGraph:::fg_get_summary_index(
        fg,type=type, index, summary_meta)
    summary_meta <- base::unlist(fg@summary_desc[[type]][index,])
    if (!base::grepl("SpecEnr",base::unlist(summary_meta["feat"])))
        filter_adjust0 <- 1
    qvals_ <- fg_get_summary(
        fg, type, index, summary_meta, default_p_thres=p_thres,
        filter_adjust0=filter_adjust0, filter_es=filter_es,
        filter_btwn_tpthres=filter_btwn_tpthres, filter_btwn_es=filter_btwn_es,
        adjust_custom=adjust_custom
    )
    qvals <- qvals_$values

    if (base::is.null(main))
        main <- base::paste0(
            "(", sum(qvals<p_thres), "/", base::length(qvals), ") ",
            base::ifelse(logged,"logged ",""), "qq plot.\n",
            "- ", type, " feature: ", summary_meta["feat"], ".\n",
            "- p-value: ", summary_meta["test_name"],
            "\n- comparing ", summary_meta["class"], ", labels: ", summary_meta["label1"], " & ", summary_meta["label2"], ".")

    qo <- base::order(qvals)
    # m1 <- fg_get_feature_means(fg, type=type, feature="prop", id=qvals_$id1)

    uni <- seq_len(base::length(qvals))/(base::length(qvals)+1)

    p_thres_ <- p_thres
    if (logged) {
        qvals <- log(qvals)
        if (any(qvals==-Inf)) {
            warning("converting log(0) p values to minimum value for display")
            qvals[qvals==-Inf] <- min( qvals[is.finite(qvals)])
        }
        uni <- log(uni)
        p_thres <- log(p_thres)
    }



    df <- data.frame(y=qvals, cohensd_size=qvals_$cohensd_size,
                     phenotype=base::names(qvals),
                     phenogroup=fg@graph$v$phenogroup,
                     d_size=colMeans(as.matrix(fg@feat$node$count)))
    df <- df[qo,]
    df$x <- uni
    qp <- ggplot2::ggplot(df, ggplot2::aes(
        y=y, x=x, colour=cohensd_size, alpha=.3, stroke=1)) +
        ggplot2::geom_abline(intercept=0, slope=1) +
        ggplot2::ggtitle(main) +
        ggplot2::labs(x="uniform distribution", y="p-value",
                      col="cohen's d", size="count mean") +
        ggplot2::geom_hline(yintercept=p_thres)

    if (interactive) {
        if (shiny_plot) {
            qp <- qp + ggiraph::geom_point_interactive(alpha=.4,
                ggplot2::aes(tooltip=phenotype, data_id=phenotype, size=d_size))
        } else {
            qp <- qp + ggiraph::geom_point_interactive(shape=1,
                ggplot2::aes(tooltip=phenotype, data_id=phenogroup, size=d_size))
        }
        if (!shiny_plot) {
            qp <- ggiraph::girafe(ggobj=qp)
            if (!base::is.null(path))
                htmlwidgets::saveWidget(
                    qp, ifelse(base::grepl("[.]html$",path, ignore.case=TRUE),
                               path, base::paste0(path, ".html")))
        }
    } else {
        qp <- qp + ggplot2::geom_point(shape=1, ggplot2::aes(size=d_size))
        if (!base::is.null(path))
            suppressMessages({
                ggplot2::ggsave(
                    ifelse(base::grepl("[.]png$",path, ignore.case=TRUE),
                           path, base::paste0(path, ".png")),
                    plot=qp, scale=1, width=5, height=5,
                    units="in", dpi=500, limitsize=TRUE)
            })
    }
    return(qp)
}


#' @title Creates a boxplot of the values of one node/edge
#' @description Creates a boxplot comparing the
#'  features of samples belonging to different classes corresponding
#'  to an existing summary statistic using ggplot2.
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
#'  \code{feature} (feature name), \code{test_name} (summary statistic name),
#'  \code{class} (class), \code{label1}, and \code{label2} (class labels compared).
#'  See \code{\link[flowGraph]{fg_get_summary_desc}} for details.
#' @param node_edge An integer/index of or a string of the cell population (node) /
#'  edge name (edge) the user wants to plot.
#' @param adjust_custom A function or a string indicating the
#' test adjustment method to use.
#'  If a string is provided, it should be one of
#'  \code{c("holm", "hochberg", "hommel",
#'  "bonferroni", "BH", "BY", "fdr", "none")} (see \code{p.adjust.methods}).
#'  If a function is provided, it should take as input
#'  a numeric vector and output the
#'  same vector adjusted.
#' @param p_thres A numeric variable indicating a p-value threshold
#' @param filter_adjust0 A numeric variable indicating what percentage of
#'  SpecEnr values compared (minimum) should be not close to 0.
#'  Set to 1 to not conduct filtering.
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
#' @param paired A logical indicating whether the summary is paired.
#' @param dotplot A logical indicating whether or not to plot sample points.
#' @param outlier A logical indicating whether or not outliers should be plotted.
#' @param all_labels A logical indicating whether or not to plot samples of all
#'  classes outside of just those used in the summary statistic test.
#' @param main A string or the title of the plot; if left as \code{NULL},
#'  a default title will be applied.
#' @param path A string indicating the path to where the function should save
#'  the plot; leave as \code{NULL} to not save the plot. Static plots are saved
#'  as PNG.
#' @return A static boxplot.
#' @details The plot is made using the \code{ggplot2} package. The interactive
#'  version is the same as the static version, it is only here to support the
#'  shiny app.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos15)
#'  fg <- flowGraph(fg_data_pos15$count, class=fg_data_pos15$meta$class,
#'                  no_cores=no_cores)
#'
#'  fg_plot_box(fg, type="node", summary_meta=NULL, adjust_custom="byLayer", index=1, node_edge=10)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_plot}}
#'  \code{\link[flowGraph]{plot_gr}}
#'  \code{\link[flowGraph]{fg_get_feature}}
#'  \code{\link[flowGraph]{fg_get_summary}}
#'  \code{\link[flowGraph]{fg_plot_qq}}
#' @rdname fg_plot_box
#' @export
#' @importFrom htmlwidgets saveWidget
#' @importFrom ggplot2 ggsave ggplot ggtitle geom_point geom_dotplot stat_summary geom_line aes position_dodge labs
#' @importFrom ggiraph girafe geom_boxplot_interactive
#' @importFrom stringr str_split
#' @importFrom grDevices boxplot.stats
#' @importFrom utils head tail
fg_plot_box <- function(
    fg, type="node", index=1, summary_meta=NULL, node_edge=1,
    adjust_custom="byLayer", p_thres=.05,
    filter_adjust0=.5, filter_es=.5,
    filter_btwn_tpthres=.05, filter_btwn_es=.5,
    paired=FALSE, dotplot=TRUE, outlier=TRUE, all_labels=FALSE,
    main=NULL, path=NULL
) {
    type <- base::match.arg(type, c("node", "edge"))

    index <- fg_get_summary_index(fg,type="node",index,summary_meta)
    summary_meta <- base::unlist(fg@summary_desc[[type]][index,])
    pp  <- fg_get_summary(
        fg, type, index, summary_meta, default_p_thres=p_thres,
        filter_adjust0=filter_adjust0, filter_es=filter_es,
        filter_btwn_tpthres=filter_btwn_tpthres, filter_btwn_es=filter_btwn_es,
        adjust_custom=adjust_custom
    )
    p = pp$values
    node_edge <-
        ifelse(base::is.character(node_edge),
               base::which(base::names(p)==node_edge), node_edge)

    feature <- summary_meta["feat"]
    features <- se_feats(feature)
    if (all_labels) {
        id_ord <- base::order(fg@meta[,summary_meta["class"]])
        class_ <- fg@meta[id_ord,summary_meta["class"]]
        dta <- data.frame(
            val=fg@feat[[type]][[feature]][id_ord,node_edge],
            class=class_,
            ID=base::unlist(purrr::map(unique(class_), function(x)
                base::seq_len(sum(class_==x))))) # paired assum ordered!

    } else {
        a = fg@feat[[type]][[feature]][pp$id1,node_edge]
        b = fg@feat[[type]][[feature]][pp$id2,node_edge]

        classes <- c(summary_meta["label1"], summary_meta["label2"])
        class_ <- base::append(base::rep(classes[1], base::length(a)),
                               base::rep(classes[2], base::length(b)))
        dta <- data.frame(
            val=base::append(a,b),
            class=class_,
            ID=base::append(base::seq_len(base::length(a)),
                            base::seq_len(base::length(b)))) # paired
    }

    if (base::is.null(main))
        main <- base::paste0(
            "boxplot (",
            ifelse(base::is.numeric(node_edge),
                   base::names(p)[node_edge], node_edge),").\n",
            "- ", type, " feature: ", feature ," (",summary_meta["feat"],").\n",
            "- p-value (p=", base::round(p[node_edge],3),"): ",
            summary_meta["test_name"], "\n- comparing ",
            summary_meta["class"], ", labels: ",
            "\n    ", summary_meta["label1"], " (mean=",
            round(pp$m1[node_edge],3),") & ", summary_meta["label2"],
            " (mean=",round(pp$m2[node_edge],3),").")

    # plot
    gp <- ggplot2::ggplot(
        dta, ggplot2::aes(x=class, y=val, fill=class)) +
        ggplot2::labs(y=base::paste0(
            type, " feature values (",ifelse(
                base::is.numeric(node_edge),
                base::names(p)[node_edge], node_edge),")"))

    gp <- gp + ggplot2::geom_boxplot(outlier.shape=ifelse(outlier,19,NA))

    if (paired)
        gp <- gp + ggplot2::geom_line(
            ggplot2::aes(group=ID), colour="black", linetype="11")

    if (dotplot)
        gp <- gp + ggplot2::geom_dotplot(
            binaxis='y', stackdir='center',
            position=ggplot2::position_dodge(1), dotsize=.5)

    gp <- gp +
        ggplot2::stat_summary(fun=mean, geom="point", shape=23, size=3) +
        ggplot2::ggtitle(main)

    if (!outlier) {
        bstats <- purrr::map(base::unique(class_), function(x)
            grDevices::boxplot.stats(dta$val[class_==x])$stats)
        ymax <- max(purrr::map_dbl(bstats, utils::tail, 1))
        ymin <- min(purrr::map_dbl(bstats, utils::head, 1))

        gp <- gp + ggplot2::ylim(ymin, ymax)
    }

    if (!base::is.null(path)) {
        base::suppressMessages({
            ggplot2::ggsave(
                ifelse(base::grepl("[.]png$",path, ignore.case=TRUE),
                       path, base::paste0(path, ".png")),
                plot=gp, scale=1, width=5, height=5,
                units="in", dpi=500, limitsize=TRUE)
            })
    }
    gp
}


# plot a set of boxplots for SpecEnr
#' @importFrom gridExtra arrangeGrob
#' @importFrom ggplot2 element_blank coord_cartesian
fg_plot_box_set <- function(
    fg, type="node", index=1, summary_meta=NULL,
    adjust_custom="byLayer", node_edge=1,
    filter_adjust0=.5, filter_es=.5,
    filter_btwn_tpthres=.05, filter_btwn_es=.5,
    paired=FALSE, dotplot=TRUE, outlier=TRUE,
    main=NULL, path=NULL
) {
    type <- match.arg(type, c("node", "edge"))

    index <- fg_get_summary_index(fg,type="node",index,summary_meta)
    summary_meta <- base::unlist(fg@summary_desc[[type]][index,])
    feature <- summary_meta["feat"]
    if (!base::grepl("SpecEnr",feature)) {
        fg_plot_box(
            fg, type=type, index=index, summary_meta=summary_meta,
            adjust_custom=adjust_custom, node_edge=node_edge,
            paired=paired, dotplot=dotplot, outlier=outlier, all_labels=FALSE,
            main=main, path=path,
            filter_adjust0=filter_adjust0, filter_es=filter_es,
            filter_btwn_tpthres=filter_btwn_tpthres,
            filter_btwn_es=filter_btwn_es
        )
    }
    node_edge <- ifelse(is.character(node_edge),
                        base::which(fg@graph$v$phenotype==node_edge), node_edge)
    if (base::length(node_edge)==0) {
        warning("node/edge not found")
        return(NULL)
    }
    pp <- fg_get_summary(
        fg, type="node", index, summary_meta, adjust_custom=adjust_custom,
        filter_adjust0=filter_adjust0, filter_es=filter_es,
        filter_btwn_tpthres=filter_btwn_tpthres, filter_btwn_es=filter_btwn_es)
    p <- pp$values
    dfb <- pp$btwn
    base::rownames(dfb) <- dfb$phenotype

    features <- flowGraph:::se_feats(feature)
        a1 <- fg@feat[[type]][[features[2]]][pp$id1,node_edge]
        a2 <- fg@feat[[type]][[features[2]]][pp$id2,node_edge]
        b1 <- fg@feat[[type]][[features[3]]][pp$id1,node_edge]
        b2 <- fg@feat[[type]][[features[3]]][pp$id2,node_edge]

    a = fg@feat[[type]][[feature]][pp$id1,node_edge]
    b = fg@feat[[type]][[feature]][pp$id2,node_edge]
    al <- base::length(a)
    bl <- base::length(b)

    classes <- c(summary_meta["label1"], summary_meta["label2"])
    dta <- data.frame(
        val=c(a,b, a1,a2, b1,b2),
        feature=base::rep(features, each=al+bl),
        class=base::rep(append(rep(classes[1],length(a)),
                               rep(classes[2],length(b))),3),
        ID=base::rep(base::append(base::seq_len(base::length(a)),
                                  base::seq_len(base::length(b))),3) # paired
    )


    if (base::is.null(main))
        main <- base::paste0(
            "boxplot (",
            ifelse(base::is.numeric(node_edge),
                   base::names(p)[node_edge], node_edge),").\n",
            "- ", type, " feature: ", feature ," (",summary_meta["feat"],").\n",
            "- p-value (p=", base::round(p[node_edge],3),"): ",
            summary_meta["test_name"], "\n- comparing ",
            summary_meta["class"], ", labels: ",
            "\n    ", summary_meta["label1"], " (mean=",
            base::round(pp$m1[node_edge],3),") & ",
            summary_meta["label2"], " (mean=",
            base::round(pp$m2[node_edge],3),").")

    # specenr/actual/expect class 1 vs class 2
    gp <- ggplot2::ggplot(dta[dta$feature==feature,],
                            ggplot2::aes(x=class, y=val, fill=class)) +
        ggplot2::labs(y=base::paste0(type, " SpecEnr feature value")) +
        ggplot2::geom_boxplot(outlier.shape=ifelse(outlier,19,NA)) +
        ggplot2::stat_summary(fun=mean, geom="point", shape=23, size=3) +
        ggplot2::ggtitle(base::paste0(
            "(",classes[2],"-",base::length(a)," vs ",
            classes[2],"-",base::length(b),")\n",main))

    gpae <- ggplot2::ggplot(dta[dta$feature!=features[1],],
                                ggplot2::aes(x=class, y=val, fill=class)) +
        ggplot2::labs(y=base::paste0(type, " raw feature values")) +
        ggplot2::ggtitle(base::paste0("feature: ", features[2])) +
        ggplot2::geom_boxplot(outlier.shape=ifelse(outlier,19,NA)) +
        ggplot2::stat_summary(fun=mean, geom="point", shape=23, size=3) +
        ggplot2::theme(legend.position="none") +
        ggplot2::facet_grid(cols=ggplot2::vars(feature))

    # gpe <- ggplot2::ggplot(dta[dta$feature==features[3],],
    #                         ggplot2::aes(x=class, y=val, fill=class)) +
    #     ggplot2::ggtitle(paste0("feature: ", features[3])) +
    #     ggplot2::geom_boxplot(outlier.shape=ifelse(outlier,19,NA)) +
    #     ggplot2::stat_summary(fun=mean, geom="point", shape=23, size=3) +
    #     ggplot2::theme(legend.position="none",
    #                    axis.title.y=ggplot2::element_blank())

    # class 1/2 actual vs expected
    gp12 <- ggplot2::ggplot(dta[dta$feature!=feature,],
                           ggplot2::aes(x=feature, y=val))+
        ggplot2::labs(y=base::paste0(type, " raw feature values")) +
        ggplot2::ggtitle(base::paste0(
            "diff btwn diff (p=",
            base::signif(dfb$btp[node_edge],3),
            ", CohenD=",signif(dfb$bcd[node_edge],3),")")) +
        ggplot2::geom_boxplot(outlier.shape=ifelse(outlier,19,NA)) +
        ggplot2::stat_summary(fun=mean, geom="point", shape=23, size=3) +
        ggplot2::theme(legend.position="none") +
        ggplot2::facet_grid(cols=ggplot2::vars(class))

    # gp2 <- ggplot2::ggplot(dta[dta$feature!=feature & dta$class==classes[2],],
    #                        ggplot2::aes(x=feature, y=val))+
    #     ggplot2::ggtitle(paste0("class label: ", classes[2],
    #                             " (p=",signif(dfb$tpv2[node_edge],3),")")) +
    #     ggplot2::geom_boxplot(outlier.shape=ifelse(outlier,19,NA)) +
    #     ggplot2::stat_summary(fun=mean, geom="point", shape=23, size=3) +
    #     ggplot2::ggtitle(paste0("class label: ", classes[2])) +
    #     ggplot2::theme(legend.position="none",
    #                    axis.title.y=ggplot2::element_blank())

    if (paired) {
        gp <- gp + ggplot2::geom_line(
            ggplot2::aes(group=ID), colour="black", linetype="11")
        gpae <- gpae + ggplot2::geom_line(
            ggplot2::aes(group=ID), colour="black", linetype="11")
    }

    if (dotplot) {
        gp <- gp + ggplot2::geom_dotplot(
            binaxis='y', stackdir='center',
            position=ggplot2::position_dodge(1), dotsize=.3)
        gpae <- gpae + ggplot2::geom_dotplot(
            binaxis='y', stackdir='center',
            position=ggplot2::position_dodge(1), dotsize=.3)
        # gpe <- gpe + ggplot2::geom_dotplot(
        #     binaxis='y', stackdir='center',
        #     position=ggplot2::position_dodge(1), dotsize=.3)
        gp12 <- gp12 + ggplot2::geom_dotplot(
            binaxis='y', stackdir='center',
            position=ggplot2::position_dodge(1), dotsize=.3)
        # gp2 <- gp2 + ggplot2::geom_dotplot(
        #     binaxis='y', stackdir='center',
        #     position=ggplot2::position_dodge(1), dotsize=.3)
    }

    if (!outlier) {
        ylim1 = grDevices::boxplot.stats(a)$stats[c(1,5)]
        ylim2 = grDevices::boxplot.stats(b)$stats[c(1,5)]
        gp <- gp + ggplot2::coord_cartesian(
            ylim=c(min(ylim1[1],ylim2[1]), max(ylim1[2],ylim2[2])))
        ylim1 = grDevices::boxplot.stats(a1)$stats[c(1,5)]
        ylim2 = grDevices::boxplot.stats(a2)$stats[c(1,5)]
        ylim1_ = grDevices::boxplot.stats(b1)$stats[c(1,5)]
        ylim2_ = grDevices::boxplot.stats(b2)$stats[c(1,5)]
        gpae <- gpae + ggplot2::coord_cartesian(
            ylim=c(min(ylim1[1],ylim2[1],ylim1_[1],ylim2_[1]),
                   max(ylim1[2],ylim2[2],ylim1_[2],ylim2_[2])))
        # ylim1 = grDevices::boxplot.stats(b1)$stats[c(1,5)]
        # ylim2 = grDevices::boxplot.stats(b2)$stats[c(1,5)]
        # gpe <- gpe + ggplot2::coord_cartesian(
        #     ylim=c(min(ylim1[1],ylim2[1]), max(ylim1[2],ylim2[2])))
        ylim1 = grDevices::boxplot.stats(a1)$stats[c(1,5)]
        ylim2 = grDevices::boxplot.stats(b1)$stats[c(1,5)]
        ylim1_ = grDevices::boxplot.stats(a2)$stats[c(1,5)]
        ylim2_ = grDevices::boxplot.stats(b2)$stats[c(1,5)]
        gp12 <- gp12 + ggplot2::coord_cartesian(
            ylim=c(min(ylim1[1],ylim2[1],ylim1_[1],ylim2_[1]),
                   max(ylim1[2],ylim2[2],ylim1_[2],ylim2_[2])))
        # ylim1 = grDevices::boxplot.stats(a2)$stats[c(1,5)]
        # ylim2 = grDevices::boxplot.stats(b2)$stats[c(1,5)]
        # gp2 <- gp2 + ggplot2::coord_cartesian(
        #     ylim=c(min(ylim1[1],ylim2[1]), max(ylim1[2],ylim2[2])))
    }

    gp_grid <- gridExtra::arrangeGrob(
        gp, gpae, gp12,# gp1, gp2,
        layout_matrix=cbind(c(1,1,2,3),c(1,1,2,3)))

    # # histogram
    # gp <- ggplot2::ggplot() + ggplot2::ggtitle(main) +
    #     ggplot2::geom_density(ggplot2::aes(x=x), colour="red",
    #                           data=data.frame(x=a)) +
    #     ggplot2::geom_density(ggplot2::aes(x=x), colour="blue",
    #                           data=data.frame(x=b))

    # if (interactive & !shiny_plot)
    #     gp <- ggiraph::girafe(ggobj=gp)

    if (!base::is.null(path)) {
        # if (interactive & !shiny_plot) {
        #     htmlwidgets::saveWidget(
        #         gp, ifelse(grepl("[.]html$",path, ignore.case=TRUE),
        #                    path, paste0(path, ".html")))
        # } else if (!interactive) {
        suppressMessages({
            ggplot2::ggsave(ifelse(base::grepl(
                "[.]png$",path, ignore.case=TRUE),
                path, base::paste0(path, ".png")),
                plot=gp_grid, scale=1, width=5, height=10,
                units="in", dpi=500, limitsize=TRUE)
        })
        # }
    }
    gp_grid
}



#' @title Creates a p value vs feature difference plot
#' @description Creates a p value vs feature difference plot where the
#'  difference is that of the
#'  features of samples belonging to different classes corresponding
#'  to an existing summary statistic.
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
#'  \code{feature} (feature name), \code{test_name} (summary statistic name),
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
#' @param logged A logical indicating whether or not to log the summary
#'  statistic p value.
#' @param label_max An integer indicating the maximum number of max difference
#'  and/or min p value nodes/edges that should be labelled.
#' @param p_thres A numeric variable indicating a p-value threshold; a line will
#'  be plotted at this threshold.
#' @param filter_adjust0 A numeric variable indicating what percentage of
#'  SpecEnr values compared (minimum) should be not close to 0.
#'  Set to 1 to not conduct filtering.
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
#' @param shiny_plot A logical indicating whether this plot is made for shiny;
#'  users don't need to change this.
#' @param nodes_max An integer indicating maximum number of nodes to plot;
#'  this limit is set for interactive plots only.
#' @param main A string or the title of the plot; if left as \code{NULL},
#'  a default title will be applied.
#' @param interactive A logical variable indicating whether the plot should be
#'  an interactive plot; see package \code{ggiraph}.
#' @param path A string indicating the path to where the function should save
#'  the plot; leave as \code{NULL} to not save the plot. Static plots are saved
#'  as PNG.
#' @return A static or interactive p value vs difference plot.
#' @details The interactive plot is made using the \code{ggiraph} package.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos15)
#'  fg <- flowGraph(fg_data_pos15$count, class=fg_data_pos15$meta$class,
#'                  no_cores=no_cores)
#'
#'  gp <- fg_plot_pVSdiff(fg, type="node", summary_meta=NULL,
#'                        adjust_custom="byLayer", index=1, label_max=10)
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_plot}}
#'  \code{\link[flowGraph]{plot_gr}}
#'  \code{\link[flowGraph]{fg_get_feature}}
#'  \code{\link[flowGraph]{fg_get_summary}}
#'  \code{\link[flowGraph]{fg_plot_qq}}
#' @rdname fg_plot_pVSdiff
#' @export
#' @importFrom ggplot2 ggsave ggplot ggtitle geom_point labs geom_vline
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggiraph girafe geom_point_interactive
#' @importFrom htmlwidgets saveWidget
#' @importFrom utils head tail
fg_plot_pVSdiff <- function(
    fg, type="node", index=1, summary_meta=NULL, adjust_custom="byLayer",
    logged=TRUE, label_max=5, p_thres=.05,
    filter_adjust0=1, filter_es=0,
    filter_btwn_tpthres=1, filter_btwn_es=0,
    shiny_plot=FALSE, nodes_max=30,
    main=NULL, interactive=FALSE, path=NULL
) {
    type <- match.arg(type, c("node", "edge"))

    if (shiny_plot) interactive <- TRUE

    index <- flowGraph:::fg_get_summary_index(fg, type, index, summary_meta)
    summary_meta <- base::unlist(fg@summary_desc[[type]][index,])
    if (!base::grepl("SpecEnr",base::unlist(summary_meta["feat"])))
        filter_adjust0 <- 1
    pp <- fg_get_summary(
        fg, type="node", index, summary_meta, default_p_thres=p_thres,
        filter_adjust0=filter_adjust0, filter_es=filter_es,
        filter_btwn_tpthres=filter_btwn_tpthres, filter_btwn_es=filter_btwn_es,
        adjust_custom=adjust_custom
    )
    p <- pp$values
    mse_c = pp$m1
    mse_m = pp$m2
    mse_ = mse_m - mse_c

    p_ <- p
    if (logged) {
        p_ = -log(p)
        p_[is.infinite(p_)] = max(p_[!is.infinite(p_)])
    }
    dta <- data.frame(log_p_value=p_, difference=mse_,
                      cohensd_size=pp$cohensd_size)
    dta <- cbind(dta, fg@graph$v)
    dta$phenotype_ <- base::names(p)
    dta$avg_count <- base::colMeans(base::as.matrix(fg@feat$node$count))

    if (base::is.null(main))
        main <- base::paste0(
            "(",base::sum(p<p_thres),"/",base::length(p),") ",
            "-ln(p-value) vs feature difference plot.\n",
            "- ", type, " feature: ", summary_meta["feat"], ".\n",
            "- p-value: ", summary_meta["test_name"],
            "\n- comparing ",
            summary_meta["class"], ", labels: ", summary_meta["label1"], " & ",
            summary_meta["label2"], ".")


    if (!interactive) {
        label_ind <- p<p_thres &
            (base::names(p) %in%
                 base::names(utils::head(base::sort(p),label_max)) |
            base::names(p) %in%
                 base::names(utils::tail(
                     base::sort(base::abs(mse_)),label_max)))

        dta$phenotype <- ifelse(label_ind, base::names(p), "")

        gp <- ggplot2::ggplot(
            dta, ggplot2::aes(x=log_p_value, y=difference,
                     color=cohensd_size)) +
            ggplot2::geom_point(shape=1, ggplot2::aes(size=avg_count)) +
            ggplot2::ggtitle(main) +
            ggplot2::labs(x=ifelse(logged,"-ln(p-value)","p-value"),
                            y="difference between mean feature values",
                            col="cohen's d", size="mean count") +
            ggplot2::geom_vline(xintercept=ifelse(logged,-log(p_thres),p_thres))
        if (any(label_ind))
            gp <- gp + ggrepel::geom_text_repel(
                ggplot2::aes(label=phenotype),
                nudge_x      = -0.35,
                direction    = "y",
                hjust        = 1,
                segment.size = 0.2
            )
        if (!base::is.null(path))
            suppressMessages({
                ggplot2::ggsave(ifelse(base::grepl(
                    "[.]png$",path, ignore.case=TRUE),
                    path, base::paste0(path, ".png")), plot=gp)
            })
    } else {
        # i was trying out the speed of different interactive plots here, don't mind teh comments.
        dta$phenotype <- base::names(p)
        dta$phenogroup <- fg@graph$v$phenogroup
        if (shiny_plot & nodes_max!=Inf) {
            node_ind <-
                (base::names(p) %in%
                     base::names(utils::head(base::sort(p),nodes_max)) |
                     base::names(p) %in%
                     base::names(utils::tail(base::sort(abs(mse_)),nodes_max)))
            if (sum(node_ind)==0) return(NULL)
            dta <- dta[node_ind,]
        }

        # if (gg) {
        start1 <- Sys.time()
        gp <- ggplot2::ggplot(
            dta, ggplot2::aes(x=log_p_value,y=difference,color=cohensd_size)) +
            ggplot2::ggtitle(main) +
            ggplot2::labs(x="-ln(p-value)",
                          y="difference between mean feature values",
                          col="cohen's d", size="mean count") +
            ggplot2::geom_vline(xintercept=-log(p_thres))
        if (shiny_plot) {
            gp <- gp + ggiraph::geom_point_interactive(
                ggplot2::aes(tooltip=phenotype, data_id=phenotype,
                                       size=avg_count), alpha=.3)
        } else {
            gp <- gp + ggiraph::geom_point_interactive(
                ggplot2::aes(tooltip=phenotype, data_id=phenogroup,
                             size=avg_count))
            gp <- ggiraph::girafe(ggobj=gp)
        }
        gp
    #     flowGraph:::time_output(start1)
    # } else {
    #     start1 <- Sys.time()
    #     gp <- highcharter::hchart(
    #         dta, "scatter",
    #         highcharter::hcaes(x=log_p_value, y=difference,
    #                            group=phenogroup, size=avg_count,
    #                            color=cohensd_size))
    #     gp
    #     flowGraph:::time_output(start1)
    # }
    # start1 <- Sys.time()
    # gp <- plotly::plot_ly(
    #     dta, x = ~log_p_value, y = ~difference,
    #     # Hover text:
    #     text = ~phenotype_,
    #     color = ~cohensd_size, size = ~avg_count, group=~fg@graph$v$phenogroup
    # )
    # gp
    # flowGraph:::time_output(start1)

        if (!base::is.null(path))
            htmlwidgets::saveWidget(
                gp, ifelse(base::grepl("[.]html$", path, ignore.case=TRUE),
                       path, base::paste0(path, ".html")))
    }
    return(gp)

}





#' @title Saves numerous plots for all summary statistics to a folder.
#' @description Saves numerous plots for all summary statistics in a given
#'  flowGraph object to a user specified folder.
#' @param fg flowGraph object.
#' @param plot_path A string indicating the folder path to where the function
#'  should save the plots.
#' @param plot_types A string or a vector of strings indicating what feature
#'  types and their summaries the function should plot for: 'node' or 'edge'.
#' @param interactive A logical indicating whether the QQ plot, p-value vs
#'  difference plot, and the
#'  cell hierarchy plots should be interactive; see functions
#'  \code{fg_plot} and \code{fg_plot_qq}.
#' @param adjust_custom A function or a string indicating the
#' test adjustment method to use.
#'  If a string is provided, it should be one of
#'  \code{c("holm", "hochberg", "hommel",
#'  "bonferroni", "BH", "BY", "fdr", "none")} (see \code{p.adjust.methods}).
#'  If a function is provided, it should take as input
#'  a numeric vector and output the
#'  same vector adjusted.
#' @param label_max An integer indicating how many labels should be shown
#'  in the functions \code{fg_plot_pVSdiff} and \code{fg_plot}.
#' @param box_no An integer indicating the maximum number of boxplots to save;
#'  used in function \code{fg_plot_box}.
#' @param paired A logical indicating whether the summary is paired; used in
#'  function \code{fg_plot_box}.
#' @param logged A logical indicating whether or not to log the summary
#'  statistic p value in the qq plots.
#' @param filter_adjust0 A numeric variable indicating what percentage of
#'  SpecEnr values compared (minimum) should be not close to 0.
#'  Set to 1 to not conduct filtering. This parameter is used for the QQ and the
#'  pVSdifference plots.
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
#' @param overwrite A logical variable indicating whether or not to replace
#'  old plots if they exist under the same folder name.
#' @param ... Other parameters for the \code{fg_plot} function.
#' @return No return; plots are saved to file.
#' @details The interactive plots are made using the \code{ggiraph} package.
#' @examples
#'
#'  no_cores <- 1
#'  data(fg_data_pos15)
#'  fg <- flowGraph(fg_data_pos15$count, class=fg_data_pos15$meta$class,
#'                  no_cores=no_cores)
#'
#'  fg_save_plots(fg, "temp")
#'
#' @seealso
#'  \code{\link[flowGraph]{flowGraph-class}}
#'  \code{\link[flowGraph]{fg_plot}}
#'  \code{\link[flowGraph]{plot_gr}}
#'  \code{\link[flowGraph]{fg_get_feature}}
#'  \code{\link[flowGraph]{fg_get_summary}}
#'  \code{\link[flowGraph]{fg_plot_qq}}
#'  \code{\link[flowGraph]{fg_plot_pVSdiff}}
#'  \code{\link[flowGraph]{fg_plot_box}}
#' @rdname fg_save_plots
#' @importFrom stringr str_pad
#' @importFrom utils head tail
#' @export
fg_save_plots <- function(
    fg, plot_path, plot_types="node", interactive=FALSE,
    adjust_custom="byLayer",
    label_max=10, box_no=20, paired=FALSE, logged=FALSE,
    filter_adjust0=1, filter_es=0,
    filter_btwn_tpthres=1, filter_btwn_es=0, overwrite=TRUE, ...
) {
    for (type in plot_types) {
        if (base::is.null(fg@summary[[type]])) next
        for (index in base::seq_len(base::length(fg@summary[[type]]))) {
            sm <- unlist(fg@summary_desc[[type]][index,])
            sm[2] <- base::paste0(sm[2], "-", ifelse(
                    is.function(adjust_custom),"adjusted",adjust_custom))
            plot_path_ <- base::paste0(
                plot_path, "/", type, "/", base::paste0(sm, collapse="_"))
            while (dir.exists(plot_path_) & !overwrite)
                plot_path_ <- base::paste0(plot_path_,"_")
            dir.create(plot_path_, recursive=TRUE, showWarnings=FALSE)

            try ({
                ## P value vs SpecEnr difference
                gp <- fg_plot_pVSdiff(
                    fg, type=type, summary_meta=NULL,
                    adjust_custom=adjust_custom, index=index,
                    filter_adjust0=filter_adjust0, filter_es=filter_es,
                    filter_btwn_tpthres=filter_btwn_tpthres,
                    filter_btwn_es=filter_btwn_es,
                    label_max=label_max, interactive=interactive,
                    path=base::paste0(plot_path_, "/pVSdifference.png"))
            })

            try ({
                ## boxplots for each node
                if (box_no>0) {
                pp = fg_get_summary(
                    fg, type="node", index=index, adjust_custom=adjust_custom,
                    filter_adjust0=filter_adjust0, filter_es=filter_es,
                    filter_btwn_tpthres=filter_btwn_tpthres,
                    filter_btwn_es=filter_btwn_es)
                seda = abs(pp$m1 - pp$m2)
                p = pp$values
                node_edges = union(
                    base::names(utils::head(base::sort(p),box_no)),
                    base::names(utils::tail(base::sort(seda),box_no)))
                rdir_ <- base::paste0(plot_path_, "/boxplots")
                dir.create(rdir_, recursive=TRUE, showWarnings=FALSE)
                for (node_edgei in base::seq_len(base::length(node_edges))) {
                    node_edge = node_edges[node_edgei]
                    if (node_edge=="") next
                    # if (p[node_edge]==1) next
                    # node_edge = "A+B+C+"
                    path <- base::paste0(
                        rdir_,"/",stringr::str_pad(node_edgei,width=3,pad="0"),
                        "_",node_edge,".png")
                    gp <- fg_plot_box_set(
                        fg, type="node", index=index, node_edge=node_edge,
                        adjust_custom=adjust_custom,
                        path=path, paired=paired,
                        filter_adjust0=filter_adjust0, filter_es=filter_es,
                        filter_btwn_tpthres=filter_btwn_tpthres,
                        filter_btwn_es=filter_btwn_es)
                }
                }
            })

            try ({
                ## qq plot
                gp <- fg_plot_qq(
                    fg, type=type, index=index, adjust_custom=adjust_custom,
                    filter_adjust0=filter_adjust0, filter_es=filter_es,
                    filter_btwn_tpthres=filter_btwn_tpthres,
                    filter_btwn_es=filter_btwn_es,
                    path=base::paste0(plot_path_,"/qq.png"),
                    interactive=interactive, logged=logged)
            })
#
#             try ({
#                 ## ecdf plot
#                 p <- fg_get_summary(fg, type=type, index=index)$values
#                 fun_ecdf <- ecdf(sort(p)); py <- fun_ecdf(sort(p))
#                 # df <- data.frame(y=sort(p),x=py)
#                 # a <- elbow(df)
#                 png(paste0(plot_path_,"/eCDF.png"))
#                 plot.ecdf(p, main=paste0("p-value eCDF plot"),
#                           xlab="p-values")
#                 # abline(v=a[1])
#                 dev.off()
#             })

            if (type=="node")
                try ({
                    ## cell hierarchy plot
                    # p <- fg_get_summary(fg, type=type, index=index)$values
                    gr <- fg_plot(
                        fg, type="node", index=index,
                        adjust_custom=adjust_custom,
                        path=base::paste0(plot_path_,"/cell_hierarchy.png"),
                        label_max=label_max, interactive=interactive, ...)

                })
        }
    }
}
