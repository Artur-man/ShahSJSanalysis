#' plot_violin_seurat - modification of function from Garber Lab SignallingSingleCell package for input of Seurat Object
#' visualizing multiple genes
#' 
#' @description Create violin plot.
#' @param input SatijaLab’s Seurat Class, with normalized expression values in assay data slot. Supports RNA and ADTs via Seurat::FetchData()
#' @param title Title of the graph. Would be the gene name if not specified
#' @param gene Feature for which to plot the expression level.
#' @param color_by a meta.data column variable.
#' @param log_scale If true, transform UMIs by log2(UMI + 1).
#' @param colors What colors to utilize for categorical data. Be sure it is of the proper length.
#' @param color_order a vector of the order of variables to display on the x axis. Ensure that factors have been re-leveled if dataset has been subsetted
#' @param facet_by a vector with one or two meta.data column variables. If two, the first variable as columns and the second as rows.
#' @param facet_order a vector with the desired order of the first facet column variable. Only one column variable is currently supported.
#' @param wrapfacet_cols an integer for how many columns to wrap plots. Requires facet_by to be of order 1.
#' @param wrapfacet_scales a string to indicate scale for axes in wrapped plots. Default "free_y". May be "fixed", "free", "free_x", "free_y"
#' @param spread e.g. Healthy category is unique in Disease and Skin. To use Healthy only as skin but not Disease, that is adding Healthy skin to each disease, spread = c("Disease", "Healthy").
#' @param title_size Size of title, default is 20
#' @param axis_title_size Size of axis title, default is 10
#' @param axis_text_size Size of text of tick mark on the axis, default is 5
#' @param legend_title_size Size of legend title, default is 10
#' @param legend_text_size Size of text of legend, default is 5
#' @param facet_text_size Size of text of facet box, default is 5
#' @param number_label_text_size Size of the numbers below the violin plots, default is 2
#' @param theme the plot theme. Default to be "classic" if not set to "bw".
#' @param number_labels show the total cell numbers and cell fraction with non-zero expression values under each bar.
#' @param plot_mean plot the mean value as black dot with second y-axis on the right.
#' @param size the size of dots.
#' @param sig the number of digits after the decimal point for cell fraction value.
#' @param contour_line_width the thickness of the violin contour line
#' @details
#' Utilize information stored in meta.data to control the plot display. Each point_by as a dot with a bar showing the weighted mean of all point_by dots.
#' @importFrom reshape2 melt
#' @examples
#' plot_violin_seurat(ex_sc, gene = "CD8A", color_by = "Skin", facet_by = c("Disease", "CellType"), log_scale = F)
#' plot_violin_seurat(ex_sc, gene = "adt_CD8", color_by = "Skin", facet_by = c("Disease", "CellType"), log_scale = F)
#' plot_violin_seurat(ex_sc, gene = "CXCL13", color_by = "Skin", facet_by = c("CellType", "Disease"), spread = T, log_scale = T)
#' @export
plot_violin_seurat <- function (input, title = "", gene, color_by, log_scale = F, color_order = NULL,
                                colors = NULL, facet_by = NULL, facet_order = NULL, wrapfacet_cols = NULL, wrapfacet_scales = 'free_y', 
                                spread = NULL, jitter_pts = T,
                                plot_mean = T, plot_mean_dot_size = 2, size = 1, sig = 3, number_labels = T,
                                title_size = 20, axis_title_size = 10, axis_text_size = 5,
                                legend_title_size = 10, legend_text_size = 5, facet_text_size = 5,
                                number_label_text_size = 2, alpha = 0.5, theme = "classic",
                                contour_line_width = 0.3)
{
  text_sizes = c(title_size, axis_title_size, axis_text_size, legend_title_size, legend_text_size,
                 2, number_label_text_size,facet_text_size)
  #title_size, axis_title, axis_text, legend_title, legend_text, facet_text, number_label_text_size,
  df <- input@meta.data[,colnames(input@meta.data) %in% c(color_by, facet_by), drop = F]
  df <- cbind(df, FetchData(input,gene))
  colnames(df) <- gsub("-", "", colnames(df))
  gene <- gsub("-", "", gene)
  if (any(!is.null(spread))) {
    others <- setdiff(unique(df[,spread[1]]), spread[2])
    ind <- which(df[, spread[1]] == spread[2])
    rmdf <- df[ind,]
    df <- df[-ind,]
    for (i in 1:length(others)) {
      rmdf[,spread[1]] <- others[i]
      df <- rbind(df, rmdf)
    }
  }
  
  if(length(gene) > 1){
    df <- melt(df, measure.vars = gene) 
    facet_by <- c("variable", facet_by)
  } else {
    colnames(df)[colnames(df) %in% gene] <- "value"
  }
  
  if (log_scale) {df$plot <- log2(df$value + 1)}else{df$plot <- df$value}
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  cols <- gg_color_hue(length(unique(df[, color_by])))
  
  if(all(!is.null(color_order))) {
    df <- df %>% arrange(factor(.data[[color_by]], levels = color_order, ordered = TRUE))
  }
  
  #facet order currently only supported for 1 dimension of facet
  if(all(!is.null(facet_order))) {
    df <- df %>% mutate(across(facet_by, factor, levels=facet_order))
  }
  
  
  g <- ggplot(df)
  if (all(!is.null(colors))) {
    g <- g + scale_color_manual(values = c(colors))
    g <- g + scale_fill_manual(values = c(colors))
  }
  if (theme == "bw") {
    g <- g + theme_bw()
  }else{
    g <- g + theme_classic()
  }
  if (title == "") title <- gene
  g <- g + labs(title = title, y = "Expression")
  g <- g + theme(plot.title = element_text(size = text_sizes[1]),
                 axis.title = element_text(size = text_sizes[2]), axis.text = element_text(size = text_sizes[3]),
                 legend.title = element_text(size = text_sizes[4]), legend.text = element_text(size = text_sizes[5]), 
                 strip.text = element_text(size = text_sizes[8]))
  g <- g + theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.text.x = element_blank(),
                 axis.ticks.x = element_blank())
  if (jitter_pts) g <- g + geom_jitter(aes_string(x = color_by, y = "plot", col = color_by), width = 0.2, size = size)
  g <- g + geom_violin(aes_string(x = color_by, y = "plot", fill = color_by), col = "black", trim = T, scale = "width", alpha = alpha, size=contour_line_width) 
  if (number_labels) {
    g <- g + stat_summary(aes_string(x = color_by, y = "value"), fun.data = function(x) {return(c(y = -max(df$plot)/25, label = length(x)))}, colour = "black",
                          geom = "text", size = text_sizes[7])
    g <- g + stat_summary(aes_string(x = color_by, y = "value"), fun.data = function(x) {return(c(y = -max(df$plot)/10, label = round(mean(as.numeric(x > 0)), sig)))}, colour = "black",
                          geom = "text", size = text_sizes[7])
  }
  if (plot_mean) {
    scale <- max(df$plot)/max(tapply(df$value, INDEX = as.list(df[, colnames(df) %in% c(color_by, facet_by), drop = F]), FUN=mean), na.rm = T)
    g <- g + suppressWarnings(stat_summary(aes_string(x = color_by, y = "value"), fun.y = function(x) mean(x)*(scale * 0.5), colour = "black", geom = "point", size = plot_mean_dot_size))
    g <- g + scale_y_continuous(sec.axis = sec_axis(~./(scale * 0.5), name = "Mean Expression"))
  }
  
  if (length(facet_by) == 1) {
    if(!is.null(wrapfacet_cols)){
      g <- g + facet_wrap(facets = reformulate(facet_by), ncol = wrapfacet_cols, scales = wrapfacet_scales)
    }
    else{
      g <- g + facet_grid(facets = reformulate(facet_by), scales = "free_y")
    }
  }else if (length(facet_by) == 2) {
    g <- g + facet_grid(facets = reformulate(facet_by[1], facet_by[2]), scales = "free_y")
  }else if (length(facet_by) > 2) {stop("Parameter facet_by needs to be a string with equal or less than two variables.")}
  return(g)
}