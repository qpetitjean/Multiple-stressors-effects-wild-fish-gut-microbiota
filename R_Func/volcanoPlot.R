#' Volcano Plot for LINDA Analysis Results
#'
#' @description
#' This function creates a volcano plot to visualize differential abundance 
#' results (e.g., from LINDA analysis) by plotting log2 fold changes versus 
#' -log10 adjusted p-values.
#'
#' @param x Numeric vector of log2 fold change values.
#' @param y Numeric vector of -log10 adjusted p-values.
#' @param xlim Numeric vector of length 2 specifying x-axis limits; if NULL, 
#'   computed automatically.
#' @param ylim Numeric vector of length 2 specifying y-axis limits; if NULL, 
#'   computed automatically.
#' @param group Vector indicating group membership for each point (for legend).
#' @param col Color specification for points (can be a vector).
#' @param main Character string for the main title of the plot.
#' @param xlab Character string for the x-axis label.
#' @param ylab Character string for the y-axis label.
#' @param yalign Character string ("left" or "right") to set the alignment of 
#'   the y-axis labels.
#' @param cex.lab Numeric value specifying the text size for axis labels.
#' @param cex.main Numeric value specifying the text size for the main title.
#' @param cex.axis Numeric value specifying the text size for axis tick labels.
#' @param cex.pts Numeric value specifying the size of the plotted points.
#' @param pch Integer or character symbol to use for plotting points.
#' @param srt Numeric value specifying the rotation angle for text labels.
#' @param legend Logical indicating whether to draw a legend.
#' @param cex.leg Numeric value specifying the text size in the legend.
#' @param leg.order Optional vector specifying the order of legend items.
#' @param hlines Numeric vector specifying horizontal line positions (y values).
#' @param vlines Numeric vector specifying vertical line positions (x values).
#' @param vlinesW Numeric value or vector for vertical line widths.
#' @param vlinesCol Color or vector of colors for vertical lines.
#' @param hlinesW Numeric value or vector for horizontal line widths.
#' @param hlinesCol Color or vector of colors for horizontal lines.
#' @param vlinesLty Numeric value or vector specifying the line type for vertical lines.
#' @param hlinesLty Numeric value or vector specifying the line type for horizontal lines.
#' @param labels Optional character vector of labels for all data points.
#' @param selLabels Optional character vector of labels to be displayed.
#' @param cex.labels Numeric value specifying the size of the text labels.
#' @param pos.label Numeric vector of length 2 controlling the position adjustment 
#'   of the text labels.
#'
#' @details
#' If \code{xlim} or \code{ylim} are not provided, they are computed by adding 
#' a 5% margin to the minimum and maximum values of \code{x} and \code{y}, respectively.
#' The function also supports drawing custom reference lines (both vertical and 
#' horizontal) with adjustable widths, colors, and line types. Optionally, selected 
#' data points can be labeled using a helper function from the \code{basicPlotteR} 
#' package.
#' 
#' **Dependencies:**
#' This function depends on several packages. It internally calls functions from:
#' \itemize{
#'   \item \strong{basicPlotteR} (To add text label)
#'   }
#'   
#' @importFrom basicPlotteR addTextLabels
#'
#' @return
#' This function produces a volcano plot.
#' 
#' @author Quentin PETITJEAN [quentin.petitjean@inrae.fr]
#'
#' @date 15/06/2023
#'
#' @export

volcanoPlot <- function(x = NULL,
                        y = NULL,
                        xlim = NULL,
                        ylim = NULL,
                        group = NULL,
                        col = "black",
                        main = "",
                        xlab = "x",
                        ylab = "y",
                        yalign = "left",
                        cex.lab = 1.5,
                        cex.main = 1,
                        cex.axis = 1.2,
                        cex.pts = 2,
                        pch = 19,
                        srt = 0,
                        legend = TRUE,
                        cex.leg = 1,
                        leg.order = NULL,
                        hlines = NULL,
                        vlines = NULL,
                        vlinesW = 0.2,
                        vlinesCol = adjustcolor("#808080", alpha.f = 0.5),
                        hlinesW = 0.2,
                        hlinesCol = adjustcolor("#808080", alpha.f = 0.5),
                        vlinesLty = 2,
                        hlinesLty = 2,
                        labels = NULL,
                        selLabels = NULL,
                        cex.labels = 1,
                        pos.label = c(1, 1)) {
  # if xlim is unspecified retrieve it approximately using the maximum value in x coordinates
  if (is.null(xlim)) {
    xlim <-
      c(signif(min(x, na.rm = T) + 5 * min(x, na.rm = T) / 100, 2), signif(max(x, na.rm = T) + 5 * max(x, na.rm = T) / 100, 2))
  }
  if (is.null(ylim)) {
    ylim <-
      c(signif(min(y, na.rm = T) + 5 * min(y, na.rm = T) / 100, 2), signif(max(y, na.rm = T) + 5 * max(y, na.rm = T) / 100, 2))
  }
  
  
  # initialize the plot window
  ScaleY <- pretty(ylim, n = 5)
  ScaleX <- pretty(xlim, n = 5)
  
  plot(
    NA,
    xlim = c(ScaleX[1], ScaleX[length(ScaleX)]),
    ylim = c(ScaleY[1], ScaleY[length(ScaleY)]),
    ylab = "",
    xlab = "",
    axes = FALSE
  )
  
  # add plot title
  graphics::mtext(
    main,
    cex = cex.main,
    side = 3,
    line = 1.5,
    adj = c(0.5, 0.5)
  )
  
  # add axis labels
  for (i in c("ylab", "xlab")) {
    graphics::mtext(
      get(i),
      cex = cex.lab,
      side = ifelse(i == "xlab", 1, 2),
      line = 1,
      adj = c(0.5, 0.5)
    )
  }
  
  # draw axes
  for (j in c("ScaleX", "ScaleY")) {
    graphics::segments(
      x0 = ifelse(j == "ScaleY", ifelse(yalign == "left", min(ScaleX), 0), get(j)[1]),
      y0 = ifelse(j == "ScaleY", get(j)[1], 0),
      x1 = ifelse(j == "ScaleY", ifelse(yalign == "left", min(ScaleX), 0), max(get(j))),
      y1 = ifelse(j == "ScaleY", max(get(j)), 0)
    )
    graphics::text(
      ifelse(j == rep("ScaleY", length(get(
        j
      ))), rep(
        ifelse(yalign == "left", min(ScaleX), 0), length(get(j))
      ), get(j)),
      ifelse(j == rep("ScaleY", length(get(
        j
      ))), get(j), rep(0, length(get(
        j
      )))),
      get(j),
      xpd = TRUE,
      srt = ifelse(length(srt) > 1, srt[2], srt),
      cex = cex.axis,
      pos = ifelse(j == "ScaleY", 2, 1)
    )
    graphics::text(
      ifelse(j == rep("ScaleY", length(get(
        j
      ))), rep(
        ifelse(yalign == "left", min(ScaleX), 0), length(get(j))
      ), get(j)),
      ifelse(j == rep("ScaleY", length(get(
        j
      ))), get(j), rep(0, length(get(
        j
      )))),
      "-",
      xpd = TRUE,
      srt = ifelse(j == "ScaleY", 0, 90),
      adj = c(0.8, 0.25)
    )
  }
  
  # draw horizontal
  if (!is.null(hlines)) {
    if (!length(hlinesW) == length(hlines)) {
      hlinesW <- rep(hlinesW, length(hlines))
    }
    if (!length(hlinesCol) == length(hlines)) {
      hlinesCol <- rep(hlinesCol, length(hlines))
    }
    if (!length(hlinesLty) == length(hlines)) {
      hlinesLty <- rep(hlinesLty, length(hlines))
    }
    for (h in seq_along(hlines)) {
      graphics::segments(
        x0 = min(ScaleX),
        y0 = hlines[[h]],
        x1 = max(ScaleX),
        y1 = hlines[[h]],
        lwd = hlinesW[[h]],
        col = hlinesCol[[h]],
        lty = hlinesLty[[h]]
      )
    }
  }
  # draw vertical lines
  if (!is.null(vlines)) {
    if (!length(vlinesW) == length(vlines)) {
      vlinesW <- rep(vlinesW, length(vlines))
    }
    if (!length(vlinesCol) == length(vlines)) {
      vlinesCol <- rep(vlinesCol, length(vlines))
    }
    if (!length(vlinesLty) == length(vlines)) {
      vlinesLty <- rep(vlinesLty, length(vlines))
    }
    for (v in seq_along(vlines)) {
      graphics::segments(
        x0 = vlines[[v]],
        y0 = min(ScaleY),
        x1 = vlines[[v]],
        y1 = max(ScaleY),
        lwd = vlinesW[[v]],
        col = vlinesCol[[v]],
        lty = vlinesLty[[v]]
      )
    }
  }
  
  # draw the points
  graphics::points(
    x = x,
    y = y,
    pch = pch,
    col = col,
    cex = cex.pts
  )
  
  # add labels
  if (!is.null(labels)) {
    if (is.null(selLabels)) {
      selLabels <- labels
    }
    LabPos <- which(labels %in% selLabels)
    xlabel <- x[LabPos]
    ylabel <- y[LabPos]
    labelCol <- col[LabPos]
    
    basicPlotteR::addTextLabels(
      xlabel,
      ylabel,
      selLabels,
      cex.label = cex.labels,
      col.label = labelCol,
      col.line = labelCol,
      col.background = adjustcolor("#FFFFFF", alpha.f = 0.8),
      border = adjustcolor(labelCol, alpha.f = 0.4)
    )
  }
  # add legend
  if (isTRUE(legend)) {
    GC <- data.frame(group, col)
    GC <- GC[!duplicated(GC), ]
    if(is.null(leg.order)){
    GC <- GC[order(GC$group),]
    }else{
    GC <- GC[match(leg.order, GC$group),]
    GC <- na.omit(GC)
    }
    
    legend_image <- graphics::legend(
      x = min(ScaleX) + 5 * max(ScaleX, na.rm = T) / 100,
      y = max(ScaleY) + 5 * max(ScaleY, na.rm = T) / 100,
      horiz = T,
      legend = GC[["group"]],
      pch = 19,
      col = GC[["col"]],
      pt.bg = GC[["col"]],
      bty = "n",
      pt.cex = 2.5,
      text.font = 1,
      cex = cex.leg,
      xpd = TRUE
    )
  }
  
}
