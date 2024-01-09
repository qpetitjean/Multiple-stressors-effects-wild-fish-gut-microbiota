################################################################################
# This is function aims to draw volcano plot                                   #
# Title:                                                                       #
# author:  "Quentin PETITJEAN[q.petitjean1@gmail.com]                          #
# date: "21/06/2023"                                                           #
################################################################################


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
