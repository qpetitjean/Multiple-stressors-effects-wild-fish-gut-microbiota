#' Customized Boxplot with Density Overlay
#'
#' @description
#' \code{Bxplt} generates a customized boxplot for visualizing the distribution of a numeric variable 
#' across different groups. In addition to drawing standard boxplots, the function overlays jittered data 
#' points and (optionally) density curves for each group, and annotates group means with dashed lines.
#'
#' @param x A vector or factor of group labels.
#' @param y A numeric vector of data values corresponding to each element of \code{x}.
#' @param fill A vector of colors used to fill the boxplots. Defaults to "white" if not specified.
#' @param colpts A vector of colors for the jittered points. Defaults to "darkgrey" if not specified.
#' @param boxwex A numeric value specifying the width of the boxplots. Default is 0.5.
#' @param xlab A character string specifying the label for the x-axis. Defaults to an empty string.
#' @param ylab A character string specifying the label for the y-axis. Defaults to an empty string.
#' @param las An integer controlling the style of axis labels. Default is 1 (horizontal).
#' @param cex.axis A numeric value specifying the character expansion factor for the axis annotation. Default is 0.6.
#' @param Jittering A numeric value controlling the amount of jitter applied to x-axis positions. Default is 1.
#' @param main A character string specifying the plot title. Defaults to an empty string.
#' @param ymin A numeric value for the minimum y-axis limit. If \code{NULL}, computed from the data.
#' @param yscale Either a single numeric value (number of tick marks) or a vector of numeric tick mark positions for the y-axis. Default is 5.
#' @param Density A logical value indicating whether to overlay density curves on the boxplot. Default is \code{TRUE}.
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
#' This function display a customized boxplot
#' 
#' @author Quentin PETITJEAN [quentin.petitjean@inrae.fr]
#'
#' @date 15/06/2023
#' @export

Bxplt <-
  function(x,
           y,
           fill = NULL,
           colpts = NULL,
           boxwex = 0.5,
           xlab = NULL,
           ylab = NULL,
           las = 1,
           cex.axis = 0.6, 
           Jittering = 1,
           main = NULL,
           ymin = NULL,
           yscale = 5,
           Density = T) {
    if (is.null(xlab)) {
      xlab <- ""
    }
    if (is.null(ylab)) {
      ylab <- ""
    }
    if (is.null(fill)) {
      fill <- "white"
    }
    if (is.null(colpts)) {
      colpts <- "darkgrey"
    }
    if (is.null(ymin)) {
      ymin <- min(y, na.rm = T)
    }
    Mean <- tapply(y, x, function(x) mean(x, na.rm = T))
    ymax <- max(y, na.rm = T)
    iter <- 0

    if(!is.null(yscale) && length(yscale) == 1){
      while(ymax < yscale) {
        iter <- iter + 1
        ymax <- ymax * 10^1
      } 
      while(ymax %% yscale != 0) {
        ymax <- ymax + ifelse(ymax %% 1 == 0, 1, ceiling(ymax) - ymax)
      } 
      Yscale <- seq(ifelse(ymin < 0, floor(ymin), ceiling(ymin)), ceiling(ymax), length.out = yscale)
      Yscaleval <- sapply(ifelse(rep(Yscale[1] < 0, length(Yscale)), Yscale, Yscale[-1]), function(x) {
      while(x %% yscale != 0) {
        if(x >= 0){
        x <- x + ifelse(x %% 1 == 0, 1, ceiling(x) - x)
        }else{
          x <- x + ifelse(x %% 1 == 0, -1, x - ceiling(x) )
        }
      } 
        return(x)
      })
      Yscaleval <- append(ifelse(ymin < 0, floor(ymin), ceiling(ymin)), Yscaleval)
      if(iter != 0){
        Yscaleval <- Yscaleval * 10^-iter
      }
      if(ymin < Yscaleval[1]){
        Yscaleval[1] <- floor(ymin)
      }
      if(ymin > Yscaleval[1]){
        Yscaleval[1] <- ymin
      }
      Yscaleval <- pretty(Yscaleval)
    }else if (!is.null(yscale) && length(yscale) > 1){
      Yscaleval = yscale
    }
    
    XscaleLab <- unique(x)[!is.na(unique(x))]
    XscaleLab <- sort(XscaleLab)
    Xscaleval <- seq_along(XscaleLab)
    
  bp <- boxplot(
      y~x,
      staplewex = 0,
      boxwex = boxwex,
      whisklty = 1,
      las = las,
      cex.axis = cex.axis,
      xlab = "",
      col = fill,
      ylab = "",
      outline = F,
      frame = F,
      axes = F,
      main = main,
      ylim = c(Yscaleval[1], Yscaleval[length(Yscaleval)]),
      xlim = c(0, length(Xscaleval) + 1)
      #yaxt = ifelse(is.null(yscale),"s", "n")
    )
    
  mtext(ylab, side = 2, line = 1.5)
  mtext(xlab, side = 1, line = 1.5)
  
    # draw x axis
    segments(
      x0 = 0,
      x1 = length(Xscaleval) + 1,
      y0 = Yscaleval[1],
      y1 = Yscaleval[1]
    )
    
    text(
      Xscaleval,
      rep(Yscaleval[1], length(Xscaleval)),
      "-",
      xpd = TRUE,
      srt = 90,
      cex = 1,
      adj = c(1,0.25))
    
    text(
      Xscaleval,
      rep(Yscaleval[1], length(Xscaleval)),
      XscaleLab,
      xpd = TRUE,
      srt = ifelse(las == 1, 0, 90),
      pos = 1,
      offset = 1,
      cex = cex.axis
    )
    
    # draw Y axis
    segments(
      x0 = 0,
      x1 = 0,
      y0 = Yscaleval[1],
      y1 = Yscaleval[length(Yscaleval)] + diff(Yscaleval)[1]
    )
    
    text(
      rep(0, length(Yscaleval)),
      Yscaleval,
      "-",
      xpd = TRUE,
      srt = 0,
      cex = 1,
      adj = c(1, 0.25))
    
    text(
      rep(0 - (1/10 * diff(Xscaleval)[1]), length(Yscaleval)),
      Yscaleval,
      Yscaleval,
      xpd = TRUE,
      srt = ifelse(las == 1, 0, 90),
      pos = 2,
      offset = 0,
      cex = cex.axis
    )
 
    JitNum <-
      match(x, XscaleLab)
    Jit <- jitter(JitNum, Jittering)
    for (k in seq(length(x))) {
      points(Jit[k],
             y[k],
             pch = 19,
             col = colpts[k],
             cex = 1)
    }
    for (l in Xscaleval) {
      segments(
        x0 = l - boxwex / 2,
        x1 = l + boxwex / 2,
        y0 = Mean[[l]],
        y1 = Mean[[l]],
        lty = "dashed"
      )
    }
    
    # add density lines 
    if(isTRUE(Density)){
      par(new = TRUE) 
      plot(
        0,
        frame = F,
        type = "n",
        ylim = c(Yscaleval[1], Yscaleval[length(Yscaleval)]),
        xlim = c(0, length(Xscaleval) + 1),
        xaxt = "n",
        yaxt = "n",
        xlab = "",
        ylab = ""
      )
      
      df <- data.frame(x,y)
      df <- df[!is.na(df$y),]
      
      densities <- lapply(split(df$y, df$x), density)
      Xmax <- max(unlist(lapply(densities, function(x) max(x$y)))) 
      
      for(i in seq_along(densities)) {
        d <- densities[[i]]
        xVal <- (d$y / Xmax) * 0.25 + i + boxwex / 2 
        yVal <- d$x
        yVal <- pmax(pmin(yVal, max(bp$stats[5,])), min(bp$stats[1,]))

        graphics::polygon(x = c(i + boxwex / 2, xVal, i + boxwex / 2),
                          y = c(min(yVal), yVal, max(yVal)),
                          col = fill[i], border = "black")
      }
      
    }
  }

