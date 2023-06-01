
# -------------------------------------------------------------------------------------------

#' Plots the trend series and the (fitted) second observation equation and gives diagnostic
#' plots based on standardized residuals.
#'
#' @param tsl A list with two multiple time series objects for the first and second plor,
#'   respectively.
#' @param legend A list with two character vectors. The first contains the legend names for
#'   the first plot and so on.
#' @param title A list with the titles for the first three plots.
#' @param boundName The legend name of the confidence bounds.
#' @param res The residual series as time series. If \code{res = NULL}, all graphs realted
#'   to the residual series will not be plotted.
#' @param namesPrint A character vector containing two names for the first two plots. The
#'   remaining names are creates automatically if \code{combine = FALSE}.
#' @inheritParams plot.NAWRUfit
#'
#' @importFrom stats start end window ts lag frequency time na.pass density acf qnorm
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @keywords internal
plotSSresults <- function(tsl, legend, title, boundName, res = NULL, namesPrint, bounds,
                          combine, path, device, width, height) {

  # to suppress R CMD check notes.
  trend <- orig <- anchor <- lb <- ub <- fitted <- E2 <- residuals <- ..density.. <- NULL

  # settings
  setPrint <- list(width = width, height = height)
  setPrintComb <- list(width = width, height = height * 2)
  set <- list()
  set$titleFontsize <- 10
  set$labelFontsize <- 8
  set$legendFontsize <- 8
  set$textSize <- 3.5
  set$colors <- c("black", "black", "black", "grey12", "grey24", "grey36", "grey48", "grey60")
  set$linetype <- c(2, 1, 3, 4, 5, 6)
  set$alpha <- 0.2
  set$freqYear <- 1 + floor(length(tsl[[1]][, 1]) / frequency(tsl[[1]][, 1]) / 15)
  set$freqYearSmall <- ceiling(set$freqYear * 2)
  if (!combine) {
    set$freqYearSmall <- set$freqYear
  }

  # ----- trend
  date <- zoo::as.Date(time(tsl[[1]]))
  df <- data.frame(date, tsl[[1]])
  color <- set$colors[1:(ncol(df) - 3)]
  color <- color[order(legend[[1]][1:(ncol(df) - 3)])]
  linetype <- set$linetype[1:(ncol(df) - 3)]
  linetype <- linetype[order(legend[[1]][1:(ncol(df) - 3)])]
  p0 <- ggplot(df, aes(x = date)) +
    theme_classic() +
    coord_cartesian(xlim = c(date[1], date[length(date)]), expand = FALSE) +
    theme(
      plot.title = element_text(size = set$titleFontsize),
      axis.title = element_text(size = set$labelFontsize),
      panel.grid.major.y = element_line(linetype = "solid")
    ) +
    geom_line(aes(y = trend, col = legend[[1]][1], linetype = legend[[1]][1])) +
    geom_line(aes(y = orig, col = legend[[1]][2], linetype = legend[[1]][2]))
  if (ncol(df) > 5) {
    p0 <- p0 + geom_line(aes(y = anchor, col = legend[[1]][3], linetype = legend[[1]][3]))
  }
  if (bounds) {
    p0 <- p0 +
      geom_ribbon(aes(ymin = lb, ymax = ub, fill = boundName), alpha = set$alpha) +
      scale_fill_manual(values = "grey12")
  }
  p0 <- p0 +
    scale_color_manual(name = NULL, values = color) +
    scale_linetype_manual(name = NULL, values = linetype) +
    scale_x_date(
      date_labels = "%Y", date_minor_breaks = "1 year",
      date_breaks = paste0(set$freqYear, " year")
    ) +
    labs(title = title[[1]], x = "year", y = "") +
  theme(
    legend.position="bottom",
    legend.box = "horizontal",
    legend.box.margin=margin(-15,0,0,0),
    legend.justification="left",
    legend.text = element_text(size = set$legendFontsize),
    legend.key.size = unit(0.75, "lines"),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.spacing.y = unit(0, "pt"),
    plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")
  ) +
    guides(color = guide_legend(order=1, ncol = 1, byrow = TRUE),
           linetype = guide_legend(order=1, ncol = 1, byrow = TRUE),
           fill = guide_legend(order=2)) 
  
  if (length(tsl) > 2) {
    # ----- trend growth
    date <- zoo::as.Date(time(tsl[[3]]))
    df <- data.frame(date, tsl[[3]])
    color <- set$colors[1:(ncol(df) - 3)]
    color <- color[order(legend[[3]][1:(ncol(df) - 3)])]
    linetype <- set$linetype[1:(ncol(df) - 3)]
    linetype <- linetype[order(legend[[3]][1:(ncol(df) - 3)])]
    p02 <- ggplot(df, aes(x = date)) +
      theme_classic() +
      coord_cartesian(xlim = c(date[1], date[length(date)]), expand = FALSE) +
      theme(
        plot.title = element_text(size = set$titleFontsize),
        axis.title = element_text(size = set$labelFontsize),
        panel.grid.major.y = element_line(linetype = "solid")
      ) +
      geom_line(aes(y = trend, col = legend[[3]][1], linetype = legend[[3]][1])) +
      geom_line(aes(y = orig, col = legend[[3]][2], linetype = legend[[3]][2]))
    if (ncol(df) > 5) {
      p02 <- p02 + geom_line(aes(y = anchor, col = legend[[3]][3], linetype = legend[[3]][3]))
    }
    if (bounds) {
      p02 <- p02 +
        geom_ribbon(aes(ymin = lb, ymax = ub, fill = boundName), alpha = set$alpha) +
        scale_fill_manual(values = "grey12")
    }
    p02 <- p02 +
      scale_color_manual(name = NULL, values = color) +
      scale_linetype_manual(name = NULL, values = linetype) +
      scale_x_date(
        date_labels = "%Y", date_minor_breaks = "1 year",
        date_breaks = paste0(set$freqYear, " year")
      ) +
      labs(title = title[[4]], x = "year", y = "") +
      theme(
        legend.position="bottom",
        legend.box = "horizontal",
        legend.box.margin=margin(-15,0,0,0),
        legend.justification="left",
        legend.text = element_text(size = set$legendFontsize),
        legend.key.size = unit(0.75, "lines"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.spacing.y = unit(0, "pt"),
        plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")
      ) +
      guides(color = guide_legend(order=1, ncol = 1, byrow = TRUE),
             linetype = guide_legend(order=1, ncol = 1, byrow = TRUE),
             fill = guide_legend(order=2)) 
  }
  
  # ----- second observation equation
  date <- zoo::as.Date(time(tsl[[2]]))
  df <- data.frame(date, tsl[[2]])
  color <- set$colors[1:2]
  color <- color[order(legend[[1]][1:2])]
  linetype <- set$linetype[1:(ncol(df) - 3)]
  linetype <- linetype[order(legend[[1]][1:(ncol(df) - 3)])]
  p1 <- ggplot(df, aes(x = date)) +
    theme_classic() +
    coord_cartesian(xlim = c(date[1], date[length(date)]), expand = FALSE) +
    theme(
      plot.title = element_text(size = set$titleFontsize),
      axis.title = element_text(size = set$labelFontsize),
      panel.grid.major.y = element_line(linetype = "solid")
    ) +
    geom_line(aes(y = fitted, col = legend[[2]][1], linetype = legend[[2]][1])) +
    geom_line(aes(y = E2, col = legend[[2]][2], linetype = legend[[2]][2]))
  if (bounds) {
    p1 <- p1 +
      geom_ribbon(aes(ymin = lb, ymax = ub, fill = boundName), alpha = set$alpha) +
      scale_fill_manual(values = set$colors[4])
  }
  p1 <- p1 +
    scale_color_manual(name = NULL, values = color) +
    scale_linetype_manual(name = NULL, values = linetype) +
    scale_x_date(
      date_labels = "%Y", date_minor_breaks = "1 year",
      date_breaks = paste0(set$freqYearSmall, " year")
    ) +
    labs(x = "year", y = "") +
  theme(
    legend.position="bottom",
    legend.box = "horizontal",
    legend.box.margin=margin(-15,0,0,0),
    legend.justification="left",
    legend.text = element_text(size = set$legendFontsize),
    legend.key.size = unit(0.75, "lines"),
    legend.title = element_blank(),
    legend.background = element_rect(fill = "transparent"),
    legend.spacing.y = unit(0, "pt"),
    plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")
  ) +
    guides(color = guide_legend(order=1, ncol = 1, byrow = TRUE),
           linetype = guide_legend(order=1, ncol = 1, byrow = TRUE),
           fill = guide_legend(order=2)) 

  # residual diagnostics
  if (!is.null(res)) {

    # ----- residuals
    dfRes <- data.frame(
      "residuals" = as.numeric(res),
      "date" = zoo::as.Date(time(res))
    )
    p2 <- ggplot(dfRes, aes(x = date)) +
      theme_classic() +
      coord_cartesian(xlim = c(dfRes$date[1], dfRes$date[length(dfRes$date)]), expand = FALSE) +
      theme(
        plot.title = element_text(size = set$titleFontsize),
        axis.title = element_text(size = set$labelFontsize),
        panel.grid.major.y = element_line(linetype = "solid")
      ) +
      geom_line(aes(y = residuals, col = "recursive residuals")) +
      scale_color_manual(name = NULL, values = set$colors[1]) +
      scale_x_date(
        date_labels = "%Y", date_minor_breaks = "1 year",
        date_breaks = paste0(set$freqYearSmall, " year")
      ) +
      labs(x = "year", y = "") +
    theme(
      legend.position="bottom",
      legend.box = "horizontal",
      legend.box.margin=margin(-15,0,0,0),
      legend.justification="left",
      legend.text = element_text(size = set$legendFontsize),
      legend.key.size = unit(0.75, "lines"),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent"),
      legend.spacing.y = unit(0, "pt"),
      plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")
    ) 

    # ----- acf
    ciline <- qnorm(0.05 / 2) / sqrt(length(na.trim(res)))
    bacf <- acf(na.trim(res), plot = FALSE)
    bacfdf <- with(bacf, data.frame(lag, acf))
    p3 <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
      theme_classic() +
      theme(
        plot.title = element_text(size = set$titleFontsize),
        axis.title = element_text(size = set$labelFontsize)
      ) +
      geom_hline(aes(yintercept = 0)) +
      geom_segment(mapping = aes(xend = lag, yend = 0)) +
      geom_hline(aes(yintercept = ciline), linetype = 3, color = set$colors[4]) +
      geom_hline(aes(yintercept = -ciline), linetype = 3, color = set$colors[4]) +
      labs(title = "Autocorrelation of recursive residuals")

    # ----- histogram
    data <- data.frame(na.trim(as.numeric(res)))
    colnames(data) <- "residuals"
    p4 <- ggplot(data = data, aes(x = residuals)) +
      theme_classic() +
      theme(
        plot.title = element_text(size = set$titleFontsize),
        axis.title = element_text(size = set$labelFontsize)
      ) +
      geom_histogram(aes(y = ..density..), bins = 20, color = set$colors[7], fill = set$colors[7], alpha = set$alpha) +
      geom_density() +
      labs(title = "Histogram", x = "recursive residuals", y = "count")

    # print
    if (combine) {
      if (length(tsl) > 2) {
        p0 <- suppressWarnings(gridExtra::grid.arrange(p0, p02, nrow = 1))
      } else {
        suppressWarnings(print(p0))
      }
      pCombine <- suppressWarnings(gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2, top = title[[2]]))
    } else {
      suppressWarnings(print(p0))
      if (length(tsl) > 2) {
        suppressWarnings(print(p02))
      }
      p1 <- p1 + labs(title = title[[2]])
      p2 <- p2 + labs(title = title[[3]])
      suppressWarnings(print(p1))
      suppressWarnings(print(p2))
      suppressWarnings(print(p3))
      suppressWarnings(print(p4))
    }

    # save files
    if (!is.null(path)) {
      do.call(ggsave, c(list(filename = file.path(path, paste0(namesPrint[1], ".", device)), plot = p0), setPrint))
      if (combine) {
        do.call(ggsave, c(list(filename = file.path(path, paste0(namesPrint[2], ".", device)), plot = pCombine), setPrintComb))
      } else {
        do.call(ggsave, c(list(filename = file.path(path, paste0(namesPrint[2], "_fitted.", device)), plot = p1), setPrint))
        do.call(ggsave, c(list(filename = file.path(path, paste0(namesPrint[2], "_residuals.", device)), plot = p2), setPrint))
        do.call(ggsave, c(list(filename = file.path(path, paste0(namesPrint[2], "_residuals_acf.", device)), plot = p3), setPrint))
        do.call(ggsave, c(list(filename = file.path(path, paste0(namesPrint[2], "_residuals_histogram.", device)), plot = p4), setPrint))
      }
    }
  } else {
    # no residuals

    # print
    if (combine && (length(tsl) > 2)) {
      p0 <- suppressWarnings(gridExtra::grid.arrange(p0, p02, nrow = 1))
    } else {
      suppressWarnings(print(p0))
    }
    p1 <- p1 + labs(title = title[[2]])
    suppressWarnings(print(p1))

    # save files
    if (!is.null(path)) {
      do.call(ggsave, c(list(filename = file.path(path, paste0(namesPrint[1], ".", device)), plot = p0), setPrint))
      do.call(ggsave, c(list(filename = file.path(path, paste0(namesPrint[2], ".", device)), plot = p1), setPrint))
    }
  }
}


# -------------------------------------------------------------------------------------------

#' Plots the diagnostic plots of the posterior distribution.
#'
#' @param draws the number of draws.
#' @param parName the name of the parameter.
#' @param burnin length of thw burnin phase.
#' @param mu mean for normal distribution
#' @param prec precision for normal distribution
#' @param shape shape for Gamma distribution
#' @param scale scale for Gamma distribution
#' @param shape1 shape1 for Beta distribution
#' @param shape2 shape2 for Beta distribution
#' @param lb lower bound
#' @param ub upper bound
#' @inheritParams plot.TFPfit
#'
#' @importFrom stats density qbeta dnorm dgamma
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @keywords internal
plot_gibbs_output <- function(path = NULL, draws, parName, burnin, mu = NULL, prec = NULL,
                              shape = NULL, scale = NULL, shape1 = NULL, shape2 = NULL,
                              ub = NULL, lb = NULL, prefix = NULL, device, width, height) {
  posterior <- prior <- iter <- th <- NULL

  drawsBurnin <- draws
  draws <- draws[(burnin + 1):length(draws)]

  # settings
  n <- 500
  # printing settings
  width <- width
  height <- height * 2
  pointsize <- 12
  set <- list()
  set$titleFontsize <- 10
  set$labelFontsize <- 8
  set$legendFontsize <- 8
  set$textSize <- 3.5
  set$colors <- c("black", "grey12", "grey24", "grey36", "grey48", "grey60")
  set$linetype <- c(2, 1, 3, 4, 5, 6)
  set$alpha <- 0.2

  # get limits for densities
  if (!is.null(mu)) {
    from <- min(min(draws), qnorm(0.2, mean = mu, sd = sqrt(1 / prec)))
    to <- max(max(draws), qnorm(0.8, mean = mu, sd = sqrt(1 / prec)))
  } else if (!is.null(shape)) {
    from <- min(min(draws), qgamma(0.2, shape = shape, rate = 1 / scale))
    to <- max(max(draws), qgamma(0.8, shape = shape, rate = 1 / scale))
  } else if (!is.null(shape1)) {
    const <- ifelse(lb > -Inf & ub < Inf, lb, 0)
    scale <- ifelse(lb > -Inf & ub < Inf, (ub - lb), 1)
    from <- min(min(draws), const + scale * qbeta(0.2, shape1 = shape1, shape2 = shape2))
    to <- max(max(draws), const + scale * qbeta(0.8, shape1 = shape1, shape2 = shape2))
  }

  # prior density
  grid <- seq(from, to, length = n)
  if (!is.null(mu)) {
    dprior <- dnorm(grid, mu, 1 / sqrt(prec))
    dprior[grid < lb] <- 0
    dprior[grid > ub] <- 0
  } else if (!is.null(shape)) {
    dprior <- dgamma(grid, shape = shape, rate = 1 / scale)
    dprior[grid < lb] <- 0
    dprior[grid > ub] <- 0
  } else if (!is.null(shape1)) {
    dprior <- ((grid - lb)^(shape1 - 1) * (ub - grid)^(shape2 - 1)) / ((ub - lb)^(shape1 + shape2 - 1) * beta(shape1, shape2))
  }

  # posterior density
  dpost <- density(draws, n = n, from = from, to = to)$y
  dpost[grid < lb] <- 0
  dpost[grid > ub] <- 0

  # ----- histogram
  data <- data.frame(draws)
  colnames(data) <- parName
  p1 <- ggplot(data = data, aes(x = get(parName))) +
    theme_classic() +
    theme(
      plot.title = element_text(size = set$titleFontsize),
      axis.title = element_text(size = set$labelFontsize)
    ) +
    geom_histogram(bins = 30, color = set$colors[1], fill = set$colors[4], alpha = set$alpha) +
    geom_vline(aes(xintercept = mean(draws)), color = set$colors[5], linetype = set$linetype[1]) +
    geom_vline(aes(xintercept = median(draws)), color = set$colors[3], linetype = set$linetype[3]) +
    labs(title = "Histogram of posterior distributions", x = parName, y = "count")
  # get position of text labels
  locText <- rep(0, 2)
  locText[1] <- ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range[2] * (1 - 0.005)
  locText[2] <- ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range[2] * (1 - 20 * 0.005)
  locX <- diff(ggplot_build(p1)$layout$panel_scales_x[[1]]$range$range) * 0.005
  # add text labels
  p1 <- p1 +
    annotate("text",
      x = mean(draws) + locX, y = locText[1], label = paste0("mean ", sprintf(mean(draws), fmt = "%#.3g")),
      hjust = 0, color = set$colors[5], size = set$textSize
    ) +
    annotate("text",
      x = median(draws) + locX, y = locText[2], label = paste0("median ", sprintf(median(draws), fmt = "%#.3g")),
      hjust = 0, color = set$colors[3], size = set$textSize
    )

  # ----- prior and posterior
  grid_dat <- data.frame(
    grid = grid,
    prior = dprior,
    posterior = dpost
  )
  p2 <- ggplot(grid_dat, aes(x = grid)) +
    theme_classic() +
    theme(
      plot.title = element_text(size = set$titleFontsize),
      axis.title = element_text(size = set$labelFontsize)
    ) +
    geom_area(aes(y = posterior, col = "posterior", fill = "posterior")) +
    geom_area(aes(y = prior, col = "prior", fill = "prior")) +
    scale_fill_manual(name = NULL, values = alpha(set$colors[c(1, 5)], set$alpha)) +
    scale_color_manual(name = NULL, values = set$colors[c(1, 5)]) +
    labs(title = "Prior and posterior distribution", x = parName, y = "density", col = "") +
    theme(
      legend.position = c(0.8, 0.95),
      legend.text = element_text(size = set$legendFontsize),
      legend.key.size = unit(0.75, "lines"),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent")
    ) +
    guides(fill = guide_legend(override.aes = list(alpha = set$alpha)))

  # ----- trace plot
  data <- data.frame("iter" = 1:length(drawsBurnin), "th" = drawsBurnin)
  p3 <- ggplot(data, aes(x = iter, y = th)) +
    theme_classic() +
    coord_cartesian(xlim = c(0, 1.02 * length(drawsBurnin)), expand = FALSE) +
    theme(
      plot.title = element_text(size = set$titleFontsize),
      axis.title = element_text(size = set$labelFontsize)
    ) +
    geom_line() +
    geom_vline(aes(xintercept = burnin), color = set$colors[2], linetype = set$linetype[1]) +
    annotate("text",
      x = burnin + length(drawsBurnin) / 100, y = range(draws)[1] + 0.95 * diff(range(draws)),
      label = "burnin", hjust = 0, color = set$colors[2], size = set$textSize
    ) +
    labs(title = "Trace plot", x = "iterations", y = parName)

  # ----- acf
  p4 <- NULL
  if (var(draws) != 0) {
    ciline <- qnorm(0.05 / 2) / sqrt(length(draws))
    bacf <- acf(draws, plot = FALSE)
    bacfdf <- with(bacf, data.frame(lag, acf))
    p4 <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
      theme_classic() +
      theme(
        plot.title = element_text(size = set$titleFontsize),
        axis.title = element_text(size = set$labelFontsize)
      ) +
      geom_hline(aes(yintercept = 0)) +
      geom_segment(mapping = aes(xend = lag, yend = 0)) +
      geom_hline(aes(yintercept = ciline), linetype = 3, color = set$colors[2]) +
      geom_hline(aes(yintercept = -ciline), linetype = 3, color = set$colors[2]) +
      labs(title = "Autocorrelation")
  }

  # add graphs
  if (!is.null(p4)) {
    pAll <- suppressWarnings(gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2, top = parName))
  } else {
    pAll <- suppressWarnings(gridExtra::grid.arrange(p1, p2, p3, nrow = 2, top = parName))
  }

  # save file
  if (!is.null(path)) {
    ggsave(filename = file.path(path, paste0(paste(prefix, "posterior_", sep = "_"), parName, ".", device)), plot = pAll, width = width, height = height)
  }

  # readline(prompt="Press [enter] to continue")
}

# -------------------------------------------------------------------------------------------

#' Plots the trend series and the (fitted) second observation equation and gives diagnostic
#' plots based on standardized residuals.
#'
#' @param tsl A list with two multiple time series objects for the first and second plot,
#'   respectively.
#' @param legend A list with two character vectors. The first contains the legend names for
#'   the first plot and so on.
#' @param title A list with the titles for the first three plots.
#' @param boundName The legend name of the confidence bounds.
#' @param res The residual series as time series. If \code{res = NULL}, all graphs related
#'   to the residual series will not be plotted.
#' @param namesPrint A character vector containing two names for the first two plots. The
#'   remaining names are creates automatically if \code{combine = FALSE}.
#' @inheritParams plot.NAWRUfit
#' @inheritParams plot.gap
#'
#' @importFrom stats start end window ts lag frequency time na.pass density acf qnorm reshape
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @keywords internal
plotGap <- function(tsl, legend, title, boundName, contribution, res = NULL, namesPrint, bounds,
                    combine, path, device, width, height) {

  # to suppress R CMD check notes
  potential <- gdp <- lb <- ub <- gap <- value <- variable <- NULL

  # settings
  setPrint <- list(width = width, height = height)
  set <- list()
  set$titleFontsize <- 10
  set$labelFontsize <- 8
  set$legendFontsize <- 8
  set$textSize <- 3.5
  set$colors_contr <- c("darkolivegreen4", "darkorange2", "deepskyblue3", "cyan4", "darkred", "darkgrey")
  set$colors_contr <- c("darkolivegreen4", "darkorange2", "deepskyblue3", "darkred", "darkgrey")
  set$colors <- c("black", "grey12", "grey24", "grey36", "grey48", "grey60")
  set$linetype <- c(2, 1, 3, 4, 5, 6)
  set$alpha <- 0.2
  set$freqYear <- 1 + floor(length(tsl[[1]][, 1]) / frequency(tsl[[1]][, 1]) / 15)
  set$freqYearSmall <- ceiling(set$freqYear * 2)
  if (!combine) {
    set$freqYearSmall <- set$freqYear
  }

  # ----- potential
  date <- zoo::as.Date(time(tsl[[1]]))
  df <- data.frame(date, tsl[[1]])

  if (contribution) {
    tmp <- na.trim(tsl[[1]])
    date <- zoo::as.Date(tmp)
    freq <- frequency(tsl[[1]])

    df <- t(tmp)
    colnames(df) <- as.character(time(tmp))
    df <- cbind(data.frame("variable" = rownames(df)), df)
    df <- reshape(df, direction = "long", varying = list(names(df)[-1]), idvar = "variable", times = date, v.names = "value")
    df$variable <- factor(df$variable, levels = colnames(tsl[[1]]))

    p0 <- ggplot(df) +
      geom_bar(aes(y = value, x = time, fill = variable), stat = "identity") +
      theme_classic() +
      coord_cartesian(xlim = c(date[1] - 365 - (365 + 30) / (freq * 2), date[length(date)] + (365 + 30) / (freq * 2)), expand = FALSE) +
      theme(
        plot.title = element_text(size = set$titleFontsize),
        axis.title = element_text(size = set$labelFontsize)
      ) +
      geom_line(data = data.frame("value" = apply(tmp, 1, sum), "time" = date), aes(y = value, x = time)) +
      scale_color_manual(name = NULL, values = "black") +
      scale_fill_manual(values = set$colors_contr)
    p0 <- p0 +
      scale_x_date(
        date_labels = "%Y", date_minor_breaks = "1 year",
        date_breaks = paste0(set$freqYearSmall, " year")
      ) +
      labs(title = title[[1]], x = "year", y = "") +
    theme(
      legend.position="bottom",
      legend.box = "horizontal",
      legend.box.margin=margin(-15,0,0,0),
      legend.justification="left",
      legend.text = element_text(size = set$legendFontsize),
      legend.key.size = unit(0.75, "lines"),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent"),
      legend.spacing.y = unit(0, "pt"),
      plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")
    ) +
      guides(color = guide_legend(order=1, ncol = 1, byrow = TRUE),
             # linetype = guide_legend(order=1, ncol = 1, byrow = TRUE),
             fill = guide_legend(order=2, ncol = 1)) 
    
  } else {
    color <- set$colors[1:(ncol(df) - 3)]
    color <- color[order(legend[[1]][1:(ncol(df) - 3)])]
    linetype <- set$linetype[1:(ncol(df) - 3)]
    linetype <- linetype[order(legend[[1]][1:(ncol(df) - 3)])]
    p0 <- ggplot(df, aes(x = date)) +
      theme_classic() +
      coord_cartesian(xlim = c(date[1] - 365, date[length(date)]), expand = FALSE) +
      theme(
        plot.title = element_text(size = set$titleFontsize),
        axis.title = element_text(size = set$labelFontsize),
        panel.grid.major.y = element_line(linetype = "solid")
      ) +
      geom_line(aes(y = potential, col = legend[[1]][1], linetype = legend[[1]][1])) +
      geom_line(aes(y = gdp, col = legend[[1]][2], linetype = legend[[1]][2]))
    if (bounds) {
      p0 <- p0 +
        geom_ribbon(aes(ymin = lb, ymax = ub, fill = boundName), alpha = 0.2) +
        scale_fill_manual(values = "grey12")
    }
    p0 <- p0 +
      scale_color_manual(name = NULL, values = color) +
      scale_linetype_manual(name = NULL, values = linetype) +
      scale_x_date(
        date_labels = "%Y", date_minor_breaks = "1 year",
        date_breaks = paste0(set$freqYearSmall, " year")
      ) +
      labs(title = title[[1]], x = "year", y = "") +
    theme(
      legend.position="bottom",
      legend.box = "horizontal",
      legend.box.margin=margin(-15,0,0,0),
      legend.justification="left",
      legend.text = element_text(size = set$legendFontsize),
      legend.key.size = unit(0.75, "lines"),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent"),
      legend.spacing.y = unit(0, "pt"),
      plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")
    ) +
      guides(color = guide_legend(order=1, ncol = 1, byrow = TRUE),
             linetype = guide_legend(order=1, ncol = 1, byrow = TRUE),
             fill = guide_legend(order=2)) 
  }

  # ----- gap
  if (contribution) {
    tmp <- na.trim(tsl[[2]])
    date <- zoo::as.Date(tmp)
    freq <- frequency(tsl[[2]])

    df <- t(tmp)
    colnames(df) <- as.character(time(tmp))
    df <- cbind(data.frame("variable" = rownames(df)), df)
    df <- reshape(df, direction = "long", varying = list(names(df)[-1]), idvar = "variable", times = date, v.names = "value")
    df$variable <- factor(df$variable, levels = colnames(tsl[[1]]))

    p1 <- ggplot(df) +
      geom_bar(aes(y = value, x = time, fill = variable), stat = "identity") +
      theme_classic() +
      coord_cartesian(xlim = c(date[1] - (365 + 30) / (freq * 2), date[length(date)] + (365 + 30) / (freq * 2)), expand = FALSE) +
      theme(
        plot.title = element_text(size = set$titleFontsize),
        axis.title = element_text(size = set$labelFontsize)
      ) +
      geom_line(data = data.frame("value" = apply(tmp, 1, sum), "time" = date), aes(y = value, x = time))
    p1 <- p1 +
      scale_color_manual(name = NULL, values = "black") +
      scale_fill_manual(values = set$colors_contr) +
      scale_x_date(
        date_labels = "%Y", date_minor_breaks = "1 year",
        date_breaks = paste0(set$freqYearSmall, " year")
      ) +
      labs(title = title[[2]], x = "year", y = "") +
    theme(
      legend.position="bottom",
      legend.box = "horizontal",
      legend.box.margin=margin(-15,0,0,0),
      legend.justification="left",
      legend.text = element_text(size = set$legendFontsize),
      legend.key.size = unit(0.75, "lines"),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent"),
      legend.spacing.y = unit(0, "pt"),
      plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")
    ) +
      guides(color = guide_legend(order=1, ncol = 1, byrow = TRUE),
             # linetype = guide_legend(order=1, ncol = 1, byrow = TRUE),
             fill = guide_legend(order=2, ncol = 1)) 
    
  } else {
    date <- zoo::as.Date(time(tsl[[2]]))
    df <- data.frame(date, tsl[[2]])
    color <- set$colors[3]
    color <- color[order(legend[[1]][1])]
    p1 <- ggplot(df, aes(x = date)) +
      theme_classic() +
      coord_cartesian(xlim = c(date[1], date[length(date)]), expand = FALSE) +
      theme(
        plot.title = element_text(size = set$titleFontsize),
        axis.title = element_text(size = set$labelFontsize),
        panel.grid.major.y = element_line(linetype = "solid")
      ) +
      geom_line(aes(y = gap, col = legend[[2]][1]))
    if (bounds) {
      p1 <- p1 +
        geom_ribbon(aes(ymin = lb, ymax = ub, fill = boundName), alpha = 0.2) +
        scale_fill_manual(values = "grey12")
    }
    p1 <- p1 +
      scale_color_manual(name = NULL, values = color) +
      scale_x_date(
        date_labels = "%Y", date_minor_breaks = "1 year",
        date_breaks = paste0(set$freqYearSmall, " year")
      ) +
      labs(title = title[[2]], x = "year", y = "") +
    theme(
      legend.position="bottom",
      legend.box = "horizontal",
      legend.box.margin=margin(-15,0,0,0),
      legend.justification="left",
      legend.text = element_text(size = set$legendFontsize),
      legend.key.size = unit(0.75, "lines"),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent"),
      legend.spacing.y = unit(0, "pt"),
      plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")
    ) +
      guides(color = guide_legend(order=1, ncol = 1, byrow = TRUE),
             linetype = guide_legend(order=1, ncol = 1, byrow = TRUE),
             fill = guide_legend(order=2)) 
    
  }

  if (combine) {
    pCombine <- suppressWarnings(gridExtra::grid.arrange(p0, p1, nrow = 1))
  } else {
    suppressWarnings(print(p0))
    suppressWarnings(print(p1))
  }

  # save files
  if (!is.null(path)) {
    if (combine) {
      do.call(ggsave, c(list(filename = file.path(path, paste0(namesPrint[1], "_", namesPrint[2], ".", device)), plot = pCombine), setPrint))
    } else {
      do.call(ggsave, c(list(filename = file.path(path, paste0(namesPrint[1], ".", device)), plot = p0), setPrint))
      do.call(ggsave, c(list(filename = file.path(path, paste0(namesPrint[2], ".", device)), plot = p1), setPrint))
    }
  }
}


# -------------------------------------------------------------------------------------------

#' Plots the trend series and the (fitted) second observation equation and gives diagnostic
#' plots based on standardized residuals.
#'
#' @param tsl A list with two multiple time series objects for the first and second plor,
#'   respectively.
#' @param legend A list with two character vectors. The first contains the legend names for
#'   the first plot and so on.
#' @param title A list with the titles for the first three plots.
#' @param n.ahead Integer indicating the forecast horizon.
#' @param boundName The legend name of the confidence bounds.
#' @param res The residual series as time series. If \code{res = NULL}, all graphs realted
#'   to the residual series will not be plotted.
#' @param namesPrint A character vector containing two names for the first two plots. The
#'   remaining names are creates automatically if \code{combine = FALSE}.
#' @inheritParams plot.NAWRUfit
#'
#' @importFrom stats start end window ts lag frequency time na.pass density acf qnorm
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @keywords internal
plotSSprediction <- function(tsl, legend, title, n.ahead, boundName, res = NULL, namesPrint, bounds,
                             combine, path, device, width, height) {
  
  # to suppress R CMD check notes.
  trend <- cycle <- orig <- anchor <- lb <- lb_fc <- ub <- ub_fc <- fitted <- E2 <- residuals <- ..density.. <- NULL
  
  # settings
  setPrint <- list(width = width, height = height)
  setPrintComb <- list(width = width, height = height * 2)
  set <- list()
  set$titleFontsize <- 10
  set$labelFontsize <- 8
  set$legendFontsize <- 8
  set$textSize <- 3.5
  set$colors <- c("black", "black", "black", "grey12", "grey24", "grey36", "grey48", "grey60")
  set$linetype <- c(2, 1, 3, 4, 5, 6)
  set$alpha <- 0.2
  set$freqYear <- 1 + floor(length(tsl[[1]][, 1]) / frequency(tsl[[1]][, 1]) / 15)
  set$freqYearSmall <- ceiling(set$freqYear * 2)
  if (!combine) {
    set$freqYearSmall <- set$freqYear
  }
  
  # ----- trend
    date <- zoo::as.Date(time(tsl[[1]]))
    xmin_fc_window <- date[length(date) - n.ahead + 1]
    tsl[[1]][1:(length(date) - n.ahead), c("lb_fc", "ub_fc")] <- NA
    df <- data.frame(date, tsl[[1]])
    color <- c(set$colors[1:2], rep(set$colors[7], 2))
    linetype <- c(set$linetype[1:2], rep(set$linetype[5], 2))
    p1 <- ggplot(df, aes(x = date)) +
      theme_classic() +
      coord_cartesian(xlim = c(date[1], date[length(date)]), expand = FALSE) +
      theme(
        plot.title = element_text(size = set$titleFontsize),
        axis.title = element_text(size = set$labelFontsize),
        panel.grid.major.y = element_line(linetype = "solid")
      ) +
      geom_line(aes(y = trend, color = legend[[1]][1], linetype = legend[[1]][1])) +
      geom_line(aes(y = orig, color = legend[[1]][2], linetype = legend[[1]][2])) +
      geom_line(aes(y = lb_fc, color = legend[[1]][3], linetype = legend[[1]][3])) +
      geom_line(aes(y = ub_fc, color = legend[[1]][4], linetype = legend[[1]][4])) +
      scale_color_manual(name = NULL, values = color, limits = legend[[1]]) +
      scale_linetype_manual(name = NULL, values = linetype, limits = legend[[1]]) +
      annotate("rect", xmin = xmin_fc_window, xmax = date[length(date)], 
               ymin = -Inf, ymax = Inf, alpha = .1)
    # add ribbon for state
    p1 <- p1 +
      geom_ribbon(aes(ymin = lb, ymax = ub, fill =  boundName), alpha = set$alpha) +
      scale_fill_manual(values = "grey12") 
    # color and legend settings
    p1 <- p1 +
      scale_x_date(
        date_labels = "%Y", date_minor_breaks = "1 year",
        date_breaks = paste0(set$freqYear, " year")
      ) +
      labs(title = title[[1]], x = "year", y = "") +
      theme(
        legend.position="bottom",
        legend.box = "horizontal",
        legend.box.margin=margin(-15,0,0,0),
        legend.justification="left",
        legend.text = element_text(size = set$legendFontsize),
        legend.key.size = unit(0.75, "lines"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.spacing.y = unit(0, "pt"),
        plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")
      ) +
      guides(color = guide_legend(order=1, ncol = 1, byrow = TRUE),
             linetype = guide_legend(order=1, ncol = 1, byrow = TRUE),
             fill = guide_legend(order=2)) 

  # ----- trend growth
  if (length(tsl) == 4) {
    date <- zoo::as.Date(time(tsl[[4]]))
    tsl[[4]][1:(length(date) - n.ahead), c("lb_fc", "ub_fc")] <- NA
    df <- data.frame(date, tsl[[4]])
    color <- c(set$colors[1:2], rep(set$colors[7], 2))
    linetype <- c(set$linetype[1:2], rep(set$linetype[5], 2))
    p4 <- ggplot(df, aes(x = date)) +
      theme_classic() +
      coord_cartesian(xlim = c(date[1], date[length(date)]), expand = FALSE) +
      theme(
        plot.title = element_text(size = set$titleFontsize),
        axis.title = element_text(size = set$labelFontsize),
        panel.grid.major.y = element_line(linetype = "solid")
      ) +
      geom_line(aes(y = trend, color = legend[[4]][1], linetype = legend[[4]][1])) +
      geom_line(aes(y = orig, color = legend[[4]][2], linetype = legend[[4]][2])) +
      geom_line(aes(y = lb_fc, color = legend[[4]][3], linetype = legend[[4]][3])) +
      geom_line(aes(y = ub_fc, color = legend[[4]][4], linetype = legend[[4]][4])) +
      scale_color_manual(name = NULL, values = color, limits = legend[[4]]) +
      scale_linetype_manual(name = NULL, values = linetype, limits = legend[[4]]) +
      annotate("rect", xmin = xmin_fc_window, xmax = date[length(date)], 
               ymin = -Inf, ymax = Inf, alpha = .1)
    # add ribbon for state
    p4 <- p4 +
      geom_ribbon(aes(ymin = lb, ymax = ub, fill =  boundName), alpha = set$alpha) +
      scale_fill_manual(values = "grey12") 
    # color and legend settings
    p4 <- p4 +
      scale_x_date(
        date_labels = "%Y", date_minor_breaks = "1 year",
        date_breaks = paste0(set$freqYear, " year")
      ) +
      labs(title = title[[4]], x = "year", y = "") +
      theme(
        legend.position="bottom",
        legend.box = "horizontal",
        legend.box.margin=margin(-15,0,0,0),
        legend.justification="left",
        legend.text = element_text(size = set$legendFontsize),
        legend.key.size = unit(0.75, "lines"),
        legend.title = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.spacing.y = unit(0, "pt"),
        plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")
      ) +
      guides(color = guide_legend(order=1, ncol = 1, byrow = TRUE),
             linetype = guide_legend(order=1, ncol = 1, byrow = TRUE),
             fill = guide_legend(order=2)) 
  }
  
  # ----- 2nd obs equation
  date <- zoo::as.Date(time(tsl[[2]]))
  tsl[[2]][1:(length(date) - n.ahead), c("lb_fc", "ub_fc")] <- NA
  df <- data.frame(date, tsl[[2]])
  color <- c(set$colors[1], rep(set$colors[7], 2))
  linetype <- c(set$linetype[2], rep(set$linetype[5], 2))
  p2 <- ggplot(df, aes(x = date)) +
    theme_classic() +
    coord_cartesian(xlim = c(date[1], date[length(date)]), expand = FALSE) +
    theme(
      plot.title = element_text(size = set$titleFontsize),
      axis.title = element_text(size = set$labelFontsize),
      panel.grid.major.y = element_line(linetype = "solid")
    ) +
    geom_line(aes(y = E2, color = legend[[2]][1], linetype = legend[[2]][1])) +
    geom_line(aes(y = lb_fc, color = legend[[2]][2], linetype = legend[[2]][2]), show.legend = FALSE) +
    geom_line(aes(y = ub_fc, color = legend[[2]][3], linetype = legend[[2]][3]), show.legend = FALSE) +
    annotate("rect", xmin = xmin_fc_window, xmax = date[length(date)], 
             ymin = -Inf, ymax = Inf, alpha = .1)
  # color and legend settings
  p2 <- p2 +
    scale_color_manual(name = NULL, values = color, limits = legend[[2]]) +
    scale_linetype_manual(name = NULL, values = linetype, limits = legend[[2]]) +
    scale_x_date(
      date_labels = "%Y", date_minor_breaks = "1 year",
      date_breaks = paste0(set$freqYear, " year")
    ) +
    labs(title = title[[2]], x = "year", y = "") +
    theme(
      legend.position="bottom",
      legend.box = "horizontal",
      legend.box.margin=margin(-15,0,0,0),
      legend.justification="left",
      legend.text = element_text(size = set$legendFontsize),
      legend.key.size = unit(0.75, "lines"),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent"),
      legend.spacing.y = unit(0, "pt"),
      plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")
    ) +
    guides(color = guide_legend(order=1, ncol = 1, byrow = TRUE),
           linetype = guide_legend(order=1, ncol = 1, byrow = TRUE),
           fill = guide_legend(order=2)) 

  # ----- cycle
  date <- zoo::as.Date(time(tsl[[3]]))
  df <- data.frame(date, tsl[[3]])
  color <- set$colors[1]
  linetype <- set$linetype[2]
  p3 <- ggplot(df, aes(x = date)) +
    theme_classic() +
    coord_cartesian(xlim = c(date[1], date[length(date)]), expand = FALSE) +
    theme(
      plot.title = element_text(size = set$titleFontsize),
      axis.title = element_text(size = set$labelFontsize),
      panel.grid.major.y = element_line(linetype = "solid")
    ) +
    geom_line(aes(y = cycle, color = legend[[3]][1], linetype = legend[[3]][1])) 
  p3 <- p3 +
    annotate("rect", xmin = xmin_fc_window, xmax = date[length(date)], 
             ymin = -Inf, ymax = Inf, alpha = .1)
  # add ribbon for state
  p3 <- p3 +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill =  boundName), alpha = set$alpha) +
    scale_fill_manual(values = "grey12")
  # color and legend settings
  p3 <- p3 +
    scale_color_manual(name = NULL, values = color, limits = legend[[3]]) +
    scale_linetype_manual(name = NULL, values = linetype, limits = legend[[3]]) +
    scale_x_date(
      date_labels = "%Y", date_minor_breaks = "1 year",
      date_breaks = paste0(set$freqYear, " year")
    ) +
    labs(title = title[[3]], x = "year", y = "") +
    theme(
      legend.position="bottom",
      legend.box = "horizontal",
      legend.box.margin=margin(-15,0,0,0),
      legend.justification="left",
      legend.text = element_text(size = set$legendFontsize),
      legend.key.size = unit(0.75, "lines"),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent"),
      legend.spacing.y = unit(0, "pt"),
      plot.margin = unit(c(0.1, 0.5, 0.1, 0.1), "cm")
    ) +
    guides(color = guide_legend(order=1, ncol = 1, byrow = TRUE),
           linetype = guide_legend(order=1, ncol = 1, byrow = TRUE),
           fill = guide_legend(order=2)) 

  # print
  if (combine) {
    if (length(tsl)==4) {
      pCombine <- suppressWarnings(gridExtra::grid.arrange(p1, p4, p2, p3, nrow = 2))
    } else {
      pCombine <- suppressWarnings(gridExtra::grid.arrange(p1, p2, p3, nrow = 2))
    }
  } else {
    suppressWarnings(print(p1))
    readline(prompt="Press [enter] to continue")
    if (length(tsl)==4) {
      suppressWarnings(print(p4))
      readline(prompt="Press [enter] to continue")
    }
    suppressWarnings(print(p2))
    readline(prompt="Press [enter] to continue")
    suppressWarnings(print(p3))
  }
  

  # save files
  if (!is.null(path)) {
    
    if (combine) {
      do.call(ggsave, c(list(filename = file.path(path, paste0(namesPrint[1], "_etc_prediction.", device)), plot = pCombine), setPrintComb))
    } else {
      do.call(ggsave, c(list(filename = file.path(path, paste0(namesPrint[1], ".", device)), plot = p1), setPrint))
      do.call(ggsave, c(list(filename = file.path(path, paste0(namesPrint[2], ".", device)), plot = p2), setPrint))
      do.call(ggsave, c(list(filename = file.path(path, paste0(namesPrint[3], ".", device)), plot = p3), setPrint))
      if (length(tsl)==4) {
        do.call(ggsave, c(list(filename = file.path(path, paste0(namesPrint[4], ".", device)), plot = p4), setPrint))
      }      
    }
  }
  
}

