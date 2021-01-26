# AUTHORS: Yulia Rubanova and Nil Sahin
# Modified for package trackSig by Cait Harrigan

## Round vector of number to percentages. This function is from CRAN package MESS,
## replicated in full here to avoid unecessary dependancy for single internal funciton.
##
## Rounds a vector of numeric values to percentages ensuring that they add up to 100%
##
## Returns a vector of numeric values.
##
## @param x A numeric vector with non-negative values.
## @param decimals An integer giving the number of decimals that are used
## @param ties A string that is either 'random' (the default) or 'last'. Determines how to break ties. Random is random, last prefers to break ties at the last position
## @return Returns a numeric vector of the same length as x
## @author Claus Ekstrom \email{claus@@rprimer.dk}
##
round_percent <- function(x, decimals=0L, ties=c("random", "last")) {

  ties <- match.arg(ties)

  ## Do a few sanity checks
  if (!is.numeric(x)) { stop("only works on numeric vectors") }
  if (min(x)<0) { stop("only works on non-negative vectors") }
  if (decimals<0) { stop("number of decimals should be a non-negative integer") }
  decimals <- as.integer(decimals)

  multiplier <- 10^(2+decimals)

  x <- x/sum(x)*multiplier  # Standardize result
  res <- floor(x)           # Find integer bits
  rsum <- sum(res)          # Find out how much we are missing
  if(rsum<multiplier) {
    ## Distribute points based on remainders and a random tie breaker
    tiebreaker <- switch(ties,
                         random = sample(length(x)),
                         last = seq(length(x))
    )
    o <- order(x%%1, tiebreaker, decreasing=TRUE)
    res[o[1:(multiplier-rsum)]] <- res[o[1:(multiplier-rsum)]]+1
  }
  res/(10^decimals)
}


#' Display a histogram of mutation phi and signature exposure for a signature trajectory
#'
#' @param trajectory a list containing named elements "mixtures", "changepoints",
#'   and "binData". See @seealso \link{TrackSig}.
#' @param trajPlot a ggplot object generated using plotTrajectory @seealso [plotTrajectory()]
#' @param truncateStrategy strategy to handle a truncated axis (for example, for a trajectory plotted with anmac = F). \cr
#' One of "exclude", in which mutations lying off the axis are not displayed, or "stack, in which mutaitons lying off the axis are added to the first histogram bin. \cr
#' It is advised to use "exclude", unless there are very few clonal mutations.
#' @return tableGrob
#'
#' @name addPhiHist
#' @export

# TODO: create phiHist function as standalone, and add two plots with cowplot.
addPhiHist <- function(trajectory, trajPlot, truncateStrategy = c("exclude", "stack")){
  # TODO: phiHist plot should be able to stack or exclude >1 ccf if x range is truncated.
  # create phi histogram ggplot and add it on top of cpPlot

  vcaf <- trajectory$binData
  truncateStrategy <- match.arg(truncateStrategy)

  # input checking
  assertthat::assert_that((trajPlot$plot_env$linearX == FALSE),
                          msg = "Can't add phi histogram to linearly scaled axis\n")


  sigs <- base::unique(trajPlot$data$Signatures)
  inBin <- stats::aggregate(vcaf, by = list(vcaf$bin), FUN = length)$bin

  trajPlot$data$bin <- rep((dim(trajPlot$data)[1]/length(sigs)):1, each = length(sigs))

  # truncate vcaf if necessary (anmac = F may recude bins)
  lowerLim = max(vcaf$bin) - max(trajPlot$data$bin)
  if(truncateStrategy == "exclude"){

    vcaf <- vcaf[vcaf$bin > lowerLim,]
    vcaf$bin <- vcaf$bin - lowerLim

  } else if (truncateStrategy == "stack"){

    vcaf$bin <- vcaf$bin - lowerLim
    vcaf$phi[vcaf$bin <= 0] <- vcaf$phi[vcaf$bin == 1][[1]]
    vcaf$bin[vcaf$bin <= 0] <- 1

  }

  vcaf$sigAssignment <- factor(NA, levels = sigs)

  # display check
  if(any(trajPlot$data$exposure * inBin[trajPlot$data$bin] < 1)){
    warning("Signatures with activity level less than 1% will not be displayed in the phi histogram\n")
  }

  for (bin in vcaf$bin){

    # get signatures for each mutation
    freq = trajPlot$data$exposure[trajPlot$data$bin == bin] * inBin[bin]
    sigFreq <- rep(sigs, times = round_percent(freq))

    # sample order for pretty histograms
    sigFreq <- sample(sigFreq)

    # set signature assignment
    vcaf$sigAssignment[vcaf$bin == bin] <- sigFreq

  }

  # plot stacked phi histogram
  phiHist <- ( ggplot2::ggplot(vcaf, ggplot2::aes(x = .data$phi, fill = .data$sigAssignment))
               + ggplot2::geom_histogram(binwidth = 0.02, position = "stack")
               + ggplot2::scale_fill_hue(limits=levels(vcaf$sigAssignment))
               + ggplot2::xlab("")
               + ggplot2::ylab("")
               + ggplot2::theme_bw()
               + ggplot2::theme_void()
               + ggplot2::theme(legend.position="none")
               + ggplot2::scale_x_reverse()

  )


  # insert the histogram
  plotHat <- cowplot::insert_xaxis_grob(trajPlot, phiHist, position = "top", height = ggplot2::unit(0.15, "null"))

  grid::grid.draw(plotHat)

  return(plotHat)
}

#' Plot the exolutionary trajectory of a tumour
#'
#' For each bin in a set of signature mixtures, the mixture is plotted accross
#' pseudo-time. Provided changepoints will be highlighted.
#'
#'
#' @param trajectory a list containing named elements "mixtures", "changepoints",
#'   and "binData". See @seealso \link{TrackSig}.
#' @param linearX logical whether to plot with a linearly spaced x-axis grid, or
#'   with binned phi values
#' @param anmac logical whether to plot x-axis restricted to ccf space, or use
#'   estimated average number of mutant alleles per cell (anmac)
#' @param show logical whether to print the plot as a side effect when called
#' @return ggplot object
#'
#' @name plotTrajectory
#' @import rlang
#' @export

plotTrajectory <- function(trajectory, linearX = F, anmac = F, show = F){

  if(!is.null(trajectory)){
    mixtures <- trajectory[["mixtures"]]
    changepoints <- trajectory[["changepoints"]]
    binData <- trajectory[["binData"]]
  }

  # input checking
  assertthat::assert_that(!is.null(mixtures), msg = "Could not find mixtures for timeline, please supply through results or mixtures paramter.\n")

  # set the phis to colnames(mixtures) - note: used when anmac = T
  phis <- as.numeric(colnames(mixtures))

  # mixtures and phis are binned the same way
  assertthat::assert_that(length(phis) == dim(mixtures)[2],
                          msg = "The mixtures object is mal-specified. Column names should correspond to binned phis.\n")

  # phis should be decreasings
  assertthat::assert_that(all(order(phis, decreasing = T) == 1:length(phis)),
                          msg = "The mixtures object is mal-specified. Binned phis (column names) should be in decreasing order.\n")

  if(!anmac){ # take x-axis as ccf scale

    # ccf is min(1, anmac)
    # truncate x-axis at phi = 1
    truncateSel <- which(phis <= 1)
    phis <- phis[truncateSel]
    mixtures <- mixtures[,truncateSel,drop = FALSE]

    # change x-axis lable
    xAx <- "Cancer cell fraction"

    # adjust changepoint indexing
    if (!is.null(changepoints)){
      changepoints <- which(truncateSel %in% changepoints)
    }

  }else{ xAx <- "Average number of mutant alleles per cell" }

  # Plotting the change of mutational signature weights during evolution specified as the order of phi
  colnames(mixtures) <- dim(mixtures)[2]:1
  timeline <- reshape2::melt(mixtures)
  colnames(timeline) <- c("Signatures", "xBin", "exposure")
  timeline$xBin <- as.numeric(timeline$xBin)
  timeline$exposure <- as.numeric(timeline$exposure)


  if(!linearX){ # ggplot formatting specific for non-linear scale

    # non-linear scale shows ccf densities

    # place labels in a way that depends on bin density
    # take 8 times the smallest spacing (%)
    spacing <- 800 * min(c(NA, phis) - c(phis, NA), na.rm = T)

    ticSel <- seq(1, length(phis), by = spacing)
    ticLab <- rep("", length(phis))
    ticLab[ticSel] <- round(phis, 2)[ticSel]

    # increasing phi by bin
    timeline$xBin <- phis[timeline$xBin]
    timeline$xBin <- timeline$xBin[length(timeline$xBin) : 1]

    g <- (  ggplot2::ggplot(data = timeline)
          + ggplot2::geom_vline(xintercept = phis, alpha = 0.3)
          + ggplot2::aes(x = .data$xBin, y = .data$exposure * 100,
                         group = .data$Signatures, color = .data$Signatures)
          + ggplot2::scale_x_reverse(breaks = phis, labels = ticLab)
         )

    # slice changepoints (reverse axis means max to min)
    cpPos <- base::cbind(phis[changepoints], phis[changepoints + 1])

  }else{ # ggplot formatting specific for linear scale

    ticSel <- seq(1, length(phis), length.out = min(length(phis), 25))
    ticLab <- rep("", length(phis))
    ticLab[ticSel] <- round(phis, 2)[ticSel]

    g <- (  ggplot2::ggplot(data = timeline)
          + ggplot2::geom_vline(xintercept = 0:(length(phis) + 1), alpha = 0.3)
          + ggplot2::aes(x = .data$xBin, y = .data$exposure * 100, group = .data$Signatures, color = .data$Signatures)
          + ggplot2::scale_x_reverse(breaks = length(phis):1, labels = ticLab)
    )

    # slice changepoints (reverse axis means max to min)
    cpPos <- base::cbind((length(phis):1)[changepoints], (length(phis):1)[changepoints + 1])

  }

  # TODO: have truncate x range as option
  # TODO: adjust text element size, and alpha for repear lines

  # general ggplot formatting
  g <- (   g
           + ggplot2::geom_point()
           + ggplot2::geom_line()
           + ggplot2::theme_bw()
           + ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                            panel.grid.minor.x = ggplot2::element_blank())
           + ggplot2::ylab("Signature Exposure (%)")
           + ggplot2::xlab(xAx)
  )

  # add changepoints
  if (!is.null(changepoints)) {

    for (i in 1:dim(cpPos)[1]) {
      g <- g + ggplot2::annotate("rect", xmax=cpPos[i,1], xmin=cpPos[i,2],
                                 ymin=-Inf, ymax=Inf, alpha=0.3, fill = "black")

    }
  }

  if (show){print(g)}

  return(g)
}


#TODO: ggplot version of this funciton that can be passed to addPhiHist
plotChangepointChoice <- function(trajectory){

  nMut <- dim(trajectory$binData)[1]
  nBin <- trajectory$binData$bin[nMut]
  binSize <- sum(trajectory$binData$bin == 1)

  potentialCps <- binSize * 1:(nMut/binSize)

  graphics::plot(trajectory$binData$phi, ylab = "Empirical Phi", xlab = "Mutation Index",
       main = "Potential changepoints in black, selected changepoints in red")

  graphics::abline(v = potentialCps)

  if(!is.null(trajectory$changepoints)){
    graphics::abline(v = potentialCps[trajectory$changepoints], col = 2)
  }

}




