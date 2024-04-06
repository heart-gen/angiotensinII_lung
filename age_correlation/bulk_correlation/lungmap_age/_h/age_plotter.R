## This is the ageplotter function adapted for caudate BSP3.

require(rafalib)
require(RColorBrewer)

agePlotter_v2 <- function(y, age, mod = matrix(rep(1, length(y)), ncol = 1),
                          mainText, smoothIt = TRUE, jitter = TRUE, ageLabel = "bottom",
                          orderByAge = TRUE, ylim = NULL, ageBreaks = c(-1, 0, 1, 10, 100),
                          ylab = "Adjusted Expression", pointColor = 2, lineColor = 1,
                          alreadyFitted = NULL, fetal_present = TRUE, ...)
{
    stopifnot(length(ageBreaks) >= 4)
    stopifnot(ageLabel %in% c("bottom", "top"))
    stopifnot(length(lineColor) == length(unique(pointColor)))

    if (orderByAge) {
        oo <- order(age, decreasing = FALSE)
        y <- y[oo]
        age <- age[oo]
        mod <- mod[oo, , drop = FALSE]
        if (length(pointColor) == length(age)){pointColor <- pointColor[oo]}
        if (!is.null(alreadyFitted)){alreadyFitted <- alreadyFitted[oo]}
    }

    nfits <- length(unique(pointColor))
    if (is.null(alreadyFitted)) {
        if (nfits == 1) {
            fit <- fitted(lm(y ~ mod - 1))
        } else {
            fit <- lapply(unique(pointColor), function(col) {
                fitted(lm(y[pointColor == col] ~ mod[pointColor == col, , drop = FALSE] - 1))
            })
        }
    } else {
        if (nfits > 1) {warning('Only one line will be draw. Try again specifying "mod"')}
        nfits <- 1
        fit <- alreadyFitted
    }
                                        # Make breaks
    idx <- rafalib::splitit(cut(age, breaks = ageBreaks, lab = FALSE))
    make_line <- function(i, case0 = FALSE) {
        if (nfits == 1) {
            if (case0) {
                lines(age[idx[[i]]][-1], fit[idx[[i]]][-1], col = lineColor, lwd = 6)
            } else {
                lines(age[idx[[i]]], fit[idx[[i]]], col = lineColor, lwd = 6)
            }
        } else {
            for (j in seq_len(nfits)) {
                nfit_split <- rafalib::splitit(cut(age[pointColor == unique(pointColor)[j]],
                    breaks = ageBreaks, lab = FALSE
                ))
                if (case0) {
                    lines(age[pointColor == unique(pointColor)[j]][nfit_split[[i]]][-1],
                        fit[[j]][nfit_split[[i]]][-1],
                        col = lineColor[j], lwd = 6
                    )
                } else {
                    lines(age[pointColor == unique(pointColor)[j]][nfit_split[[i]]], fit[[j]][nfit_split[[i]]],
                        col = lineColor[j], lwd = 6
                    )
                }
            }
        }
    }

    nBreaks <- length(ageBreaks) - 1
    layout(matrix(rep(seq_len(nBreaks), c(5, rep(2, nBreaks - 2), 4)),
        nrow = 1, byrow = TRUE
    ))
    palette(RColorBrewer::brewer.pal(8, "Set1"))
    par(mar = c(4, 5, 3, 0.45))
    if (is.null(ylim)) ylims <- range(y, na.rm = TRUE) else ylims <- ylim
    if (jitter) xx <- jitter(age, amount = 0.005) else xx <- age
                                        # Make fetal tissue optional parameter
    if(fetal_present){
        plot(y ~ xx, subset = idx[[1]],
             main = "", ylab = ylab, xlab = "",
             ylim = ylims, cex.axis = 1.4, cex.lab = 1.75,
             pch = 21, cex = 1.4, bg = pointColor, xaxt = "n",
             xlim = c(range(age[idx[[1]]]) + c(-0.01, 0.01)), ...)
        if (smoothIt) {make_line(1)}
        fetal_rang <- round(range(age[idx[[1]]]) * 52 + 40, 0)
        fetal_ax <- seq(fetal_rang[1], fetal_rang[2], round(diff(fetal_rang) / 4, 0))
        axis(1, at = (fetal_ax - 40) / 52, labels = fetal_ax, 1, cex.axis = 1.5)
        if (ageLabel == "bottom") {
            text(x = quantile(age[idx[[1]]], 0.33), y = min(ylims), "PCW",cex = 1.5)
        } else if (ageLabel == "top") {
            text(x = quantile(age[idx[[1]]], 0.33), y = max(ylims), "PCW",cex = 1.5)
        }
                                        # infant + child
        par(mar = c(4, 0.25, 3, 0.25))
        for (j in 2:(nBreaks - 1)) {
            plot(y ~ age, subset = idx[[j]],
                 main = "", ylab = "", xlab = "", yaxt = "n", cex = 1.4,
                 xlim = range(age[idx[[j]]]) + c(-0.02, 0.02),
                 ylim = ylims, cex.axis = 1.4, pch = 21, bg = pointColor, ...)
        if (ageBreaks[2] == 0 & smoothIt) {make_line(j)}
        if (ageBreaks[2] < 0 & smoothIt) {make_line(j, case0 = TRUE)}
        }
                                        # adults
        par(mar = c(4, 0.25, 3, 1))
        plot(y ~ age, subset = idx[[length(idx)]],
             main = "", ylab = "", xlab = "", yaxt = "n", cex = 1.4,
             xlim = range(age[idx[[length(idx)]]]) + c(-0.01, 0.01),
             ylim = ylims, cex.axis = 1.4, pch = 21, bg = pointColor, ...)
        if (smoothIt) {make_line(length(idx))}
    } else {
                                        # infant + child
        par(mar = c(4, 0.25, 3, 0.25))
        for (j in 1:(nBreaks - 1)) {
            if(j == 1){
                par(mar=c(4,5,3,0.45))
                plot(y ~ age, subset = idx[[j]], cex.lab = 1.75,
                 main = "", ylab = ylab, xlab = "", cex = 1.4,
                 xlim = range(age[idx[[j]]]) + c(-0.01, 0.01),
                 ylim = ylims, cex.axis = 1.4, pch = 21, bg = pointColor, ...)
            } else {
                par(mar=c(4,0.25,3,0.25))
                plot(y ~ age, subset = idx[[j]],
                 main = "", ylab = "", xlab = "", yaxt = "n", cex = 1.4,
                 xlim = range(age[idx[[j]]]) + c(-0.02, 0.02),
                 ylim = ylims, cex.axis = 1.4, pch = 21, bg = pointColor, ...)
            }
        if (ageBreaks[1] == 0 & smoothIt) {make_line(j)}
        if (ageBreaks[1] < 0 & smoothIt) {make_line(j, case0 = TRUE)}
        }
        # adults
        par(mar = c(4, 0.25, 3, 1))
        plot(y ~ age, subset = idx[[length(idx)]],
             main = "", ylab = "", xlab = "", yaxt = "n", cex = 1.4,
             xlim = range(age[idx[[length(idx)]]]) + c(-0.01, 0.01),
             ylim = ylims, cex.axis = 1.4, pch = 21, bg = pointColor, ...)
        if (smoothIt) {make_line(length(idx))}
    }
    mtext(mainText, outer = TRUE, line = -2.5, cex = 1.5)
    mtext("Age", side = 1, outer = TRUE, line = -1.25, cex = 1.35)
    ##mtext("Residualized Expression", side = 2, outer = TRUE, line=0, cex=1.5)
}

agePlotter_v3 <- function(y, age, mod = matrix(rep(1, length(y)), ncol = 1),
                          mainText, smoothIt = TRUE, jitter = TRUE, ageLabel = "bottom",
                          orderByAge = TRUE, ylim = NULL, ageBreaks = c(-1, 0, 1, 10, 100),
                          ylab = "Adjusted Expression", pointColor = 2, lineColor = 1,
                          alreadyFitted = NULL, fetal_present = TRUE, ...)
{
    stopifnot(length(ageBreaks) >= 3)
    stopifnot(ageLabel %in% c("bottom", "top"))
    stopifnot(length(lineColor) == length(unique(pointColor)))

    if (orderByAge) {
        oo <- order(age, decreasing = FALSE)
        y <- y[oo]
        age <- age[oo]
        mod <- mod[oo, , drop = FALSE]
        if (length(pointColor) == length(age)){pointColor <- pointColor[oo]}
        if (!is.null(alreadyFitted)){alreadyFitted <- alreadyFitted[oo]}
    }

    nfits <- length(unique(pointColor))
    if (is.null(alreadyFitted)) {
        if (nfits == 1) {
            fit <- fitted(lm(y ~ mod - 1))
        } else {
            fit <- lapply(unique(pointColor), function(col) {
                fitted(lm(y[pointColor == col] ~ mod[pointColor == col, , drop = FALSE] - 1))
            })
        }
    } else {
        if (nfits > 1) {warning('Only one line will be draw. Try again specifying "mod"')}
        nfits <- 1
        fit <- alreadyFitted
    }
                                        # Make breaks
    idx <- rafalib::splitit(cut(age, breaks = ageBreaks, lab = FALSE))
    make_line <- function(i, case0 = FALSE) {
        if (nfits == 1) {
            if (case0) {
                lines(age[idx[[i]]][-1], fit[idx[[i]]][-1], col = lineColor, lwd = 6)
            } else {
                lines(age[idx[[i]]], fit[idx[[i]]], col = lineColor, lwd = 6)
            }
        } else {
            for (j in seq_len(nfits)) {
                nfit_split <- rafalib::splitit(cut(age[pointColor == unique(pointColor)[j]],
                    breaks = ageBreaks, lab = FALSE
                ))
                if (case0) {
                    lines(age[pointColor == unique(pointColor)[j]][nfit_split[[i]]][-1],
                        fit[[j]][nfit_split[[i]]][-1],
                        col = lineColor[j], lwd = 6
                    )
                } else {
                    lines(age[pointColor == unique(pointColor)[j]][nfit_split[[i]]], fit[[j]][nfit_split[[i]]],
                        col = lineColor[j], lwd = 6
                    )
                }
            }
        }
    }

    nBreaks <- length(ageBreaks) - 1
    layout(matrix(rep(seq_len(nBreaks), c(5, rep(2, nBreaks - 2), 4)),
        nrow = 1, byrow = TRUE
    ))
    palette(RColorBrewer::brewer.pal(8, "Set1"))
    par(mar = c(4, 5, 3, 0.45))
    if (is.null(ylim)) ylims <- range(y, na.rm = TRUE) else ylims <- ylim
    if (jitter) xx <- jitter(age, amount = 0.005) else xx <- age
                                        # Make fetal tissue optional parameter
    if(fetal_present){
        plot(y ~ xx, subset = idx[[1]],
             main = "", ylab = ylab, xlab = "",
             ylim = ylims, cex.axis = 1.4, cex.lab = 1.75,
             pch = 21, cex = 1.4, bg = pointColor, xaxt = "n",
             xlim = c(range(age[idx[[1]]]) + c(-0.01, 0.01)), ...)
        if (smoothIt) {make_line(1)}
        fetal_rang <- round(range(age[idx[[1]]]) * 52 + 40, 0)
        fetal_ax <- seq(fetal_rang[1], fetal_rang[2], round(diff(fetal_rang) / 4, 0))
        axis(1, at = (fetal_ax - 40) / 52, labels = fetal_ax, 1, cex.axis = 1.5)
        if (ageLabel == "bottom") {
            text(x = quantile(age[idx[[1]]], 0.33), y = min(ylims), "PCW",cex = 1.5)
        } else if (ageLabel == "top") {
            text(x = quantile(age[idx[[1]]], 0.33), y = max(ylims), "PCW",cex = 1.5)
        }
                                        # infant + child
        par(mar = c(4, 0.25, 3, 0.25))
        for (j in 2:(nBreaks - 1)) {
            plot(y ~ age, subset = idx[[j]],
                 main = "", ylab = "", xlab = "", yaxt = "n", cex = 1.4,
                 xlim = range(age[idx[[j]]]) + c(-0.02, 0.02),
                 ylim = ylims, cex.axis = 1.4, pch = 21, bg = pointColor, ...)
        if (ageBreaks[2] == 0 & smoothIt) {make_line(j)}
        if (ageBreaks[2] < 0 & smoothIt) {make_line(j, case0 = TRUE)}
        }
                                        # adults
        par(mar = c(4, 0.25, 3, 1))
        plot(y ~ age, subset = idx[[length(idx)]],
             main = "", ylab = "", xlab = "", yaxt = "n", cex = 1.4,
             xlim = range(age[idx[[length(idx)]]]) + c(-0.01, 0.01),
             ylim = ylims, cex.axis = 1.4, pch = 21, bg = pointColor, ...)
        if (smoothIt) {make_line(length(idx))}
    } else {
                                        # infant + child
        par(mar = c(4, 0.25, 3, 0.25))
        for (j in 1:(nBreaks - 1)) {
            if(j == 1){
                par(mar=c(4,5,3,0.45))
                plot(y ~ age, subset = idx[[j]], cex.lab = 1.75,
                 main = "", ylab = ylab, xlab = "", cex = 1.4,
                 xlim = range(age[idx[[j]]]) + c(-0.01, 0.01),
                 ylim = ylims, cex.axis = 1.4, pch = 21, bg = pointColor, ...)
            } else {
                par(mar=c(4,0.25,3,0.25))
                plot(y ~ age, subset = idx[[j]],
                 main = "", ylab = "", xlab = "", yaxt = "n", cex = 1.4,
                 xlim = range(age[idx[[j]]]) + c(-0.02, 0.02),
                 ylim = ylims, cex.axis = 1.4, pch = 21, bg = pointColor, ...)
            }
        if (ageBreaks[1] == 0 & smoothIt) {make_line(j)}
        if (ageBreaks[1] < 0 & smoothIt) {make_line(j, case0 = TRUE)}
        }
        # adults
        par(mar = c(4, 0.25, 3, 1))
        plot(y ~ age, subset = idx[[length(idx)]],
             main = "", ylab = "", xlab = "", yaxt = "n", cex = 1.4,
             xlim = range(age[idx[[length(idx)]]]) + c(-0.01, 0.01),
             ylim = ylims, cex.axis = 1.4, pch = 21, bg = pointColor, ...)
        if (smoothIt) {make_line(length(idx))}
    }
    mtext(mainText, outer = TRUE, line = -2.5, cex = 1.5)
    mtext("Age", side = 1, outer = TRUE, line = -1.25, cex = 1.35)
    ##mtext("Residualized Expression", side = 2, outer = TRUE, line=0, cex=1.5)
}
