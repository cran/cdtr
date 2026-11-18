## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6,
                      fig.height = 3.2, dpi = 100)

## ----setup--------------------------------------------------------------------
library(cdtr)

## ----input--------------------------------------------------------------------
x <- c(0, 1, 1, 0, 0.25, 0.75, 0.75, 0.25)
y <- c(0, 0, 1, 1, 0.25, 0.25, 0.75, 0.75)
s0 <- 1:8
s1 <- c(2, 3, 4, 1, 6, 7, 8, 5)

r <- cdt_triangulate(x, y, s0, s1)
names(r)

## ----helper-------------------------------------------------------------------
plot_mesh <- function(m, ...) {
  e <- cdt_edges(m)
  plot(m$P, asp = 1, axes = FALSE, xlab = "", ylab = "", pch = ".", ...)
  segments(m$P[e[, 1], 1], m$P[e[, 1], 2], m$P[e[, 2], 1], m$P[e[, 2], 2],
           col = "grey50")
  segments(m$P[m$S[, 1], 1], m$P[m$S[, 1], 2], m$P[m$S[, 2], 1], m$P[m$S[, 2], 2],
           lwd = 2)
}
op <- par(mar = c(0.5, 0.5, 0.5, 0.5))
plot_mesh(r)
par(op)

## ----depth--------------------------------------------------------------------
table(r$depth)

## ----erase--------------------------------------------------------------------
sapply(c("hull", "outer", "holes"), function(e) {
  nrow(cdt_triangulate(x, y, s0, s1, erase = e)$T)
})

## ----holes--------------------------------------------------------------------
rh <- cdt_triangulate(x, y, s0, s1, holes = cbind(0.5, 0.5))
table(rh$depth)

## ----refine-------------------------------------------------------------------
r2 <- cdt_triangulate(x, y, s0, s1, max_area = 0.01, min_angle = 25,
                      erase = "holes")
table(r2$origin)
op <- par(mar = c(0.5, 0.5, 0.5, 0.5))
plot_mesh(r2)
par(op)

## ----unrefined----------------------------------------------------------------
r2$min_edge_length
r2$unrefined

## ----cross--------------------------------------------------------------------
xc <- c(0, 1, 0, 1)
yc <- c(0, 1, 1, 0)
rc <- cdt_triangulate(xc, yc, c(1, 3), c(2, 4))
rc$origin

## ----sources------------------------------------------------------------------
xs <- c(0, 2, 2, 0, 1, 3, 3, 1)
ys <- c(0, 0, 2, 2, 1, 1, 3, 3)
rs <- cdt_triangulate(xs, ys, 1:8, c(2, 3, 4, 1, 6, 7, 8, 5),
                      edge_sources = TRUE)
head(rs$edge_sources)

## ----locate-------------------------------------------------------------------
q <- cbind(c(0.1, 0.5, 0.9, 2), c(0.1, 0.5, 0.9, 2))
cdt_locate(r2$P, r2$T, q[, 1], q[, 2])

## ----attr---------------------------------------------------------------------
ra <- cdt_triangulate_attr(x, y, s0, s1, PA = cbind(z = 10 * x),
                           max_area = 0.02, erase = "holes")
all.equal(ra$PA[, "z"], 10 * ra$P[, 1])

