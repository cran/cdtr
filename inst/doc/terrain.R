## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6,
                      fig.height = 4.5, dpi = 100)

## ----setup--------------------------------------------------------------------
library(cdtr)

## ----data---------------------------------------------------------------------
d <- readRDS(system.file("extdata", "tas_wellington.rds", package = "cdtr"))
contours <- d$contours
head(contours)
length(unique(contours$line))

## ----pslg---------------------------------------------------------------------
n <- nrow(contours)
same_line <- contours$line[-n] == contours$line[-1]
s0 <- which(same_line)
s1 <- s0 + 1L
c(vertices = n, segments = length(s0))

## ----mesh---------------------------------------------------------------------
m <- cdt_triangulate_attr(contours$x, contours$y, s0, s1,
                          PA = cbind(z = contours$z), erase = "hull")
c(vertices = nrow(m$P), triangles = nrow(m$T), crossings = m$n_crossing)

## ----plot-mesh----------------------------------------------------------------
plot_surface <- function(m, main = "") {
  zt <- rowMeans(matrix(m$PA[t(m$T), "z"], ncol = 3, byrow = TRUE))
  pal <- hcl.colors(100, "Terrain 2")
  col <- pal[cut(zt, 100, labels = FALSE)]
  plot(m$P, type = "n", asp = 1, axes = FALSE, xlab = "", ylab = "", main = main)
  polygon(m$P[t(cbind(m$T, NA)), 1], m$P[t(cbind(m$T, NA)), 2],
          col = col, border = NA)
}
op <- par(mar = c(0, 0, 1.5, 0))
plot_surface(m, "triangles shaded by mean height")
par(op)

## ----grid---------------------------------------------------------------------
g <- d$grid
nr <- nrow(g$z); nc <- ncol(g$z)
cell <- 40
xq <- rep(g$xmin + (seq_len(nc) - 0.5) * cell, each = nr)
yq <- rep(g$ymax - (seq_len(nr) - 0.5) * cell, times = nc)
zq <- as.vector(g$z)

inside <- !is.na(cdt_locate(m$P, m$T, xq, yq))
mean(inside)

surface_error <- function(m) {
  zz <- cdt_interpolate(m$P, m$T, m$PA, xq[inside], yq[inside])[, "z"]
  err <- zz - zq[inside]
  c(rmse = sqrt(mean(err^2)), mean_error = mean(err), max_abs = max(abs(err)))
}
round(surface_error(m), 2)

## ----quality------------------------------------------------------------------
summary(cdt_quality(m)[c("min_angle", "radius_edge")])

## ----relate-------------------------------------------------------------------
tri <- cdt_locate(m$P, m$T, xq, yq)
zz <- cdt_interpolate(m$P, m$T, m$PA, xq, yq)[, "z"]
err <- zz - zq
q <- cdt_quality(m)
angle_class <- cut(q$min_angle, c(0, 10, 20, 30, 60),
                   labels = c("<10", "10-20", "20-30", ">30"))
size_class <- cut(q$area, quantile(q$area, 0:4 / 4), include.lowest = TRUE,
                  labels = c("smallest", "2", "3", "largest"))
at_points <- data.frame(err = err, angle = angle_class[tri], size = size_class[tri])
at_points <- at_points[!is.na(tri), ]
rmse <- function(e) sqrt(mean(e^2))
by_angle <- tapply(at_points$err, at_points$angle, rmse)
by_size <- tapply(at_points$err, at_points$size, rmse)
round(rbind(share_of_points = table(at_points$angle) / nrow(at_points),
            rmse = by_angle), 3)

## ----relate-plot, fig.width = 7, fig.height = 7-------------------------------
emat <- matrix(err, nr, nc)
xs <- g$xmin + (seq_len(nc) - 0.5) * cell
ys <- g$ymax - (seq_len(nr) - 0.5) * cell
lim <- 30
op <- par(mfrow = c(2, 2), mar = c(0.5, 0.5, 2, 0.5))
image(xs, rev(ys), pmax(pmin(t(emat)[, nr:1], lim), -lim), asp = 1, axes = FALSE,
      xlab = "", ylab = "", zlim = c(-lim, lim), col = hcl.colors(21, "Blue-Red 3"),
      main = "error (m), mesh minus DEM")
for (l in unique(contours$line)) {
  i <- contours$line == l
  lines(contours$x[i], contours$y[i], col = "#00000040", lwd = 0.5)
}
shade <- c("firebrick", "darkorange", "gold", "grey90")
plot(m$P, type = "n", asp = 1, axes = FALSE, xlab = "", ylab = "",
     main = "triangles by minimum angle")
polygon(m$P[t(cbind(m$T, NA)), 1], m$P[t(cbind(m$T, NA)), 2],
        col = shade[angle_class], border = NA)
legend("topright", legend = levels(angle_class), fill = shade, title = "degrees",
       cex = 0.7, bty = "n")
par(mar = c(4, 4, 2, 1))
barplot(by_angle, xlab = "minimum angle (degrees)", ylab = "RMSE (m)",
        main = "error by triangle shape")
barplot(by_size, xlab = "triangle area (quartile)", ylab = "RMSE (m)",
        main = "error by triangle size")
par(op)

## ----refine-------------------------------------------------------------------
mr <- cdt_triangulate_attr(contours$x, contours$y, s0, s1,
                           PA = cbind(z = contours$z), erase = "hull",
                           max_area = 2000, min_angle = 25)
qa <- cdt_quality(m)
qb <- cdt_quality(mr)
rbind(
  contours_only = c(triangles = nrow(m$T), median_min_angle = median(qa$min_angle),
                    share_below_20 = mean(qa$min_angle < 20), surface_error(m)["rmse"]),
  refined = c(triangles = nrow(mr$T), median_min_angle = median(qb$min_angle),
              share_below_20 = mean(qb$min_angle < 20), surface_error(mr)["rmse"])
) |> round(3)

## ----plot-refined-------------------------------------------------------------
op <- par(mar = c(0, 0, 1.5, 0))
plot_surface(mr, "refined: area and minimum angle")
par(op)

## ----unrefined----------------------------------------------------------------
mr$unrefined[, c("criterion", "circumcenterOutside", "shortEdgeTriangles",
                 "sharpFixedCorner")]
qr <- cdt_quality(mr)
sum(qr$min_angle < 5)

