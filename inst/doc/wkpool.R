## ----include = FALSE----------------------------------------------------------
has_wkpool <- requireNamespace("wkpool", quietly = TRUE) &&
  requireNamespace("wk", quietly = TRUE) &&
  utils::packageVersion("wkpool") >= "0.5.0"
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6,
                      fig.height = 3.2, dpi = 100, eval = has_wkpool)

## ----setup--------------------------------------------------------------------
library(cdtr)
library(wkpool)

## ----polygons-----------------------------------------------------------------
x <- wk::wkt(c(
  "POLYGON ((0 0, 4 0, 4 4, 0 4, 0 0), (1 1, 3 1, 3 3, 1 3, 1 1))",
  "POLYGON ((3.5 3.5, 6 3.5, 6 6, 3.5 6, 3.5 3.5))"
))
pool <- merge_coincident(establish_topology(x))
pslg <- as_pslg(pool)
str(pslg)

## ----tri----------------------------------------------------------------------
to_tri <- function(w, ...) {
  r <- cdt_pslg(as_pslg(merge_coincident(establish_topology(w))),
                erase = "holes", ...)
  vertices <- data.frame(.vx = seq_len(nrow(r$P)), x = r$P[, 1], y = r$P[, 2])
  new_wkpool_tri(vertices, r$T[, 1], r$T[, 2], r$T[, 3])
}
a <- to_tri(x[1])
a
tri_measures(a)[c("area", "perimeter")]

## ----locate-------------------------------------------------------------------
pts <- cbind(c(0.5, 2, 3.5), c(0.5, 2, 3.5))
tri_locate(a, pts)

## ----overlay------------------------------------------------------------------
b <- to_tri(x[2])
ov <- tri_overlay(a, b)
tri_overlay_relate(ov)

## ----refine-------------------------------------------------------------------
a2 <- to_tri(x[1], max_area = 0.5, min_angle = 25)
length(a2)
tri_measures(a2)$area

