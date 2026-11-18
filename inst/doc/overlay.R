## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6,
                      fig.height = 4.5, dpi = 100)

## ----setup--------------------------------------------------------------------
library(cdtr)

## ----inputs-------------------------------------------------------------------
rings <- list(
  list(feature = 1, xy = rbind(c(0, 0), c(4, 0), c(4, 4), c(0, 4))),
  list(feature = 1, xy = rbind(c(1, 1), c(1, 3), c(3, 3), c(3, 1))),
  list(feature = 2, xy = rbind(c(2.5, 2.5), c(5.5, 2.5), c(5.5, 5.5), c(2.5, 5.5))),
  list(feature = 3, xy = rbind(c(1, -0.5), c(3.5, 5), c(-0.5, 3))),
  list(feature = 4, xy = rbind(c(-0.8, 0.5), c(0, 0.5), c(0, 3.5), c(-0.8, 3.5))),
  ## a frame around everything, feature 0, which is not an input to be counted
  list(feature = 0, xy = rbind(c(-1, -1), c(7, -1), c(7, 7), c(-1, 7)))
)

## ----pslg---------------------------------------------------------------------
nv <- vapply(rings, function(r) nrow(r$xy), 1L)
first <- cumsum(c(0L, nv[-length(nv)]))
xy <- do.call(rbind, lapply(rings, `[[`, "xy"))
s0 <- unlist(lapply(seq_along(nv), function(k) first[k] + seq_len(nv[k])))
s1 <- unlist(lapply(seq_along(nv), function(k) first[k] + c(seq_len(nv[k])[-1], 1L)))
owner <- rep(vapply(rings, `[[`, 1, "feature"), nv)
length(s0)

## ----mesh---------------------------------------------------------------------
r <- cdt_triangulate(xy[, 1], xy[, 2], s0, s1, edge_sources = TRUE)
table(r$origin)
c(fixed_edges = nrow(r$S), triangles = nrow(r$T))

## ----sources------------------------------------------------------------------
head(r$edge_sources)
table(n_inputs = table(r$edge_sources$edge))

## ----shared-------------------------------------------------------------------
multi <- as.integer(names(which(table(r$edge_sources$edge) > 1)))
r$edge_sources[r$edge_sources$edge %in% multi, ] |>
  transform(feature = owner[segment])

## ----walk---------------------------------------------------------------------
w <- cdt_winding(r, owner, outside = c(-0.9, -0.9))
anyNA(w)
head(w)

## ----check--------------------------------------------------------------------
centre <- function(r) {
  cbind(rowMeans(matrix(r$P[t(r$T), 1], ncol = 3, byrow = TRUE)),
        rowMeans(matrix(r$P[t(r$T), 2], ncol = 3, byrow = TRUE)))
}
inside_ring <- function(p, ring) {
  j <- nrow(ring); inside <- rep(FALSE, nrow(p))
  for (i in seq_len(nrow(ring))) {
    cross <- ((ring[i, 2] > p[, 2]) != (ring[j, 2] > p[, 2])) &
      (p[, 1] < (ring[j, 1] - ring[i, 1]) * (p[, 2] - ring[i, 2]) /
         (ring[j, 2] - ring[i, 2]) + ring[i, 1])
    inside <- xor(inside, cross); j <- i
  }
  inside
}
p <- centre(r)
truth <- sapply(1:4, function(f) {
  inside <- rep(FALSE, nrow(p))
  for (ring in rings) if (ring$feature == f) inside <- xor(inside, inside_ring(p, ring$xy))
  inside
})
all((w != 0) == truth)

## ----depth--------------------------------------------------------------------
table(depth = r$depth, covered_by = rowSums(w != 0))

## ----pieces-------------------------------------------------------------------
set <- apply(w != 0, 1, function(i) paste(which(i), collapse = "+"))
set[set == ""] <- "none"
area <- function(r) {
  a <- r$P[r$T[, 1], ]; b <- r$P[r$T[, 2], ]; c <- r$P[r$T[, 3], ]
  0.5 * abs((b[, 1] - a[, 1]) * (c[, 2] - a[, 2]) - (c[, 1] - a[, 1]) * (b[, 2] - a[, 2]))
}
sort(round(tapply(area(r), set, sum), 3), decreasing = TRUE)

## ----plot---------------------------------------------------------------------
sets <- sort(unique(set))
pal <- hcl.colors(length(sets), "Set 3", alpha = 0.8)
pal[sets == "none"] <- "grey92"
col <- pal[match(set, sets)]
tri <- t(cbind(r$T, NA))
op <- par(mar = c(0, 0, 0, 0))
plot(r$P, type = "n", asp = 1, axes = FALSE, xlab = "", ylab = "")
polygon(r$P[tri, 1], r$P[tri, 2], col = col, border = "grey80", lwd = 0.5)
segments(r$P[r$S[, 1], 1], r$P[r$S[, 1], 2], r$P[r$S[, 2], 1], r$P[r$S[, 2], 2], lwd = 1.5)
legend("topright", legend = sets, fill = pal,
       cex = 0.7, bty = "n", title = "covered by")
par(op)

