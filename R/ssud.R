#' @export
ssud <- function(alpha, omegas, resources, n,
               burnin = round(n / 10), xy0 = c(0, 0),
               with_attr = FALSE, more = NULL,
               ctl = controll()) {

  r <- raster::as.array(resources)
  nc <- dim(r)[2]

  r <- r[nc:1, , , drop = FALSE] # reorder
  r <- sapply(1:dim(r)[3], function(i) as.vector(t(r[, , i]))) # flatten 2d arrays to 1d, uses column major order
  dists <- c(0, 1, 1, 1, 1)

  # start must be an integer
  xy0 <- floor(xy0)

  # Calcualte the transition probabilities once
  tm1 <- tpm_func(alpha, omegas, r, dists, nc)

  if (any(is.nan(tm1))) {
    stop("Transition matrices contain NaN. Try to rescale resources.")
  }

  # simulate walk
  if (ctl$boundary == "wrap") {
    w <- ud_func(tm1, n, xy0, nc, burnin)

  } else {
    stop("boundary: only wrap is implemented")
  }
  ud <- resources[[1]]
  ud[] <- w$ud
  last <- c(x = w$lastPos %% nc + 0.5, y = w$lastPos %/% nc + 0.5)

  xy <- if (with_attr) {
    list(ud = ud,
         last_pos = last,
         alpha = alpha,
         omegas = omegas,
         resources = resources,
         n = n,
         xy0 = xy0,
         burnin = burnin,
         more = more,
         with_attr = TRUE)
  } else {
    list(ud = ud,
         last_pos = last,
         with_attr = FALSE)
  }
}
