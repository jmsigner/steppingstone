#' @export
#' @useDynLib steppingstone

controll <- function() {
  list(
    boundary = "wrap",
    on_error = "return_xy"
  )
}

#' @export
steppingstone <- function(alpha, omegas, resources, n,
                          rarify_by = 100,
                          burnin = round(n / 10), xy0 = c(0, 0),
                          with_attr = FALSE, more = NULL, random_error = 0,
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
    w <- walk_func(tm1, n, xy0, nc)$steps

  } else if (ctl$boundary == "stop") {
    w <- walk_func_boundary_stop(tm1, n, xy0, nc)
    if (w$outcome == -1) {
      if (ctl$on_boundary_error == "return_xy") {
        w <- w$steps[1:w$last_step]
        #w <- w$steps
      } else
        stop("Encountered boundary")
    } else {
      w <- w$steps
    }

  } else if (ctl$boundary == "reflective") {
    w <- walk_func_boundary_reflective(tm1, n, xy0, nc, ctl$max_try)$steps

  }

  # rm burnin
  if (burnin > 0)
    w <- w[-(1:burnin)]

  # rarify
  w <- w[seq(1, length(w), by = rarify_by)]


  xy <- cbind(x = w %% nc + 0.5, y = w %/% nc + 0.5)

  if (random_error != 0) {
    xy[, 1] <- xy[, 1] + runif(nrow(xy), -random_error, random_error)
    xy[, 2] <- xy[, 2] + runif(nrow(xy), -random_error, random_error)

  }

  xy <- if (with_attr) {
    list(xy = xy,
         alpha = alpha,
         omegas = omegas,
         n = n,
         xy0 = xy0,
         rarify_by = rarify_by,
         burnin = burnin,
         more = more,
         with_attr = TRUE)
  } else {
    list(xy = xy,
         with_attr = FALSE)
  }

  class(xy) <-  c("sstone_walk", "list")
  xy
}



ss_ud <- function(alpha, omegas, resources, n,
                            burnin = round(n / 10), xy0 = c(0, 0),
                            with_attr = FALSE, more = NULL) {
  r <- as.array(resources)
  nc <- dim(r)[2]

  #  r <- r[nc:1, , ] # bring in the right order
  r <- sapply(1:dim(r)[3], function(i) as.vector(t(r[, , i]))) # flatten 2d arrays to 1d
  dists <- c(0, 1, 1, 1, 1)

  # tranisition porbabilities
  tm1 <- tpm_func(alpha, omegas, r, dists, nc)

  # obtain UD
  w <- ud_func(tm1, n, xy0, nc, burnin)
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

  xy
}
