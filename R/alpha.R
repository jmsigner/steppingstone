#' Movement Probability to Alpha
#'
#' @param p Probability of moving.
#'
#' @return numeric vector.
#' @export
#'
#' @examples
#'
pm2alpha <- function(p) {
  log((-4 * (p - 1))/p)
}

