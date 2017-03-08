#' @export
ba <- function(x, y) {
  r1 <- x[]
  r2 <- y[]
  r1 <- r1 / sum(r1)
  r2 <- r2 / sum(r2)
  ## bhattacharyya's afinity
  sum(sqrt(r1 * r2))
}
