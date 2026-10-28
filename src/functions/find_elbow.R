find_elbow <- function(y) {
  .vec <- sort(y)
  x <- seq_along(.vec)
  line_vec <- c(x[length(.vec)] - x[1], .vec[length(.vec)] - .vec[1])
  
  d <- sapply(seq_along(.vec), function(i) {
    px <- x[i] - x[1]
    py <- .vec[i] - .vec[1]
    abs(line_vec[1] * py - line_vec[2] * px) / sqrt(sum(line_vec^2))
  })
  
  idx <- which.max(d)
  .vec[idx]
}