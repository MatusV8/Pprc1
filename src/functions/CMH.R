# Haldane-Anscombe correction: add 0.5 to every cell of a 2x2 table if any cell is 0
haldane <- function(m) if (any(m == 0)) m + 0.5 else m

# Mantel-Haenszel OR with Robins-Breslow-Greenland (RBG) SE of log(OR_MH)
# Works for k >= 1; pseudocount applied per stratum only when a cell is 0
mh_rbg <- function(arr) {
  for (k in seq_len(dim(arr)[3])) arr[, , k] <- haldane(arr[, , k])
  a <- arr[1, 1, ]; b <- arr[1, 2, ]; c <- arr[2, 1, ]; d <- arr[2, 2, ]
  n <- a + b + c + d
  P <- (a + d) / n
  Q <- (b + c) / n
  R <- (a * d) / n
  S <- (b * c) / n
  sR <- sum(R); sS <- sum(S)
  var_logOR <- sum(P * R) / (2 * sR^2) +
    sum(P * S + Q * R) / (2 * sR * sS) +
    sum(Q * S) / (2 * sS^2)
  list(logOR = log(sR / sS), SE = sqrt(var_logOR))
}

# Build a 2 x 2 x K array from data.table group columns
build_array <- function(GFPlow_count, GFPlow_total, Input_count, Input_total) {
  K <- length(GFPlow_count)
  arr <- array(NA_real_, dim = c(2, 2, K),
               dimnames = list(c("in_gene", "not_in_gene"),
                               c("GFPlow", "Input"),
                               NULL))
  for (k in seq_len(K)) {
    arr[, , k] <- matrix(c(as.numeric(GFPlow_count[k]),
                           as.numeric(GFPlow_total[k] - GFPlow_count[k]),
                           as.numeric(Input_count[k]),
                           as.numeric(Input_total[k]  - Input_count[k])),
                         nrow = 2)
  }
  arr
}
