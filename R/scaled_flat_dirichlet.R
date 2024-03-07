#' Scaled flat Dirichlet sampling
#'
#' Scaled flat Dirichlet sampling produces n values that sum to U, with each
#' value on \[0,U\]. This takes advantage of the flat Dirichlet distribution, a
#' a special case of the Dirichlet distribution that is a multivariate uniform
#' over the (n-1) simplex. Thus, sampling from it produces a vector of sum 1
#' where all values are on \[0,1\]. A scaling parameter U can be used to change
#' the sum of the vector, but this changes the domain of the values to \[0,U\].
#'
#' The Dirichlet distribution was originally proposed in:
#'
#' I. Olkin and H. Rubin. Multivariate beta distributions and independence
#' properties of the wishart distribution. Annals of Mathematical Statistics,
#' 35(1):261–269, March 1964.
#'
#' Plentiful methods already exist in R for sampling from the Dirichlet
#' distribution. This function takes advantage of the fact that, for i.i.d.
#' random variables \eqn{X_1 ~ Gamma(a_1, b), ..., X_k \sim Gamma(a_k, b)}, the
#' vector \eqn{(\frac{X_1}{\sum_{i=1}^k(X_i)}, ..., \frac{X_k}{\sum_{i=1}^k(X_i)}) \sim Dir(a_1, ..., a_k)}.
#'
#' @param n Int, number of values in the sample.
#' @param U Numeric, sum of the sample. The default is 1.
#' @export
scaled_flat_dirichlet = function(n, U=1){
  intermediate = rgamma(n=n, shape=1, rate=1)
  return(intermediate / sum(intermediate) * U)
}
