#' UniFast sampling
#'
#' UUniFast sampling produces n values that sum to U, with each value on \[0,U\].
#' In practice, this is effectively equivalent to the scaled flat Dirichlet
#' approach.
#'
#' The algorithm was originally proposed in:
#'
#' E. Bini and G. Buttazzo. Measuring the performance of schedulability tests.
#' Journal of Real-Time Systems, 30(1-2):129–154, 2005.
#'
#' R code is adapted from this paper's MATLAB code.
#'
#' @param n Integer, number of values in the sample.
#' @param U Numeric, sum of the sample. The default is 1.
#' @export
uunifast = function(n, U=1){
  sumU = U
  vectU = rep(NA, n)
  for(i in 1:(n-1)){
    nextSumU = sumU * runif(n=1, min=0, max=1)**(1/(n-i))
    vectU[i] = sumU - nextSumU
    sumU = nextSumU
  }
  vectU[n] = sumU
  return(vectU)
}

#' UUniFast-Discard sampling
#'
#' UUniFast-Discard sampling produces n values that sum to U, with each value on
#' \[0,1\]. This is an extension of the UUniFast algorithm for cases where
#' individual values in each sample cannot exceed 1; this is achieved by
#' discarding UUniFast samples where any value > 1. Because it uses rejection
#' sampling, it is increasingly slow as \eqn{U \rightarrow n}.
#'
#' The algorithm was originally proposed in:
#'
#' R. I. Davis and A. Burns. Improved priority assignment for global fixed
#' priority pre-emptive scheduling in multiprocessor real-time systems.
#' Real-Time Systems, 47:1–40, 2010.
#'
#' R code is a natural extension of the adapted UUniFast code.
#'
#' @param n Integer, number of values in the sample.
#' @param U Numeric, sum of the sample. The default is 1.
#' @param maxiter Integer, number of attempts to make before breaking. The
#' default is 1000.
#' @export
uunifast_discard = function(n, U=1, maxiter=1000){
  k = 1
  while(k <= maxiter){
    vectU = uunifast(n=n, U=U)
    # Keep if no utilization exceeds 1
    if(all(vectU <= 1)){
      return(vectU)
    }
    # Next iter
    k = k+1
  }
  stop(
    paste(
      "UUniFast-Discard failed to find a valid solve within",
      maxiter,
      "iterations"
    )
  )
}


