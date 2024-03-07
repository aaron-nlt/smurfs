#' Randfixedsum samping
#'
#' Randfixedsum sampling produces n values that sum to U, with each value on
#' \[a,b\]. This is an improvement over the UUniFast approaches in that it allows
#' for a non-zero lower bound, but still requires that the domain of each
#' individual value is constant.
#'
#' The algorithm was originally proposed in:
#'
#' R. Stafford. Random vectors with fixed sum. Technical Report Available at
#' https://www.mathworks.com/matlabcentral/fileexchange/9700-random-vectors-with-fixed-sum,
#' MathWorks, 2006.
#'
#' And formally published with exemplified uses in sampling in:
#'
#' P. Emberson, R. Stafford, and R. I. Davis. Techniques for the synthesis of
#' multiprocessor tasksets. In Proc. International Workshop on Analysis Tools
#' and Methodologies for Embedded and Real-time Systems (WATERS), pages 6–11,
#' July 2010.
#'
#' R code is adapted from the original 2006 MATLAB code.
#'
#' @param n Integer, number of values in the sample
#' @param m Integer, number of samples. The default is 1.
#' @param s Numeric, sum of each sample. The default is 1.
#' @param a Numeric, lower bound on the values. The default is 0.
#' @param b Numeric, upper bound on the values. The default is 1.
#' @export
randfixedsum = function(n, m=1, s=1, a=0, b=1){
  # Check the arguments.
  if((m %% 1 != 0) | (n %% 1 != 0) | (m < 0) | (n < 1)){
    stop("n must be a whole number and m a non-negative integer.")
  }
  if((s < n*a)| (s > n*b)| (a >= b)){
    stop("Inequalities n*a <= s <= n*b and a < b must hold.")
  }
  # Rescale to a unit cube: 0 <= x(i) <= 1
  s = (s-n*a)/(b-a)
  # Construct the transition probability table, t.
  # t(i,j) will be utilized only in the region where j <= i + 1.
  k = max(min(floor(s),n-1),0)     #Must have 0 <= k <= n-1
  s = max(min(s,k+1),k)            #Must have k <= s <= k+1
  s1 = s - (k:(k-n+1))             #s1 will never be negative
  s2 = ((k+n):(k+1)) - s           #s2 will never be negative
  w = matrix(0, nrow=n, ncol=n+1)
  w[1,2] = .Machine$double.xmax    #Scale for full 'double' range
  t = matrix(0, nrow=n-1, ncol=n)
  tiny = .Machine$double.eps       #The smallest positive 'double' no.
  for(i in 2:n){
    tmp1 = w[i-1, 2:(i+1)] * s1[1:i] / i
    tmp2 = w[i-1, 1:i] * s2[(n-i+1):n] / i
    w[i, 2:(i+1)] = tmp1 + tmp2
    tmp3 = w[i, 2:(i+1)] + tiny                 #In case tmp1 & tmp2 are both 0,
    tmp4 = as.numeric(s2[(n-i+1):n] > s1[1:i])  #then t is 0 on left & 1 on right
    t[i-1, 1:i] <- (tmp2/tmp3)*tmp4 + (1-tmp1/tmp3)*(1-tmp4)
  }
  # Derive the polytope volume v from the appropriate element in the bottom row
  # of w.
  v = n^(3/2) * (w[n, k+2]/.Machine$double.xmax) * (b-a)^(n-1L)
  # Now compute the matrix x.
  x = matrix(0, nrow=n, ncol=m)
  rt = matrix(runif((n-1)*m, 0, 1), nrow=n-1, ncol=m)  #For random selection of simplex type
  rs = matrix(runif((n-1)*m, 0, 1), nrow=n-1, ncol=m)  #For random location within a simplex
  s = matrix(s, nrow=1, ncol=m)
  j = matrix(k+1, nrow=1, ncol=m)                      #For indexing in the t table
  sm = matrix(0, nrow=1, ncol=m)                       #Start with sum zero
  pr = matrix(1, nrow=1, ncol=m)                       #Start with product 1
  for(i in (n-1L):1L) {                      #Work backwards in the t table
    e = as.numeric(rt[n-i,] <= t[i,j])       #Use rt to choose a transition
    sx = rs[n-i,]^(1/i)                      #Use rs to compute next simplex coordinate
    sm = sm + (1-sx)*pr*s/(i+1)              #Update sum
    pr = sx * pr                             #Update product
    x[n-i,] = sm + pr*e                      #Calculate x using simplex coordinates
    s = s - e                                #Transition adjustment
    j = j - e                                #Transition adjustment
  }
  x[n,] <- sm + pr*s                         # Compute the last x
  # Randomly permute the order in the columns of x and rescale.
  xm = apply(x, 2, function(col){
    sample(col, n, replace=FALSE)
  })
  xm = xm*(b-a) + a
  return(xm)
}
