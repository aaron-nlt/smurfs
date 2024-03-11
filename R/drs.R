# Helpers ======================================================================

# Format distance matrix to square form
# This is adapted from the package pracma, since one call to pracma::squareform
# is the only non-base dependency of smurfs. Since pracma::squareform itself
# relies only on base R, adapting it here to be used as an internal gives the
# smurfs package no dependencies, which is kinda nice! Furthermore, it will only
# be used as an internal function in this package -- users will not be able to
# use it directly -- so it's not a "double release" of the same function.
squareform = function(x){
  stopifnot(is.numeric(x) || is.complex(x))
  if(is.vector(x)){
    n = length(x)
    m = floor(sqrt(2*n))
    if(m*(m+1) != 2*n){
      stop("Argument 'x' does not correspond to a distance matrix.")
    }
    inds = c()
    k = m+1
    for(i in 1:(k-1)){
      inds = c(inds, (1+i+(i-1)*k):(i*k))
    }
    y = numeric(k*k)
    y[inds] = x
    y = matrix(y, k, k) + t(matrix(y, k, k))
  } else if(is.matrix(x)){
    m = nrow(x)
    n = ncol(x)
    if(m != n){
      stop("Argument 'x' must be a vector or a square matrix.")
    }
    if(any(diag(x) != 0)){
      stop("Argument 'x' can only have 0s on the diagonal.")
    }
    if(n == 1){
      return(c())
    }
    inds = c()
    for(i in 1:(n-1)){
      inds = c(inds, (1+i+(i-1)*n):(i*n))
    }
    y = x[inds]
  } else{
    stop("Argument 'x' must be a numeric or complex vector or matrix.")
  }
  return(y)
}

# In the remainder of this section, the function descriptions are copied from
# the original DRS Python code. These functions are internally called by the DRS
# function, and are never used as stand-alones by the user.

# A Python implementation of the Standard Dirichlet. This is faster than using
# scipy.dirichlet for single values
scaled_dirichlet = function(n, u){
  if(n <= 0){
    stop("n must be > 0")
  }
  if(n == 1){
    return(u)
  } else{
    intermediate = replicate(n, -log(1-runif(n=1, min=0, max=1))) #i think this is equivalent to rgamma(n,1,1)?
    divisor = sum(intermediate)
    return(intermediate * u / divisor)
  }
}

# This computes the Cayley-Megner Matrix determinant, and then normalises its
# sign based on what the normal simplex volume calculation would do. It does
# not normalise the the values as this calculation tends to break floating
# points for n > 102. The resulting value can be used to determine which of two
# n-dimensional simplicies is bigger, but little else.
cm_matrix_det_ns = function(vertices){
  square_dists = as.vector(dist(vertices) ^ 2)
  number_of_vertices = nrow(vertices)
  bordered_values = c(rep(1, number_of_vertices), square_dists)
  distance_matrix = squareform(bordered_values)
  detr = det(distance_matrix)
  if(length(vertices) %% 2 == 1){
    detr = -detr
  }
  if(detr <= 0){
    stop("Degenerate or invalid simplex")
  }
  return(detr)
}

# Takes an N-length vector, embeds it into the N+1 transformation vector,
# multiplies by the N+1xN+1 transformation matrix (which represents an N
# dimensional transformation), and returns the unembedded N-length vector.
embed_mult_debed = function(mat, vec){
  emb = mat %*% c(vec, 1)
  emb = emb[1:(length(emb)-1)]
  return(emb)
}

# Converts constraints into the coodinates of a constraints simplex
cts = function(limits){
  n = length(limits)
  simplex_coords = matrix(0, nrow=n, ncol=n)
  for(index in 1:n){
    simplex_coords[index,] = limits
    simplex_coords[index,index] = 0
    simplex_coords[index,index] = 1 - sum(simplex_coords[index,])
  }
  return(simplex_coords)
}

# Given the coordinates of a simplex, constructs a matrix which rescales the
# given simplex to the standard simplex through via translate-scale-translate
rmss = function(coords){
  sz = nrow(coords)
  coords_2 = matrix(1, nrow=sz+1, ncol=sz+1)
  coords_2[1:sz,1:sz] = coords
  translate_matrix = diag(1, nrow=sz+1, ncol=sz+1)
  translate_matrix[1:sz, sz+1] = -coords[1,]
  scale_matrix = diag(1, nrow=sz+1, ncol=sz+1)
  for(n in 1:sz){
    if(n == 1){
      scale_matrix[n,n] = -1/((translate_matrix %*% coords_2[2,])[1])
    } else{
      scale_matrix[n,n] = 1/((translate_matrix %*% coords_2[n,])[n])
    }
  }
  translate_matrix2 = diag(1, nrow=sz+1, ncol=sz+1)
  translate_matrix2[1,ncol(translate_matrix2)] = 1
  return(translate_matrix2 %*% scale_matrix %*% translate_matrix)
}

# Condition used in power scaling
power_condition = function(mat, vec){
  vec_t = embed_mult_debed(mat, vec)
  return(all(vec_t > 0))
}

# Given a set of limits for a simplex, works out the maximum power to which the
# corresponding rescale matrix can be raised such that M^P.V lies within the
# standard simplex, then returns this point.
power_scale = function(limits, vec){
  mat = rmss(cts(limits))
  cache = list("1" = mat)
  power = 1
  while(power_condition(mat, vec)){
    cache[[as.character(power)]] = mat
    power = power * 2
    mat = mat %*% mat
  }
  # We have the least power of 2 such that we rescale outside of allowed area
  power = power %/% 2
  mat = cache[[as.character(power)]]  # Greatest power of 2 inside allowed area
  power = power %/% 2
  while(power >= 1){  # Get remaining powers of 2
    mat2 = mat %*% cache[[as.character(power)]]
    if(power_condition(mat2, vec)){
      mat = mat2
    }
    power = power %/% 2
  }
  return(embed_mult_debed(mat, vec))
}

# Returns the volume of the sz-dimensional standard simplex
standard_simplex_vol = function(sz){
  result = cm_matrix_det_ns(diag(1, nrow=sz, ncol=sz))
  if(is.infinite(result)){
    stop(paste0("Cannot compute volume of standard ", sz, "-simplex"))
  }
  return(result)
}

# Performs the rescale operation given limits and an initial starting
# coordinate. The variable max_rescales determines the maximum number of
# rescale operations allowed before declaring failure; This should be at least
# the number of dimensions. epsilon is a parameter used to detect if floating
# point errors have caused the algorithm to diverge from its target.
rescale = function(limits, coord, max_rescales=1000, epsilon=1e-4){
  count = 0
  sz = length(limits)
  base_limits = rep(0, sz)
  for(count in 1:max_rescales){
    overlimits = coord > limits
    if(all(overlimits == 0)){
      return(list(count, coord))
    }
    bounding_limits = rep(NA, sz)
    w_true = which(overlimits == TRUE)
    w_false = which(overlimits == FALSE)
    bounding_limits[w_true] = limits[w_true]
    bounding_limits[w_false] = base_limits[w_false]
    coord = power_scale(bounding_limits, coord)
    divergence = 1 - sum(coord)
    if(abs(divergence) > epsilon){
      return(list(Inf, NULL))  # Diverged
    }
  }
  return(list(count, NULL))
}

# Implementation of the smallest simplex first rescale
ssr = function(limits, coord){
  limits_simplex = cts(limits)
  limits_simplex_vol = tryCatch({
    cm_matrix_det_ns(limits_simplex)
  },
  error=function(e){0}
  )
  if(limits_simplex_vol < standard_simplex_vol(length(limits))){
    mat = rmss(limits_simplex)
    new_limits = embed_mult_debed(mat, rep(0,length(limits)))
    inv_matrix = rmss(cts(new_limits))
    rs = rescale(new_limits, coord)
    niterations = rs[[1]]
    result = rs[[2]]
    if(is.null(result)){
      return(rs)
    }
    true_result = embed_mult_debed(inv_matrix, result)
    return(list(niterations, true_result))
  } else{
    return(rescale(limits, coord))
  }
}

# Dirichlet-Rescale ============================================================

#' Dirichlet-Rescale sampling
#'
#' Dirichlet-Rescale (DRS) sampling produces n values that sum to U, with each
#' value allowed its own domain: \eqn{x_i \in [a_i,b_i], \forall i}. This can be
#' seen as a generalization of the scaled flat Dirichlet, UUniFast\[-Discard\]
#' and Randfixedsum approaches, and can be used as a general-purpose replacement
#' for these algorithms. You can find a basic description of DRS and
#' comparison to previous algorithms
#' \link[https://sigbed.org/2020/12/21/the-dirichlet-rescale-drs-algorithm-a-general-purpose-method-underpinning-synthetic-task-set-generation/]{here}.
#'
#' The algorithm was originally proposed in:
#'
#' D. Griffin, I. Bate and R. I. Davis. "Generating Utilization Vectors for the
#' Systematic Evaluation of Schedulability Tests," 2020 IEEE Real-Time Systems
#' Symposium (RTSS), Houston, TX, USA, 2020, pp. 76-88.
#'
#' R code is adapted from the author's Python implementation, available publicly
#' on \link[https://github.com/dgdguk/drs]{GitHub}. Code used for adaptation was
#' acquired on March 5, 2024. At this time, GitHub edit records showed a last
#' commit date of November 7, 2020, and last commit message of "update to
#' version 2.0".
#'
#' @param n Integer, number of values in the sample.
#' @param sumu Numeric, sum of the sample. The default is 1.
#' @param lower_bounds Vector of numerics, lower bound for each value. Must have
#' length=n if provided. The default is NULL, which corresponds to a lower bound
#' of 0 for each value.
#' @param upper_bounds Vector of numerics, upper bound for each value. Must have
#' length=n if provided. The default is NULL, which corresponds to an upper
#' bound of sumu for each value.
#' @param float_tolerance Numeric, tolerance for floating point error. The
#' default is 1e-10.
#' @param drs_retries Integer, number of attempts to find a point. Default is
#' 1000.
#' @export
drs = function(n, sumu=1, lower_bounds=NULL, upper_bounds=NULL,
               float_tolerance=1e-10, drs_retries=1000){
  if(n == 0){
    if(!is.null(upper_bounds) & length(upper_bounds) != 0){
      stop("length(upper_bounds) must be equal to n")
    }
    if(!is.null(lower_bounds) & length(lower_bounds) != 0){
      stop("length(lower_bounds) must be equal to n")
    }
    return(
      list(
        output_point = c(),
        initial_point = c(),
        rescales_required = 0,
        retries = 0
      )
    )
  }
  if(n == 1){
    if(!is.null(upper_bounds)){
      if(length(upper_bounds) != 1){
        stop("length(upper_bounds) must be equal to n")
      }
      if(sumu - upper_bounds[1] > float_tolerance){
        stop("Upper bounds must sum to more than sumu")
      }
    }
    if(!is.null(lower_bounds)){
      if(length(lower_bounds) != 1){
        stop("length(lower_bounds) must be equal to n")
      }
      if(lower_bounds[1] - sumu > float_tolerance){
        stop("Lower bounds bounds must sum to less than max utilization")
      }
      if(!is.null(upper_bounds) & upper_bounds[1] == lower_bounds[1]){
        # sumu is prone to floating point error (if multiple bounds removed)
        # but if upper_bounds == lower_bounds, we know the true value
        sumu = upper_bounds[1]
      }
      return(
        list(
          output_point = c(sumu),
          initial_point = c(sumu),
          rescales_required = 0,
          retries = 0
        )
      )
    }
  }
  if(is.null(upper_bounds)){
    if(is.null(lower_bounds)){
      result = scaled_dirichlet(n, sumu)
    } else if(abs(sum(lower_bounds) - sumu) < float_tolerance){
      return(
        list(
          output_point = lower_bounds,
          initial_point = lower_bounds,
          rescales_required = 0,
          retries = 0
        )
      )
    } else{
      transformed_result = scaled_dirichlet(n, sumu - sum(lower_bounds))
      result = transformed_result + lower_bounds
    }
    return(
      output_point = result,
      initial_point = result,
      rescales_required = 0,
      retries = 0
    )
  }
  if(!is.null(lower_bounds)){
    if(sum(lower_bounds) - sumu > float_tolerance){ #original DRS code has this wrong
      stop("Lower bounds bounds must sum to less than max utilization")
    }
    for(index in 1:n){
      upper_bound = upper_bounds[index]
      lower_bound = lower_bounds[index]
      if(lower_bound > upper_bound){
        stop(
          paste0(
            "Lower bound > Upper bound (",
            lower_bound, " > ", upper_bound, ")"
          )
        )
      }
      if(lower_bound == upper_bound){
        fixed_point = upper_bound
        upper_bounds = upper_bounds[-index]
        lower_bounds = lower_bounds[-index]
        drs_result = drs(
          n = n - 1,
          sumu = sumu - lower_bound,
          upper_bounds = upper_bounds,
          lower_bounds = lower_bounds
        )
        drs_result$initial_point = append(
          x = drs_result$initial_point,
          values = fixed_point,
          after = index-1
        )
        drs_result$output_point = append(
          x = drs_result$output_point,
          values = fixed_point,
          after = index-1
        )
        return(drs_result)
      }
    }
    transformed_upper_bounds = upper_bounds - lower_bounds
    if(any(transformed_upper_bounds < 0)){
      stop("Bug detected with input, please report")
    }
    transformed_problem_result = drs(
      n = n,
      sumu = sumu - sum(lower_bounds),
      upper_bounds = transformed_upper_bounds
    )
    return(
      list(
        output_point = transformed_problem_result$output_point + lower_bounds,
        initial_point = transformed_problem_result$initial_point + lower_bounds,
        rescales_required = transformed_problem_result$rescales_required,
        retries = transformed_problem_result$retries
      )
    )
  }
  if(n != length(upper_bounds)){
    stop(
      paste0(
        "n=", n,
        "but utilisation constraints has length ", length(upper_bounds)
      )
    )
  }
  if(abs(sum(upper_bounds) - sumu) < float_tolerance){
    return(
      list(
        output_point = upper_bounds,
        initial_point = upper_bounds,
        rescales_required = 0,
        retries = 0
      )
    )
  }
  if(sum(upper_bounds) < sumu){
    stop("Upper bounds must sum to more than sumu")
  }
  limits = pmin(rep(1,length(upper_bounds)), upper_bounds/sumu)
  for(count in 1:drs_retries){
    initial_point = scaled_dirichlet(n, 1)
    ssr_result = ssr(limits, initial_point)
    iterations_required = ssr_result[[1]]
    final_unit_point = ssr_result[[2]]
    if(!is.null(final_unit_point)){
      break
    }
  }
  if(is.null(final_unit_point)){
    stop(
      paste0(
        "In ", drs_retries,
        " attempts, DRS failed to find a point that converged ",
        "before floating point error crept in"
      )
    )
  }
  final_point = final_unit_point * sumu
  return(
    list(
      output_point = final_point,
      initial_point = initial_point,
      rescales_required = iterations_required,
      retries = count
    )
  )
}

