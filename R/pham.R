##
##  Implementation of the 2004 Pham et al. paper
##
##  Created by Daniel Rodríguez Pérez on 6/9/2014.
##
##  Copyright (c) 2014 Daniel Rodríguez Pérez.
##
##  This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>
##

#' Selection of K in K-means Clustering
#' 
#' Selection of k in k-means clustering based on Pham et al. paper
#' 
#' @param x numeric matrix of data, or an object that can be coerced to such a
#'        matrix.
#' @param max_centers maximum number of clusters for evaluation.
#' @param k_threshold maximum f(k) from which it can not consider the existence
#'        of more than one cluster in the data set. The default value is 0.85.
#' @param trace show a progress bar
#' @param ... arguments to be passed to the kmeans method
#' 
#' @return an object with the f(k) results
#' 
#' @details
#' This function implements the method proposed by Pham, Dimov and Nguyen for
#' selecting the number of clusters for the K-means algorithm. In this method
#' a function \eqn{f(K)} is used to evaluate the quality of the resulting
#' clustering and help decide on the optimal value of \eqn{K} for each data
#' set. The \eqn{f(K)} function is defined as
#' \deqn{f(K) = \left\{
#' \begin{array}{rl}
#'  1 & \mbox{if $K = 1$} \\
#'  \frac{S_K}{\alpha_K S_{K-1}} & \mbox{if $S_{K-1} \ne 0$, $\forall K >1$} \\
#'  1 & \mbox{if $S_{K-1} = 0$, $\forall K >1$}
#' \end{array} \right.}{f(K) =
#'  1, if K = 1;
#'  (S_K)/(\alpha_K S_{K-1}, if S_{K-1} \ne 0, forall K >1;
#'  1, if S_{K-1} = 0, forall K > 1}
#' where \eqn{S_K} is the sum of the distortion of all cluster and \eqn{\alpha_K}
#' is a weight factor which is defined as
#' \deqn{\alpha_K = \left\{
#' \begin{array}{rl}
#'  1 - \frac{3}{4 N_d}                        & \mbox{if $K = 1$ and $N_d > 1$} \\
#'  \alpha_{K-1} + \frac{1 - \alpha_{K-1}}{6}  & \mbox{if $K > 2$ and $N_d > 1$}
#' \end{array} \right.}{\alpha_K = 
#'  1 - 3/(4 * N_d), if K = 1 and N_d > 1;
#'  \alpha_{K-1} + (1 - \alpha_{K-1})/6, if K > 2 and N_d > 1}
#' where \eqn{N_d} is the number of dimensions of the data set.
#' 
#' In this definition \eqn{f(K)} is the ratio of the real distortion to the
#' estimated distortion and decreases when there are areas of concentration in
#' the data distribution.
#' 
#' The values of \eqn{K} that yield \eqn{f(K) < 0.85} can be recommended for
#' clustering. If there is not a value with \eqn{f(K) < 0.85}, it cannot be
#' considered the existence of clusters in the data set.
#' 
#' @examples
#' # Create a data set with two clusters
#' dat <- matrix(c(rnorm(100, 2, .1), rnorm(100, 3, .1),
#'                 rnorm(100, -2, .1), rnorm(100, -3, .1)), 200, 2)
#'
#' # Ejecute the method
#' sol <- kselection(dat)
#' 
#' # Get the results
#' k   <- num_clusters(sol) # optimal number of clustes
#' f_k <- get_f_k(sol)      # the f(k) vector
#' 
#' # Plot the results
#' plot(sol)
#' 
#' @author Daniel Rodriguez Perez
#' 
#' @references
#' D T Pham, S S Dimov, and C D Nguyen "Selection of k in k-means clustering",
#' Mechanical Engineering Science, 2004, pp. 103-119.
#' 
#' @seealso \code{\link{num_clusters}}, \code{\link{get_f_k}}
#' 
#' @import tools
#' 
#' @rdname kselection
#' @export kselection
kselection <- function(x,
                       max_centers = 15,
                       k_threshold = 0.85,
                       trace       = FALSE, ...) {
  if (max_centers < 2)
    stop("'max_centers' must be greater than 2")
    
  x <- as.matrix(x)
  if (!is.numeric(x)) 
    stop('x must contain numerical data')
  
  num_row <- dim(x)[1]
  num_col <- dim(x)[2]
  
  if (num_row < max_centers) {
    max_centers <- num_row
  } 
  
  f_k <- rep(1, max_centers)
  s_k <- rep(1, max_centers)
  a_k <- alpha_k(num_col, max_centers)
  
  if (trace) {
    pb <- txtProgressBar(min   = 0,
                         max   = max_centers,
                         style = 3)
  }
  
  for (k in 1:max_centers) {
    mod_info <- kmeans(x, centers = k, ...)
    s_k[k]   <- sum(mod_info$withinss)
    
    if (k == 1) {
      f_k[k] <- 1  
    } else {
      if (s_k[k - 1] == 0) {
        f_k[k] <- 1
      } else {
        f_k[k] <- s_k[k]  / (a_k[k] * s_k[k - 1])
      }      
    }

    if (trace) {
      setTxtProgressBar(pb, k)
    }
  }
  
  result <- list(k           = which_cluster(f_k, k_threshold),
                 f_k         = f_k,
                 max_centers = max_centers,
                 k_threshold = k_threshold)
  class(result) <- 'Kselection'
  
  return(result)
}

#' Get the f(k) vector
#' 
#' Get the f(k) vector
#' 
#' @param obj the output of kselection function
#' 
#' @return the vector of f(k) function
#' 
#' @examples
#' # Create a data set with two clusters
#' dat <- matrix(c(rnorm(100, 2, .1), rnorm(100, 3, .1),
#'                 rnorm(100, -2, .1), rnorm(100, -3, .1)), 200, 2)
#'               
#' # Get the f(k) vector
#' sol <- kselection(dat)
#' f_k <- get_f_k(sol)
#' 
#' @seealso \code{\link{num_clusters}},
#' 
#' @rdname get_f_k
#' @export get_f_k
get_f_k <- function(obj) {
  UseMethod("get_f_k")
}

#' @method get_f_k default
#' @export
get_f_k.default <- function(obj) {
  NULL
}

#' @method get_f_k Kselection
#' @export
get_f_k.Kselection <- function(obj) {
  obj$f_k
}

#' Get the number of clusters
#' 
#' Get the number of clusters
#' 
#' @param obj the output of kselection function
#' 
#' @return the number of clusters
#' 
#' @examples
#' # Create a data set with two clusters
#' dat <- matrix(c(rnorm(100, 2, .1), rnorm(100, 3, .1),
#'                 rnorm(100, -2, .1), rnorm(100, -3, .1)), 200, 2)
#'               
#' # Get the optimal number of clustes
#' sol <- kselection(dat)
#' k   <- num_clusters(sol)
#' 
#' @seealso \code{\link{get_f_k}}
#' 
#' @rdname num_clusters
#' @export num_clusters
num_clusters <- function(obj) {
  UseMethod("num_clusters")
}

#' @method num_clusters default
#' @export
num_clusters.default <- function(obj) {
  NULL
}

#' @method num_clusters Kselection
#' @export
num_clusters.Kselection <- function(obj) {
  obj$k
}

#' @method plot Kselection
#' @export
plot.Kselection <- function(x, ...) {
  max_y <- 1.1 * max(x$f_k)
  plot(x$f_k,
       type = 'b',
       xlab = 'Number of clusters k',
       ylab = 'f(k)',
       ylim = c(0, max_y))
  lines(k, x$f_k[x$k],
        col  = 'red',
        pch  = 19,
        type = 'p') 
}

#' @method print Kselection
#' @export
print.Kselection <- function(x, ...) {
  if (num_clusters(x) == 1)
    cat('f(k) finds', num_clusters(x), 'cluster')
  else
    cat('f(k) finds', num_clusters(x), 'clusters')
}

# Calculate the weight factor for kselection
#
# Calculate the weight factor for kselection
#
# @param  n_d the number of dimensions (attributes)
# @param k the number of iterations
#
# @return an array with the weights
# 
# @author Daniel Rodriguez Perez
alpha_k <- function(n_d, k) {
  result <- 1 - 3/(4 * n_d)
  result <- rep(result, k)
  
  for (k in 3:k) {
    result[k] <- result[k - 1] + (1 - result[k - 1]) / 6
  }
  
  return(result)
}

# Calculate the optimal cluster for kselection
#
# Calculate the optimal cluster for kselection
#
# @param f_k the f(k) array
# @param k_threshold cluster selection threshold for f(k) value.
#
# @return the number of clusters
# 
# @author Daniel Rodriguez Perez
which_cluster <- function(f_k, k_threshold) {
  k <- which(f_k == min(f_k) & f_k < k_threshold)
  if (length(k) == 0)
    return(1)
  else
    return(min(k))
}