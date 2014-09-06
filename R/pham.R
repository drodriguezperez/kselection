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
#' @param k_threshold cluster selection threshold for f(k) value.
#' @param iter.max the maximum number of iterations allowed.
#' @param nstart if centers is a number, how many random sets should be chosen?
#' @param trace show a progress bar
#' 
#' @return an array with f(k) metric
#' 
#' @author Daniel Rodriguez Perez
#' 
#' @references
#' D T Pham, S S Dimov, and C D Nguyen "Selection of k in k-means clustering",
#' Mechanical Engineering Science, 2004, pp. 103-119.
#' 
#' @import tools
#' 
#' @rdname kselection
#' @export kselection
kselection <- function(x,
                       max_centers = 15,
                       k_threshold = 0.85,
                       iter.max    = 200,
                       nstart      = 100,
                       trace       = FALSE) {
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
    mod_info <- kmeans(x,
                       centers  = k,
                       iter.max = iter.max,
                       nstart   = nstart)
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
  
  result <- list(k    = which_cluster(f_k, k_threshold),
                 f_k  = f_k,
                 s_k  = s_k)
  class(result) <- 'Kselection'
  
  return(result)
}

#' Get the number of clusters
#' 
#' Get the number of clusters
#' 
#' @param obj the output of kselection function
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
       xlab = 'Number of clusters k',
       ylab = 'f(k)',
       ylim = c(0, max_y))
}

#' @method plot Kselection
#' @export
print.Kselection <- function(obj) {
  cat('The f(k) function finds', obj$k, 'clusters')
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

which_cluster <- function(f_k, k_threshold) {
  k <- which(f_k == min(f_k) & f_k < k_threshold)
  if (length(k) == 0)
    return(1)
  else
    return(min(k))
}