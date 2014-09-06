##
##  kselection tests
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

context("Tests for kselection")

test_that("evaluate alpha_k calculations", {
  n_d    <- 3
  k      <- 4
  a_k    <- 1 - 3/(4 * n_d)
  a_k    <- rep(a_k, k)
  a_k[3] <- a_k[2] + (1- a_k[2]) / 6
  a_k[4] <- a_k[3] + (1- a_k[3]) / 6
  
  expect_that(alpha_k(n_d, k), equals(a_k))
  
  a_k[5] <- a_k[4] + (1- a_k[4]) / 6
  a_k[6] <- a_k[5] + (1- a_k[5]) / 6
  a_k[7] <- a_k[6] + (1- a_k[6]) / 6
  a_k[8] <- a_k[7] + (1- a_k[7]) / 6
  
  expect_that(alpha_k(n_d, 8), equals(a_k))
})

test_that("evaluate which_cluster calculations", {
  f_k    <- rep(1, 10)
  expect_that(which_cluster(f_k, k_threshold = 0.85), equals(1))
  
  f_k[3] <- 0.5
  expect_that(which_cluster(f_k, k_threshold = 0.85), equals(3))
  expect_that(which_cluster(f_k, k_threshold = 0.45), equals(1))
  
  f_k[5] <- 0.5
  expect_that(which_cluster(f_k, k_threshold = 0.85), equals(3))
  expect_that(which_cluster(f_k, k_threshold = 0.45), equals(1))
})

test_that("evaluate invalid inputs values", {
  expect_that(kselection(NULL),
              throws_error("'data' must be of a vector type, was 'NULL'"))
  expect_that(kselection('test_data'), 
              throws_error('x must contain numerical data'))
  
  test_data    <- as.matrix(1:300, 100, 3)
  test_data    <- as.data.frame(test_data)
  expect_that(kselection(test_data, max_centers = 1),
              throws_error("'max_centers' must be greater than 2"))
  
  test_data$V3 <- 'a'
  expect_that(kselection(test_data), 
              throws_error('x must contain numerical data'))
})

test_that("evaluate the solution", {
  x <- matrix(c(rnorm(100, 2, .1), rnorm(100, 3, .1),
                rnorm(100, -2, .1), rnorm(100, -3, .1)), 200, 2)
  k <- kselection(x)
  
  expect_that(num_clusters(x), is_null())
  
  expect_that(class(k), equals('Kselection'))
  expect_that(k$k, equals(2))
  expect_that(num_clusters(k), equals(2))
})
