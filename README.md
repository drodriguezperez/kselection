# kselection: Selection of k in k-means clustering

`kselection`implements Pham, Dimov and Nguyen from 2004. 

## Installation

To install the latest development builds of `kselection` directly from GitHub, run this instead:

```
require(devtools)
install_github('kselection', 'drodriguezperez')
```

## Features

`kselection` implements the method proposed by Pham, Dimov and Nguyen for selecting the number of clusters for the K-means algorithm. In this method a function $f(K)$ is used to evaluate the quality of the resulting clustering and help decide on the optimal value of $K$ for each data set.

## Use
```
# Create a data set with two clusters
dat <- matrix(c(rnorm(100, 2, .1), rnorm(100, 3, .1),
                rnorm(100, -2, .1), rnorm(100, -3, .1)), 200, 2)
# Ejecute the method
sol <- kselection(dat)
 
# Get the results
k   <- num_clusters(sol) # optimal number of clustes
f_k <- get_f_k(sol)      # the f(k) vector

# Plot the results
plot(sol)
```

## References

1. D T Pham, S S Dimov, and C D Nguyen "Selection of k in k-means clustering", Mechanical Engineering Science, 2004, pp. 103-119.