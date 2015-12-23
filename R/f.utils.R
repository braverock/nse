
# Function to check if x is of good dimension
f.error.multivariate = function(x){
  if (!is.vector(x) && dim(x)[2] != 1) {
    stop("this function only handle univariate time-series")
  }
}