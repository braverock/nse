
#' Function to perform AR prewithening
f.prewhite = function(x, ar.order = 1) {
  
  if (!is.null(ar.order)) {
    if (ar.order == 0){
      out = list(ar.resid = x, ar.param = NULL, ar.order = 0, scale = 1)
      return(out)
    } 
  }
  
  aic = FALSE
  if (is.null(ar.order)){
    aic = TRUE
    ar.order = min(length(x), 10)
  }
  ar.fit = stats::ar.ols(x = x, aic = aic, order.max = ar.order, demean = FALSE)
  
  if (inherits(ar.fit, "try-error")) {
    stop("AR prewhitening of estimating functions failed")
  } 
    
  ar.resid = as.vector(na.omit(ar.fit$resid))
  ar.param = as.vector(ar.fit$ar)
  ar.order = as.numeric(ar.fit$order)
  
  scale = 1 / (1 - sum(ar.param^2))
  if (length(scale) == 0) {
    scale = 1
  }
  
  out = list(ar.resid = ar.resid, ar.param = ar.param, ar.order = ar.order, scale = scale)
  return(out)
  #n = n - prewhite
}