
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
  
  v.x = var(x)
  v.e = var(ar.resid)
  
  scale = 1
  if (ar.order == 1) {
    if (abs(ar.param) < 1) {
      scale = 1 / (1 - ar.param^2)
    } else {
      scale = v.x / v.e
    }
  }
  if (ar.order > 1) {
    # DA here we should check stationarity
    scale = v.x / v.e
  }
  
  out = list(ar.resid = ar.resid, ar.param = ar.param, ar.order = ar.order, scale = scale)
  return(out)
  #n = n - prewhite
}