h.boot.sample <-
  function(x, b, type)
  {
    if (type == "stationary") {type = 0}
    else if (type == "circular") {type = 1}
    else {stop("only stationary and circular bootstrap implemented")}
    return(f_bootstrap(x, b, type))
  }

.f.bootstrap <- function(x, nb = 1, statistic = NULL, b = NULL, type, ...)
{
  y <- embed(x, 1)
  if(is.null(statistic)) {
    n = NROW(y)
    boot <- matrix(y, nrow=n, ncol=nb)
    out <- apply(boot, 2, h.boot.sample, b, type)
    return(drop(out))
  }
  else {
    
    yi <- 1:NROW(y)
    orig.statistic <- statistic(y,...)
    l.stat <- length(orig.statistic)
    stat <- matrix(0, nb, l.stat)
    for(i in 1:nb){
      stat[i,] <- statistic(as.matrix(y[h.boot.sample(yi, b, type),]),...)
    }
    out <- list(statistic = stat)
    
    return(out)
  }
}

# Bootstrap samples generation
# @description Generate bootstrap samples
# @details Two bootstrap schemes are available, the stationary bootstrap of Politis and Romano  (1994)
# and the circular bootstrap of Politis and Romano (1992)
#     @param x       A numeric vector or a matrix
#     @param nb   The number of bootstrap replication
#     @param statistic Function to be used on the bootstrap replication
#     @param b  The block length
#     @param type    The bootstrap schemes c("stationary","circular")
#     @param seed    The seed for the bootstrap simulation
#     @param ...    Aditional argument passed to the statistic function
#     @return  The bootstrap samples or the statistic if statistic is non-null
#    @references Politis, Dimitris N., and Joseph P. Romano. "A circular block-resampling procedure for stationary data." Exploring the limits of bootstrap (1992): 263-270.
#    @references Politis, Dimitris N., and Halbert White. "Automatic block-length selection for the dependent bootstrap." Econometric Reviews 23.1 (2004): 53-70.
#    @references Politis, Dimitris N., and Joseph P. Romano. "The stationary bootstrap." Journal of the American Statistical association 89.428 (1994): 1303-1313.
#@import np
f.bootstrap = compiler::cmpfun(.f.bootstrap)
