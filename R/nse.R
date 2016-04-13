#'  Variance of sample mean of functional of reversible Markov chain using methods of Geyer (1992).
#' @description Calculate Geyer (1992) NSE estimator.
#' @details  The type "iseq" is a wrapper around \link[mcmc]{initseq} from the 
#' MCMC package and gives the positive intial sequence estimator.
#'  The type "bm" is the batch mean estimator.
#'  The type "obm" is the overlapping batch mean estimator.
#'  The type "iseq.bm" is a combinaison of the two.
#' @param x A numeric vector or a matrix(only for type "bm").
#' @param type The type c("iseq","bm","obm","iseq.bm").
#' @param nbatch An optional parameter for the type m and iseq.bm.
#' @param iseq.type constraints on function, ("pos") nonnegative, ("dec") nonnegative 
#' and nonincreasing, and ("con") nonnegative, nonincreasing, and convex. The default value is "pos".
#' @import mcmc
#' @references Geyer, Charles J. "Practical markov chain monte carlo." Statistical Science (1992): 473-483.
#' @return  The variance estimator in the univariate case or the v
#' ariance-covariance matrix estimator in the multivariate case.
#' @examples
#' n    = 1000
#' ar   = 0.9
#' mean = 1
#' sd   = 1
# 
#' set.seed(1234)
#' x = as.vector(arima.sim(n = n, list(ar = ar), sd = sd) + mean)
# 
#' nse.geyer(x = x, type = "bm", nbatch = 30)
#' nse.geyer(x = x, type = "obm", nbatch = 30)
#' nse.geyer(x = x, type = "iseq", iseq.type = "pos")
#' nse.geyer(x = x, type = "iseq.bm", iseq.type = "con")
#'  
#'@export
nse.geyer = function(x, type = c("iseq", "bm", "obm", "iseq.bm"), nbatch = 30, iseq.type = c("pos", "dec", "con")) {
  
  if (is.vector(x)) {
    x = matrix(x,ncol = 1)
  }
  n = dim(x)[1]
  
  type = type[1]  
  if (type == "iseq") {
    
    f.error.multivariate(x) 
    
    if (iseq.type == "pos") {
      var.iseq = mcmc::initseq(x = x)$var.pos
    } else if(iseq.type == "dec") {
      var.iseq = mcmc::initseq(x = x)$var.dec
    } else if(iseq.type == "con") {
      var.iseq = mcmc::initseq(x = x)$var.con
    } else {
      stop("Invalid iseq.type : must be one of c('pos','dec','con')")
    }
    
    iseq = var.iseq / n # Intial sqequence Geyer (1992)
    out = iseq
    
  } else if (type == "bm") {
    
    ncol  = dim(x)[2]
    x     = as.data.frame(x)
    batch = matrix(unlist(lapply(split(x, ceiling(seq_along(x[,1]) / (n / nbatch))), FUN = function(x) colMeans(x))), ncol = ncol,byrow = TRUE)
    out   = var(x = batch) / (nbatch - 1)
    
    if (is.matrix(out) && dim(out) == c(1,1)) {
      out = as.vector(out)
    }
    
  } else if(type == "obm") {
    
    out = as.vector(mcmcse::mcse(x, method = "obm")$se^2)
    
  } else if(type == "iseq.bm"){
    
    f.error.multivariate(x)
    batch = unlist(lapply(split(x, ceiling(seq_along(x) / (n / nbatch))), FUN = mean))
    
    iseq.type = iseq.type[1]
    if(iseq.type == "pos") {
      var.iseq = mcmc::initseq(batch)$var.pos
    } else if(iseq.type == "dec") {
      var.iseq = mcmc::initseq(x = batch)$var.dec
    } else if(iseq.type == "con") {
      var.iseq = mcmc::initseq(x = batch)$var.con
    } else {
      stop("Invalid iseq.type : must be one of c('pos','dec','con')")
    }
    
    iseq.bm = var.iseq / nbatch
    out = iseq.bm
    
  } else {
    stop("Invalid type : must be of type c('iseq','bm','iseq.bm')")
  }
  out = unname(out)
  return(out)
}


#' The spectral density at zero.
#' @description Calculate the variance of the mean with the spectrum at zero estimator.
#' @details  This is a wrapper around \link[coda]{spectrum0.ar} from the CODA package, \link[sapa]{SDF} from the 
#' sapa package, \link[psd]{psdcore} from the PSD package and \link[mcmcse]{mcse}from the mcmcse package.
#' @param x A numeric vector.
#' @param type A character string denoting the method to use in estimating the spectral density function
#' @param lag.prewhite Prewhite the serie before analysis (integer or NULL, i.e. automatic selection)
#' @return The variance estimator.
#' @references Plummer, Martyn, et al. "CODA: Convergence diagnosis and output analysis 
#' for MCMC." R news 6.1 (2006): 7-11.
#' @references D.B. Percival and A. Walden "Spectral Analysis for Physical Applications: 
#' Multitaper and Conventional Univariate Techniques". Cambridge University Press (1993).
#' @references Barbour, A. J. and R. L. Parker (2014), "psd: Adaptive, sine multitaper power spectral 
#' density estimation for R", Computers & Geosciences Volume 63 February 2014 : 1-8
#' @references James M. Flegal, John Hughes and Dootika Vats. (17-08-2015). mcmcse: Monte Carlo 
#' Standard Errors for MCMC. R package version 1.1-2. Riverside, CA and Minneapolis, MN.
#' @import coda sapa psd mcmcse
#' @examples 
#' n    = 1000
#' ar   = 0.9
#' mean = 1
#' sd   = 1
#' set.seed(1234)   
#' x = as.vector(arima.sim(n = n, list(ar = ar), sd = sd) + mean)
#'  
#' nse.spec0(x = x, type = "ar", lag.prewhite = 0)
#' nse.spec0(x = x, type = "ar", lag.prewhite = 1)
#' nse.spec0(x = x, type = "ar", lag.prewhite = NULL)
#' 
#' nse.spec0(x = x, type = "glm", lag.prewhite = 0)
#' nse.spec0(x = x, type = "glm", lag.prewhite = 1)
#' nse.spec0(x = x, type = "glm", lag.prewhite = NULL)
#' 
#' nse.spec0(x = x, type = "wosa", lag.prewhite = 0)
#' nse.spec0(x = x, type = "wosa", lag.prewhite = 1)
#' nse.spec0(x = x, type = "wosa", lag.prewhite = NULL)
#' 
#' nse.spec0(x = x, type = "bartlett", lag.prewhite = 0)
#' nse.spec0(x = x, type = "bartlett", lag.prewhite = 1)
#' nse.spec0(x = x, type = "bartlett", lag.prewhite = NULL)
#' 
#' nse.spec0(x = x, type = "tukey", lag.prewhite = 0)
#' nse.spec0(x = x, type = "tukey", lag.prewhite = 1)
#' nse.spec0(x = x, type = "tukey", lag.prewhite = NULL)
#'  
#'@export
nse.spec0 = function(x, type = c("ar", "glm", "wosa", "bartlett", "tukey"), lag.prewhite = 0) {
  
  scale = 1
  if (is.vector(x)){
    x = matrix(x,ncol = 1)
  }
  f.error.multivariate(x)
  
  n = dim(x)[1]
  # !!! DA old implementation
  #browser()
  #ar.model = psd::prewhiten(as.vector(x), AR.max = 1, plot = FALSE)
  #x = ar.model$prew_ar
  #scale = 1 / (1 - ar.model$ardfit$ar)^2
  #if (length(scale) == 0) {
  #  scale = 1
  #}
  tmp   = f.prewhite(x, ar.order = lag.prewhite) 
  x     = tmp$ar.resid
  scale = tmp$scale
  
  type = type[1]
  if (type == "ar") {
    spec0 = coda::spectrum0.ar(x)$spec
  }else if (type == "glm") {
    spec0 = coda::spectrum0(x)$spec
  }else if (type == "wosa") {
    spec0 = sapa::SDF(x, method = "wosa", single.sided = TRUE)[1]
  }else if (type == "tukey") {
    out = as.vector((mcmcse::mcse(x, method = "tukey")$se)^2) * scale
    return(out)
  }else if (type == "bartlett") {
    out = as.vector((mcmcse::mcse(x, method = "bartlett")$se)^2) * scale
    return(out)
  }else {
    stop("Invalid type : must be one of c('ar','bartlett','wosa','tukey')")
  }
  spec0 = spec0 * scale
  out   = spec0 / n
  out   = unname(out)
  return(out)
}
#' Newey-West NSE estimators.
#' @description Calculate the variance of the mean with the Newey West (1987, 1994) HAC estimator.
#' @description This is a wrapper around \link[sandwich]{lrvar} from the sandwich package.
#' @param x A numeric vector or matrix.
#' @param lag.prewhite Prewhite the serie before analysis (integer or NULL, i.e. automatic selection)
#' @return The variance estimator in the univariate case or the variance-covariance matrix 
#' estimator in the multivariate case.
#' @references Andrews, Donald WK. "Heteroskedasticity and autocorrelation consistent covariance 
#' matrix estimation." Econometrica: Journal of the Econometric Society 59.03 (1991): 817-858.
#' @references Newey, Whitney K., and Kenneth D. West. "A simple, positive semi-definite, 
#' heteroskedasticity and autocorrelationconsistent covariance matrix.", Econometrica: Journal of 
#' the Econometric Society 55.03 (1987) : 703-708.
#' @references Newey, Whitney K., and Kenneth D. West. "Automatic lag selection in covariance 
#' matrix estimation." The Review of Economic Studies 61.4 (1994): 631-653.
#' @references Zeileis, Achim. "Econometric computing with HC and HAC 
#' covariance matrix estimators." (2004).
#' @import sandwich
#' @examples 
#' n    = 1000
#' ar   = 0.9
#' mean = 1
#' sd   = 1
#'
#' set.seed(1234)   
#' x = as.vector(arima.sim(n = n, list(ar = ar), sd = sd) + mean)
#'  
#' nse.nw(x = x, lag.prewhite = 0)
#' nse.nw(x = x, lag.prewhite = 1)
#' nse.nw(x = x, lag.prewhite = NULL)
#'@export
nse.nw <- function(x, lag.prewhite = 0) {
  
  tmp = f.prewhite(x, ar.order = lag.prewhite) 
  lag.prewhite = tmp$ar.order
  out = sandwich::lrvar(x = x, type = "Newey-West", prewhite = lag.prewhite, adjust = TRUE)
  out = unname(out)
  return(out)
}

#' Andrews NSE estimators.
#' @description Calculate the variance of the mean with the kernel based variance estimator indtroduced by Andrews (1991).
#' @details  This is a wrapper around \link[sandwich]{lrvar} from the sandwich package and use Andrews (1991) automatic bandwidth estimator.
#' @param x A numeric vector or matrix.
#' @param lag.prewhite Prewhite the serie before analysis (integer or NULL, i.e. automatic selection)
#' @param type The type of kernel used c("Bartlett","Parzen","Quadratic Spectral","Truncated","Tukey-Hanning").
#' @return The variance estimator in the univariate case or the variance-covariance matrix estimator in the multivariate case.
#' @references Zeileis, Achim. "Econometric computing with HC and HAC covariance matrix estimators." (2004).
#' @references Andrews, Donald WK. "Heteroskedasticity and autocorrelation consistent 
#' covariance matrix estimation." Econometrica: Journal of the Econometric Society 59.03 (1991): 817-858.
#' @references Newey, Whitney K., and Kenneth D. West. "A simple, positive semi-definite, 
#' heteroskedasticity and autocorrelationconsistent covariance matrix.", Econometrica: Journal of the Econometric Society 55.03 (1987) : 703-708.
#' @import sandwich
#' @examples 
#' n    = 1000
#' ar   = 0.9
#' mean = 1
#' sd   = 1
#'  
#' set.seed(1234)   
#' x = as.vector(arima.sim(n = n, list(ar = ar), sd = sd) + mean)
#'  
#'nse.andrews(x = x, type = "bartlett", lag.prewhite = 0)
#'nse.andrews(x = x, type = "bartlett", lag.prewhite = 1)
#'nse.andrews(x = x, type = "bartlett", lag.prewhite = NULL)
#'
#'nse.andrews(x = x, type = "parzen", lag.prewhite = 0)
#'nse.andrews(x = x, type = "parzen", lag.prewhite = 1)
#'nse.andrews(x = x, type = "parzen", lag.prewhite = NULL)
#'
#'nse.andrews(x = x, type = "tukey", lag.prewhite = 0)
#'nse.andrews(x = x, type = "tukey", lag.prewhite = 1)
#'nse.andrews(x = x, type = "tukey", lag.prewhite = NULL)
#'
#'nse.andrews(x = x, type = "qs", lag.prewhite = 0)
#'nse.andrews(x = x, type = "qs", lag.prewhite = 1)
#'nse.andrews(x = x, type = "qs", lag.prewhite = NULL)
#'
#'nse.andrews(x = x, type = "trunc", lag.prewhite = 0)
#'nse.andrews(x = x, type = "trunc", lag.prewhite = 1)
#'nse.andrews(x = x, type = "trunc", lag.prewhite = NULL)
#'@export
nse.andrews <- function(x, type = c("bartlett", "parzen", "tukey", "qs", "trunc"), lag.prewhite = 0) {
  
  tmp = f.prewhite(x, ar.order = lag.prewhite) 
  lag.prewhite = tmp$ar.order
  type.sandwich = f.type.sandwich(type.in = type)
  out = sandwich::lrvar(x = x, type = "Andrews", prewhite = lag.prewhite, adjust = TRUE, kernel = type.sandwich)
  out = unname(out)
  return(out)
}

#' Hirukawa NSE estimators.  
#' @description Calculate the variance of the mean with the kernel based variance estimator 
#' by Andrews (1991) using Hirukawa (2010) automatic bandwidth estimator.
#' @details This is a wrapper around \link[sandwich]{lrvar} from the sandwich package 
#' and use Hirukawa (2010) automatic bandwidth estimator.
#' @param x A numeric vector.
#' @param lag.prewhite Prewhite the serie before analysis (integer or NULL, i.e. automatic selection)
#' @param type The type of kernel used c("Bartlett","Parzen").
#' @references Zeileis, Achim. "Econometric computing with HC and HAC covariance matrix estimators." (2004).
#' @references Andrews, Donald WK. "Heteroskedasticity and autocorrelation consistent covariance 
#' matrix estimation." Econometrica: Journal of the Econometric Society 59.03 (1991): 817-858.
#' @references Hirukawa, Masayuki. "A two-stage plug-in bandwidth selection and its implementation 
#' for covariance estimation." Econometric Theory 26.03 (2010): 710-743.
#' @import sandwich
#' @return The variance estimator.
#' @examples
#' n    = 1000
#' ar   = 0.9
#' mean = 1
#' sd   = 1
#' 
#' set.seed(1234) 
#' x = as.vector(arima.sim(n = n, list(ar = ar), sd = sd) + mean)
#'  
#' nse.hiruk(x = x, type = "bartlett", lag.prewhite = 0)
#' nse.hiruk(x = x, type = "bartlett", lag.prewhite = 1)
#' nse.hiruk(x = x, type = "bartlett", lag.prewhite = NULL)
#' 
#' nse.hiruk(x = x, type = "parzen", lag.prewhite = 0)
#' nse.hiruk(x = x, type = "parzen", lag.prewhite = 1)
#' nse.hiruk(x = x, type = "parzen", lag.prewhite = NULL)
#' 
#'@export
nse.hiruk <- function(x, type = c("bartlett", "parzen"), lag.prewhite = 0) {
  f.error.multivariate(x)
  tmp = f.prewhite(x, ar.order = lag.prewhite) 
  lag.prewhite = tmp$ar.order
  bandwidth = f.hiruk.bandwidth.solve(x, type = type, lag.prewhite = lag.prewhite)
  type.sandwich = f.type.sandwich(type.in = type)
  out = sandwich::lrvar(x = x, type = "Andrews", prewhite = lag.prewhite, adjust = TRUE, kernel = type.sandwich, bw = bandwidth)
  out = unname(out)
  return(out)
}


#' Bootstrap NSE estimators. 
#' @description Calculate the variance of the mean with a bootstrap variance estimator.
#' @details  Use the automatic blocksize in \link[np]{b.star} from th np package which is based 
#' on Politis and White (2004) and Patton and al (2009). 
#' Two bootstrap schemes are available; The stationary bootstrap of Politis and Romano  (1994)
#' and the circular bootstrap of Politis and Romano (1992).
#' @param x  A numeric vector or a matrix.
#' @param nb The number of bootstrap replication.
#' @param type The bootstrap schemes c("stationary","circular").
#' @param lag.prewhite Prewhite the serie before analysis (integer or NULL, i.e. automatic selection)
#' @return The variance estimator in the univariate case or the variance-covariance matrix 
#' estimator in the multivariate case.
#' @references Politis, Dimitris N., and Joseph P. Romano. "A circular block-resampling procedure 
#' for stationary data." Exploring the limits of bootstrap (1992): 263-270.
#' @references Politis, Dimitris N., and Halbert White. "Automatic block-length selection for the 
#' dependent bootstrap." Econometric Reviews 23.1 (2004): 53-70.
#' @references Patton, Andrew, Dimitris N. Politis, and Halbert White. "Correction to "Automatic 
#' block-length selection for the dependent bootstrap" by D. Politis and H. White." Econometric Reviews 28.4 (2009): 372-375.
#' @references Politis, Dimitris N., and Joseph P. Romano. "The stationary bootstrap." Journal of 
#' the American Statistical association 89.428 (1994): 1303-1313.
#' @references Hayfield, Tristen, and Jeffrey S. Racine. "Nonparametric econometrics: The np 
#' package." Journal of statistical software 27.5 (2008): 1-32.
#' @import np Rcpp
#' @useDynLib nse
#' @importFrom Rcpp evalCpp   
#' @examples  
#' n    = 1000
#' ar   = 0.9
#' mean = 1
#' sd   = 1
#' 
#' set.seed(1234) 
#' x = as.vector(arima.sim(n = n, list(ar = ar), sd = sd) + mean)
#'  
#' set.seed(1234)
#' nse.boot(x = x, nb = 1000, type = "stationary", b = NULL, lag.prewhite = 0)
#' nse.boot(x = x, nb = 1000, type = "stationary", b = NULL, lag.prewhite = 1)
#' nse.boot(x = x, nb = 1000, type = "stationary", b = NULL, lag.prewhite = NULL)
#' 
#' nse.boot(x = x, nb = 1000, type = "stationary", b = 10, lag.prewhite = 0)
#' nse.boot(x = x, nb = 1000, type = "stationary", b = 10, lag.prewhite = 1)
#' nse.boot(x = x, nb = 1000, type = "stationary", b = 10, lag.prewhite = NULL)
#' 
#' nse.boot(x = x, nb = 1000, type = "circular", b = NULL, lag.prewhite = 0)
#' nse.boot(x = x, nb = 1000, type = "circular", b = NULL, lag.prewhite = 1)
#' nse.boot(x = x, nb = 1000, type = "circular", b = NULL, lag.prewhite = NULL)
#' 
#' nse.boot(x = x, nb = 1000, type = "circular", b = 10, lag.prewhite = 0)
#' nse.boot(x = x, nb = 1000, type = "circular", b = 10, lag.prewhite = 1)
#' nse.boot(x = x, nb = 1000, type = "circular", b = 10, lag.prewhite = NULL)
#'@export
nse.boot <- function(x, nb, type = c("stationary", "circular"), b = NULL, lag.prewhite = 0){
  
  x = as.vector(x)
  
  # prewhiteneing
  tmp   = f.prewhite(x, ar.order = lag.prewhite) 
  x     = tmp$ar.resid
  scale = tmp$scale
  
  # optimal block length selection
  if (is.null(b)){
    blockSize = np::b.star(data = x, round = TRUE)
    if(type == "stationary"){
      b = as.numeric(blockSize[1,1])
    } else if(type == "circular"){
      b = as.numeric(blockSize[1,2])
    }
  }
  
  out = scale * var(f.bootstrap(x = x, nb = nb, statistic = colMeans, b = b, type = type)$statistic)
  out = as.vector(out)
  out = unname(out)
  return(out)
}
