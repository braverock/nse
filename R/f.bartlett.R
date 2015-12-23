# Estimation of the generalized derivative of spectral density
# at the origin with the Bartlett kernel
#
f.bartlett = function(g, st, q) {

  size = length(g)
  if (q == 0) {
    sbt = (t(g) %*% g) / size 
  } else {
    sbt = 0
  }

  if (st <= size - 1) {
    maxlagbt = st
  } else {
    maxlagbt = size - 1
  }

  jbt = 1
  while (jbt < maxlagbt) {
    xbt  = jbt / maxlagbt
    s0bt = (t(g[1:(size - jbt)]) %*% g[(1 + jbt):size]) / size 

    if (xbt <= 1) {
      wbt = 1 - xbt
    } else {
      wbt = 0
    }

    sbt = sbt + wbt * (jbt ^ q) * (s0bt + s0bt)
    jbt = jbt + 1
  }

  return(sbt)
}
