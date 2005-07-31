#A function whi generate a random sample of a four parameter kappa distribution
rkappa <- function(n,loc,scale,shape1,shape2){
  p <- runif(n)
  
  if (scale <=0){
    stop('Invalid scale !')
  }
  if (shape2 == -1){
    #This is the special case of the generalised logistic distribution
    return(loc  -  scale / shape1 * (1 -(1/runif(n) - 1)^(-shape1)) )
  }
  if (shape2 == 0 ){
    #This is the special case of the GEV
    return(rgev(n,loc,scale,shape1))
  }
  if (shape2 == 1){
    #This is the special case of the GPD
    return(rgpd(n,loc,scale,shape1))
  }

  #And so, this is the four parameter Kappa distribution
  return(loc - scale / shape1 * ( 1 - ( (1 -p ^shape2) / shape2 ) ^ (-shape1) ))
}

rgpd <- function(n, loc = 0, scale = 1, shape = 0){
    if (min(scale) < 0) 
        stop("invalid scale")
    if (length(shape) != 1) 
        stop("invalid shape")
    if (shape == 0) 
        return(loc + scale * rexp(n))
    else return(loc + scale * (runif(n)^(-shape) - 1)/shape)
}

rgev <- function(n, loc = 0, scale = 1, shape = 0){
    if (min(scale) < 0) 
        stop("invalid scale")
    if (length(shape) != 1) 
        stop("invalid shape")
    if (shape == 0) 
        return(loc - scale * log(rexp(n)))
    else return(loc + scale * (rexp(n)^(-shape) - 1)/shape)
}

#A function who fit the four parameter kappa distribution with L-moments
kappalmom <- function(lmom){

  param <- .Fortran("pelkap",lmom,rep(0,5), PACKAGE = "RFA")[[2]]
  if (param[5] == 0){
    param <- param[-5]
    names(param) <- c('loc','scale','shape1','shape2')
    return(param)
  }
  else
    cat(paste('Error in optimizing !!!\nThe fail flag was ',
              param[5],'.\n',sep=''))
}

## A function to compute quantile for the GPD
qgpd <- function (p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE){
  
  if (min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >= 
      1) 
    stop("`p' must contain probabilities in (0,1)")
  if (min(scale) < 0) 
    stop("invalid scale")
  if (length(shape) != 1) 
    stop("invalid shape")
  if (lower.tail) 
    p <- 1 - p
  if (shape == 0) 
    return(loc - scale * log(p))
  else return(loc + scale * (p^(-shape) - 1)/shape)
  
}
