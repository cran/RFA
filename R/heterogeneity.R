#Méthode des L-moments
#On a d'abord besoin d'une fonction qui calcule les L-moments
samlmu <- function (x, nmom = 4, sort.data = T)
{
    xok <- x[!is.na(x)]
    n <- length(xok)
    if (nmom <= 0) return(numeric(0))
    if (nmom <= 2) rnames <- paste("l", 1:nmom, sep = "_")
    else rnames <- c("l_1", "t", paste("t", 3:nmom, sep = "_"))
    lmom <- rep(NA, nmom)
    names(lmom) <- rnames
    if (n == 0) return(lmom)
    if (sort.data == T) xok <- sort(xok)
    nmom.actual <- min(nmom, n)
    lmom[1] <- mean(xok)
    if (nmom.actual == 1) return(lmom)
    temp <- seq(1-n, n-1, by = 2)
    p1 <- rep(1, n)
    p <- temp/(n-1)
    lmom[2] <- mean(xok * p)
    if (nmom.actual == 2) return(lmom)
    if (xok[1] == xok[n]) {
        warning("all data values equal")
        return(lmom)
    }
    for (j in 3:nmom.actual) {
        p2 <- p1
        p1 <- p
        p <- ((2*j-3)*temp*p1 - (j-2)*(n+j-2)*p2) / ((j-1)*(n-j+1))
        lmom[j] <- mean(xok * p)/lmom[2]
    }

    lmom[2] <- lmom[2] / lmom[1]
    
    return(lmom)
}

#A function to calculate heterogeneity statistics of Hosking
heterogeneity <- function(N.sim,N.site,size,param,Vsite){
  #N.sim is the number of simulation of an homogenous region
  #N.site is the number of sites in this region
  #size is the length of each sites
  #param is a list of the parameter of the regional distribution
  #Vsite is the weigthed std of the at-site sample L-CV cf Hosking

  if (length(size)!=N.site){
    stop(paste('"size" must have a length equal to ',N.site,' !!!',sep=''))
  }
  
  V1 <- NULL
  V2 <- NULL
  V3 <- NULL
  loc <- param[1]
  scale <- param[2]
  shape1 <- param[3]
  shape2 <- param[4]

  for (i in 1:N.sim){
    tau <- NULL
    tau3 <- NULL
    tau4 <- NULL
    for (j in 1:N.site){
      data <- rkappa(size[j],loc,scale,shape1,shape2)
      temp <- samlmu(data)
      tau <- c(tau,temp[2])
      tau3 <- c(tau3,temp[3])
      tau4 <- c(tau4,temp[4])
    }
    
    tauR.sim <- sum(size*tau)/sum(size)
    tau3R.sim <- sum(size*tau3)/sum(size)
    tau4R.sim <- sum(size*tau4)/sum(size)

    V1 <- c(V1,sqrt(sum(size*(tau-tauR.sim)^2)/sum(size)))
    V2 <- c(V2,sum(size*sqrt( (tau-tauR.sim)^2 + (tau3-tau3R.sim)^2 ) ) / sum(size) )
    V3 <- c(V3,sum(size*sqrt( (tau3 - tau3R.sim)^2 + (tau4 - tau4R.sim)^2 ) ) / sum(size) )
    

  }

  m1 <- mean(V1)
  v1 <- sd(V1)
  m2 <- mean(V2)
  v2 <- sd(V2)
  m3 <- mean(V3)
  v3 <- sd(V3)

  H1 <- (Vsite[1]-m1)/v1
  H2 <- (Vsite[2]-m2)/v2
  H3 <- (Vsite[3]-m3)/v3

  cat('La valeur de la statistique H1 est :\n',H1,'\n')
  cat('La valeur de la statistique H2 est :\n',H2,'\n')
  cat('La valeur de la statistique H3 est :\n',H3,'\n')

  invisible(list(mean1=m1,std1=v1,tau=tau,V1=V1,H1=H1))
}

#A function who calculate the first four sample L-moments of the
#region and the weigthed std of the at-site sample L-CV cf Hosking
lmomreg <- function(sample.sites,index.flood=mean){
  ##sample.sites : a list giving data from each site
  ##index.flood   : a function who evaluate the index flood from the sample

  #sample.sites is a list containing the sample of each site
  n <- length(sample.sites)

  lmoments.sites <- NULL
  size <- NULL
  lmoments.regional <- rep(0,4)
  
  for (i in 1:n){
    temp <- samlmu(sample.sites[[i]] / index.flood(sample.sites[[i]]) )
    lmoments.sites <- rbind(lmoments.sites,temp)
    size <- c(size,length(sample.sites[[i]]))
    lmoments.regional <- lmoments.regional + size[i] * lmoments.sites[i,]
    
  }

  lmoments.regional <- lmoments.regional / sum(size)
  

  V1 <- 0
  V2 <- 0
  V3 <- 0
  for (i in 1:n){
    V1 <- V1 + size[i]*(lmoments.sites[i,2] - lmoments.regional[2])^2
    V2 <- V2 + size[i] * sqrt((lmoments.sites[i,2] - lmoments.regional[2])^2 + (lmoments.sites[i,3] - lmoments.regional[3])^2 )
    V3 <- V3 + size[i] * sqrt((lmoments.sites[i,3] - lmoments.regional[3])^2 + (lmoments.sites[i,4] - lmoments.regional[4])^2 )
  }

  V1 <- sqrt(V1/sum(size))
  V2 <- V2 / sum(size)
  V3 <- V3 / sum(size)
  names(V1) <- 'V1'
  names(V2) <- 'V2'
  names(V3) <- 'V3'

  return(list(site=lmoments.sites,reg=lmoments.regional,V=c(V1,V2,V3)))
}

regdist <- function(lmom, loc, main, xlab,
                    ylab, xliminf, xlimsup, draw.site = TRUE, ...){
  ##lmom : a list. Return of function 'lmomreg'
  ##loc  : Optional. Numeric specifying the value of the
  ##       location parameter
  
  lmom.reg <- lmom$reg
  lmom.site <- lmom$site
  n.site <- length(lmom.site[,1])
  
  if (missing(loc)){
    shape <- -(1 - 3*lmom.reg[3]) / (1 + lmom.reg[3])
    scale <- (1-shape) * (2-shape) * lmom.reg[2] * lmom.reg[1]
    loc <- lmom.reg[1] - (2-shape) * lmom.reg[2] * lmom.reg[1]
  }
  else{
    shape <- - (lmom.reg[1] - loc) / ( lmom.reg[2] * lmom.reg[1] ) + 2
    scale <- (1 - shape)*(lmom.reg[1] - loc)
  }
  
  param.reg <- c(loc, scale, shape)
  names(param.reg) <- c('loc', 'scale', 'shape')
  
  reg.fun <- function(p) qgpd(p,loc,scale,shape)
  
  eps <- 10^(-3)
  
  if (missing(main)) main <- 'Return Level Plot for the Regional distribution'
  if (missing(xlab)) xlab <- 'Return Period (Years)'
  if (missing(ylab)) ylab <- 'Return Level'
  if (missing(xlimsup)) xlimsup <- .999
  if (missing(xliminf)) xliminf <- .001
  
  plot(reg.fun, from = xliminf, to = xlimsup, log='x',
       xlab = xlab, ylab = ylab, main = main, ...)
  
  if (draw.site){
    for (i in 1:n.site){
      if (missing(loc)){
        shape <- -(1 - 3*lmom.site[i,3]) / (1 + lmom.site[i,3])
        scale <- (1-shape) * (2-shape) * lmom.site[i,2] * lmom.site[i,1]
        loc <- lmom.site[i,1] - (2-shape) * lmom.site[i,2] * lmom.site[i,1]
      }
      else{
        shape <- - (lmom.site[i,1] - loc) / ( lmom.site[i,2] *
                                             lmom.site[i,1] ) + 2
        scale <- (1 - shape)*(lmom.site[i,1] - loc)
      }
      
      site.fun <- function(p) qgpd(p,loc,scale,shape)

      plot(site.fun, from= eps, to = xlimsup, add = TRUE,
           lty = 2, col = 'grey')
    }
  }

  return(param.reg)
}

locdist <- function(param.reg, mu, data, main,
                    xlab, ylab, xlimsup, index.flood = mean,
                    draw.data = TRUE, ...){
  ##param.reg   : vector of the regional distribution parameters
  ##data        : the at-site sample
  ##index.flood : A function who computes the index flood
  
  loc <- param.reg[1]
  scale <- param.reg[2]
  shape <- param.reg[3]
  
  atsite.fun <- function(T){
    p <- 1 - 1 / (mu*T)
    return(index.flood(data)*qgpd(p,loc,scale,shape))
  }
  
  eps <- 10^(-3)
  
  if (missing(main)) main <- 'Return Level Plot for the at-site distribution'
  if (missing(xlab)) xlab <- 'Return Period (Years)'
  if (missing(ylab)) ylab <- 'Return Level'
  if (missing(xlimsup)) xlimsup <- 50
  
  plot(atsite.fun, from= 1 / mu + eps, to = xlimsup, log='x',
       xlab = xlab, ylab = ylab, main = main, ...)
  
  if (draw.data){
    p.emp <- ppoints(data)
    T.emp <- 1 / (mu * (1 - p.emp))
  }
  
  points(sort(T.emp), sort(data) )
  
  param <- c( c(loc, scale) * index.flood(data), shape)
  names(param) <- c('loc','scale','shape')
  
  return(param)
}
