lmomplot <- function(lmom.reg, which = 1:3,
                      ask = nb.fig < length(which)
                      && dev.interactive(), pch = 15, draw.dist = TRUE){
  ##lmom.reg : a list composed with the first 4 L-moments
  ##           of each site, the first 4 regional L-moments
  ##which    : a vector specifying which L-moment plot should be
  ##           plot
  ##ask      : logical. If TRUE, each new plot are asking for
  ##           confirmation before being drawn.

  if (!is.numeric(which) || any(which < 1) || any(which > 4)) 
        stop("`which' must be in 1:4")
  show <- rep(FALSE, 3)
  show[which] <- TRUE
  nb.fig <- prod(par("mfcol"))
  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  }
  if (show[1]) {
    plot(lmom.reg$site[,3], lmom.reg$site[,2],
         xlab = expression(paste('L-Skewness,  ',tau[3],sep='')),
         ylab = expression(paste('L-CV,  ',tau,sep='')))
    points(lmom.reg$reg[3], lmom.reg$reg[2], pch = pch)
  }
  if (show[2]) {
    plot(lmom.reg$site[,4], lmom.reg$site[,2],
         xlab = expression(paste('L-Kurtosis,  ',tau[4],sep='')),
         ylab = expression(paste('L-CV,  ',tau,sep='')))
    points(lmom.reg$reg[4], lmom.reg$reg[2], pch = pch)
  }
  if (show[3]) {
    plot(lmom.reg$site[,3], lmom.reg$site[,4],
         ylab = expression(paste('L-Kurtosis,  ',tau[4],sep='')),
         xlab = expression(paste('L-Skewness,  ',tau[3],sep='')),
         xlim = c(0,0.5), ylim = c(0,0.4) )
    points(lmom.reg$reg[3], lmom.reg$reg[4], pch = pch)
    
    if (draw.dist){
      GP <- function(tau3) 0.20196*tau3 + 0.95924*tau3^2 -
        0.20096*tau3^3 + 0.04061*tau3^4
      GEV <- function(tau3) 0.10701 + 0.11090*tau3 +
        0.84838*tau3^2 - 0.06669*tau3^3 + 0.00567*tau3^4 -
          0.04208*tau3^5 + 0.03763*tau3^6
      GL <- function(tau3) 0.16667 + 0.83333*tau3^2
      LN3 <- function(tau3) 0.12282 + 0.77518*tau3^2 +
        0.12279*tau3^4 - 0.13638*tau3^6 + 0.11368*tau3^8
      PE3 <- function(tau3) 0.12240 + 0.30115*tau3^2 +
        0.95812*tau3^4 - 0.57488*tau3^6 + 0.19383*tau3^8
      OLB <- function(tau3) -0.25 + 1.25*tau3^2
      
      plot(GP, add = TRUE)
      plot(GEV, add = TRUE)
      plot(GL, add = TRUE)
      plot(LN3, add = TRUE)
      plot(PE3, add = TRUE)
      plot(OLB, add = TRUE)
      points(c(0,1/3,log(9/8)/log(2),0),c(0,1/6,(16*log(2) - 10*log(3))/log(2),30/pi*atan(sqrt(2))-9), pch=17)
      text(rep(0.5,4),c(0.08,0.31,0.357,0.39),c('OLB','GP','GEV','GL'),srt=45)
      text(0.5,0.235,'PE3',srt=32)
      text(0.5,0.335,'LN3',srt=40)
      text(c(0,1/3,log(9/8)/log(2),0),c(0,1/6,(16*log(2) - 10*log(3))/log(2),30/pi*atan(sqrt(2))-9) + 0.015, c('U','E','G','N'))
    }
  }
}


  
