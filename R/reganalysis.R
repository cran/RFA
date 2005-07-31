reganalysis <- function(data){
  ##This is a generic function to perform a regional
  ##frequency analysis via the Index Flood model.
  ##
  ##data : a list with each arguments correspond to
  ##       sample of a site within the region plus
  ##       an argument 'record': a vector with the record
  ##       length for each site.

  ##First of all, we should apply the test of [Hosking, J. R. M.
  ##and Wallis J. R. (1997) Regional Frequency Analysis Cambridge
  ##University Press] Ch. 4, p. 63

  n.site <- length(data) - 1
  name.site <- names(data)[-(n.site+1)]

  tt <- tktoplevel()
  tktitle(tt) <- "Index Flood and Sites Selection"
  spec.frm <- tkframe(tt, relief = 'groove', borderwidth = 2)

  frm1 <- tkframe(spec.frm, relief = 'groove', borderwidth = 2)
  tkpack(tklabel(frm1, text = 'Index Flood Selection'))
  
  index.flood <- tclVar("mean")
  prob <- tclVar(0.50)

  chkbut1 <- tkradiobutton(frm1)
  tkconfigure(chkbut1, text = "Sample Mean",
              variable = index.flood, value = "mean")
  chkbut2 <- tkradiobutton(frm1)
  tkconfigure(chkbut2, text = "Sample Median",
              variable = index.flood, value = "median")
  chkbut3 <- tkradiobutton(frm1)
  tkconfigure(chkbut3, text = "Theoretical Quantile with prob. of non exceedance",
              variable = index.flood, value = "quantile")
  slider <- tkentry(frm1, textvariable = prob, width = 5)

  onOK <- function() {
    res <- 1 + as.integer(tkcurselection(listbox))
    region <- NULL
    region <<- name.site[res]
    tkgrab.release(tt)
    tkdestroy(tt)
  }
  
  ok.button <- tkbutton(tt, text = 'OK',
                        command = onOK)

  frm2 <- tkframe(spec.frm, relief = "groove", borderwidth = 2)
  tkpack(tklabel(frm2, text = "Select sites within the region."))
  scr <- tkscrollbar(frm2, repeatinterval=5,
                   command=function(...)tkyview(listbox,...))
  listbox <- tklistbox(frm2, selectmode = "extended",
                       exportselection = "TRUE",
                       yscrollcommand=function(...)tkset(scr,...))
  tkpack(listbox, scr, side = 'left', fill = 'y')
  
  for (site in name.site) tkinsert(listbox, "end", site)

  viewLmom <- function(){
    res <- 1 + as.integer(tkcurselection(listbox))
    region <- name.site[res]
    lmom <- NULL
    n.site <- length(data) - 1
    for (site in 1:n.site){
      lmom <- rbind(lmom, samlmu(data[[site]]))
    }
    rownames(lmom) <- names(data)[1:n.site]
    plot(lmom[,3], lmom[,2], xlab = expression(paste('L-Skewness,  ',tau[3],sep='')),
         ylab = expression(paste('L-CV,  ',tau,sep='')))
    points(lmom[region,3], lmom[region,2], pch = 16)
    tkgrab.release(tt)
  }

  view.button <- tkbutton(frm2, text = 'View L-moments plot',
                          command = viewLmom)
  
  
  tkpack(chkbut1, anchor = "w")
  tkpack(chkbut2, anchor = "w")
  tkpack(chkbut3, slider)
  tkpack(listbox, view.button)
  
  tkpack(frm1, frm2, anchor = 'n', side = 'left', fill = 'y')
  tkpack(spec.frm, ok.button)
  tkwait.window(tt)
  
  index.flood <- tclvalue(index.flood)
  name.index.flood <- index.flood
  if (name.index.flood == "quantile")
    name.index.flood <- paste(name.index.flood,
                              " with probability of non exceedance ",
                              tclvalue(prob), sep="")
  
  switch(index.flood,'mean' = index.flood <- mean,
         'median' = index.flood <- median,
         'quantile' = index.flood <- function(sample){
              prob <- as.numeric(tclvalue(prob))
              loc <- min(sample)
              temp <- fitgpd(sample,loc,'mle',std.err = FALSE )$param
              scale <- temp[1]
              shape <- temp[2]
              rm(temp)
              return(qgpd(prob,loc,scale,shape))
            }
         )
  
    
  n.site <- length(region)
  
  temp <- data
  record <- NULL
  data <- list()
  
  for (i in 1:n.site){
    indice.site <- which( name.site == region[i]) 
    data <- c(data, list(temp[[indice.site]]))
    record <- c(record, temp$record[indice.site])
  }

  data <- c(data, list(record) )
  names(data) <- c(region,'record')
  n.site <- length(data) - 1
  name.site <- names(data)[-(n.site+1)]
  
  rm(temp)
    
  ##Computes several caracteristics of the specified region
  reg.car <- lmomreg(data[-(n.site+1)], index.flood)
  ##Fit a four-parameter Kappa distribution
  kappa <- kappalmom(reg.car$reg)
  
  size <- rep(NA,n.site)
  for (i in 1:n.site) size[i] <- length(data[[i]])

  tt2 <- tktoplevel()
  tktitle(tt2) <- 'Aim site and Optional Selection'
  spec.frm <- tkframe(tt2, relief = 'groove', borderwidth = 2)
  
  frm1 <- tkframe(spec.frm, relief = 'groove', borderwidth = 2)
  tkpack(tklabel(frm1, text = 'Aim site Selection'))
  scr <- tkscrollbar(frm1, repeatinterval=5,
                     command=function(...)tkyview(listbox,...))
  listbox <- tklistbox(frm1, selectmode = "extended",
                       exportselection = "TRUE",
                       yscrollcommand=function(...)tkset(scr,...))
  tkpack(listbox, scr, side = 'left', fill = 'y')
  for (site in region) tkinsert(listbox, "end", site) 
  
  frm2 <- tkframe(spec.frm, relief = 'groove', borderwidth = 2)
  tkpack(tklabel(frm2, text = 'Select Options'))

  flag.heter <- flag.lmomplot <- FALSE
  
  heterflag <- tclVar("FALSE")
  lmomplotflag <- tclVar("FALSE")
  chkbut1 <- tkradiobutton(frm2)
  tkconfigure(chkbut1, text = 'Compute the heterogeneity statistic',
              variable = heterflag, value = "TRUE")
  chkbut2 <- tkradiobutton(frm2)
  tkconfigure(chkbut2, text = 'Draw L-moments plot',
              variable = lmomplotflag, value = "TRUE")

  onOK <- function() {
    res <- 1 + as.integer(tkcurselection(listbox))
    aim.site <- NULL
    aim.site <<- name.site[res]
    tkgrab.release(tt2)
    tkdestroy(tt2)
  }
  
  ok.button <- tkbutton(tt2, text = 'OK',
                        command = onOK)
  
  tkpack(chkbut1, chkbut2, anchor = 'nw')
  tkpack(frm1, frm2, anchor = 'n', side = 'left', fill = 'y')
  tkpack(spec.frm, ok.button)
  tkwait.window(tt2)

  if (tclvalue(heterflag) == "TRUE" ){
    cat('Please wait while computing the heterogeneity statistic...\n')
    ##Evaluate the heterogeneity statistic
    heterogeneity(500, n.site, size, kappa, reg.car$V)
  }
  if (tclvalue(lmomplotflag) == "TRUE") lmomplot(reg.car)
 
  n.aimsite <- which( names(data) == aim.site )
  mu <- size[n.aimsite] / data$record[n.aimsite]
  
  par(ask = TRUE)
  param.regdist <- regdist(reg.car)
  param.sitedist <- locdist(param.regdist, mu, data[[n.aimsite]],
                            index.flood = index.flood)

  res <- list(reg.param=param.regdist, sit.param = param.sitedist,
              aimsite = aim.site, index.flood = name.index.flood,
              region = region)
  
  return(res)
  
}


    
  
