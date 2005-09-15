reganalysis <- function(data){
  ##This is a generic function to perform a regional
  ##frequency analysis via the Index Flood model.
  ##
  ##data : a list with each arguments correspond to
  ##       sample of a site within the region plus
  ##       an argument 'record': a vector with the record
  ##       length for each site. The argument ``record'' is optional
  ##       if not present then it is computed.

  ##First of all, we should apply the test of [Hosking, J. R. M.
  ##and Wallis J. R. (1997) Regional Frequency Analysis Cambridge
  ##University Press] Ch. 4, p. 63

  ##First test, if the there is the ``record'' argument in data
  if (!any( names(data)  == 'record')){
    record <- NULL
    for (i in 1:length(data))
      record <- c(record, length(data[[i]]))
    data <- c(data, list(record = record))
  }

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
    select.sites <- NULL
    select.sites <<- name.site[res]
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
    select.sites <- name.site[res]
    lmom <- NULL
    n.site <- length(data) - 1
    for (site in 1:n.site){
      lmom <- rbind(lmom, samlmu(data[[site]]))
    }
    rownames(lmom) <- names(data)[1:n.site]
    plot(lmom[,3], lmom[,2], xlab = expression(paste('L-Skewness,  ',tau[3],sep='')),
         ylab = expression(paste('L-CV,  ',tau,sep='')))
    points(lmom[select.sites,3], lmom[select.sites,2], pch = 16)

    if ( tclvalue(identify.flag) == "TRUE")
      identify(lmom[,3], lmom[,2], labels = name.site)
    tkgrab.release(tt)
  }

  view.button <- tkbutton(frm2, text = 'View L-moments plot',
                          command = viewLmom)
  identify.flag <- tclVar("FALSE")
  identify.button <- tkradiobutton(frm2, text = 'Identify points',
                                   variable = identify.flag, value = "TRUE")
    
  
  tkpack(chkbut1, anchor = "w")
  tkpack(chkbut2, anchor = "w")
  tkpack(chkbut3, slider)
  tkpack(listbox, view.button)
  tkpack(identify.button, side = "left")
  
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
  
    
  n.site <- length(select.sites)
  
  temp <- data
  record <- NULL
  data <- list()
  
  for (i in 1:n.site){
    indice.site <- which( name.site == select.sites[i]) 
    data <- c(data, list(temp[[indice.site]]))
    record <- c(record, temp$record[indice.site])
  }

  data <- c(data, list(record) )
  names(data) <- c(select.sites,'record')
  n.site <- length(data) - 1
  name.site <- names(data)[-(n.site+1)]
  
  rm(temp)
    
  ##Computes several caracteristics of the specified region
  reg.car <- lmomreg(data[-(n.site+1)], index.flood)
  ##Fit a four-parameter Kappa distribution
  kappa <- kappalmom(reg.car$reg)
  
  size <- rep(NA,n.site)
  for (i in 1:n.site) size[i] <- length(data[[i]])

  tt <- tktoplevel()
  tktitle(tt) <- 'Aim site and Optional Selection'
  spec.frm <- tkframe(tt, relief = 'groove', borderwidth = 2)
  
  frm1 <- tkframe(spec.frm, relief = 'groove', borderwidth = 2)
  tkpack(tklabel(frm1, text = 'Aim site Selection'))
  scr <- tkscrollbar(frm1, repeatinterval=5,
                     command=function(...)tkyview(listbox,...))
  listbox <- tklistbox(frm1, selectmode = "extended",
                       exportselection = "TRUE",
                       yscrollcommand=function(...)tkset(scr,...))
  tkpack(listbox, scr, side = 'left', fill = 'y')
  for (site in select.sites) tkinsert(listbox, "end", site) 
  
  frm2 <- tkframe(spec.frm, relief = 'groove', borderwidth = 2)
  tkpack(tklabel(frm2, text = 'Select Options'))

  heter.flag <- lmomplot.flag <- "FALSE"
  
  heter.flag <- tclVar("FALSE")
  lmomplot.flag <- tclVar("FALSE")
  chkbut1 <- tkradiobutton(frm2)
  tkconfigure(chkbut1, text = 'Compute the heterogeneity statistic',
              variable = heter.flag, value = "TRUE")
  chkbut2 <- tkradiobutton(frm2)
  tkconfigure(chkbut2, text = 'Draw L-moments plot',
              variable = lmomplot.flag, value = "TRUE")

  onOK <- function() {
    res <- 1 + as.integer(tkcurselection(listbox))
    aim.site <- NULL
    aim.site <<- name.site[res]
    tkgrab.release(tt)
    tkdestroy(tt)
  }
  
  ok.button <- tkbutton(tt, text = 'OK',
                        command = onOK)
  
  tkpack(chkbut1, chkbut2, anchor = 'nw')
  tkpack(frm1, frm2, anchor = 'n', side = 'left', fill = 'y')
  tkpack(spec.frm, ok.button)
  tkwait.window(tt)

  if (tclvalue(heter.flag) == "TRUE" ){
    cat('Please wait while computing the heterogeneity statistic...\n')
    ##Evaluate the heterogeneity statistic
    heterogeneity(500, n.site, size, kappa, reg.car$V)
  }
  if (tclvalue(lmomplot.flag) == "TRUE") lmomplots(reg.car)
 
  n.aimsite <- which( names(data) == aim.site )
  mu <- size[n.aimsite] / data$record[n.aimsite]
  
  par(ask = TRUE)
  param.regdist <- regdist(reg.car)
  param.sitedist <- locdist(param.regdist, mu, data[[n.aimsite]],
                            index.flood = index.flood)

  res <- list(reg.param=param.regdist, sit.param = param.sitedist,
              aimsite = aim.site, index.flood = name.index.flood,
              region = select.sites)
  
  return(res)
  
}


    
  
