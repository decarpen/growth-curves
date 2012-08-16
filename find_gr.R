####################################
# FIND GROWTH RATES                #
#----------------------------------#
# DANIELLE E CARPENTER             #
# decarpen@princeton.edu		   #
####################################

#This script contains a function that can be used to find the growth rate of
#a set of points grown on the Biotek plate reader.

#x = a vector of data points, e.g. c(0.1, 0.1, 0.2, 0.25, 0.3, ...)
#t = a vector of time points, e.g. c(0, 15, 30, 45, 60, ...)
#plottitle = the title to go on the plot of the fit
#int = the number of time steps that should be used to fit the line
#r2.cutoff = how stringent the fit needs to be --> 1 is the max.

findgr = function(x, t, plottitle, int=8, r2.cutoff=0.997) {
  x = as.numeric(x)
  n = length(x)
  mat = NULL
  
  #are x and t the same length?
  if (length(x) != length(t)) {
    cat("Error: Your data and time are not the same length.\n")
    stop()
  }
  
  #is the line basically flat?
  fit = lm(x~t)
  m = abs(coefficients(fit)[[2]])
  if (m < 0.00001) {
    max=c(0,0,0,NA)
    lag=NA
    
    plot(t, log(x), pch=20, xlab="time", ylab="ln(OD600)", main=plottitle)
    mtext("no growth", side=3, line=-1, at=0, cex=0.8, adj=0)
  }
  
  #if not, find a slope
  else{
  	x[which(x <= 0)] = 0.001 	#transform values < 0
  	
    x = log(x)
    for (i in 1:(n-int)) {
      fit = lm(x[i:(i+int)]~t[i:(i+int)])		#linear regression on log transformed data.
      m = coefficients(fit)[[2]]
      b = coefficients(fit)[[1]]
      r2 = summary(fit)$r.squared
      mat = rbind(mat, c(i, b, m, r2))
    }
  mat = mat[which(mat[,4] > r2.cutoff),]  #only include slopes greater than the R2 cutoff.
  max = mat[which.max(mat[,3]),]
  par(las=1, mar=c(5, 4, 4, 4) + 0.1)
  plot(t,x, type="n", pch=20, xlab="time", ylab="ln(OD600)", main=plottitle)
  usr.old = par("usr")
    
  #how long is this in exponential growth?
  fit.line = sapply(t, function(x) max[3]*x+max[2])
  resid = fit.line-x
  residper = abs(resid/fit.line)
  resid.mat = rbind(t, fit.line, x, resid, residper)
  resid.mat = resid.mat[,which(resid.mat[5,] < 0.05)]
  lag = resid.mat[1,1]
  mtext(paste("lag =",round(lag,2)),side=3, line=-3, at=0, cex=0.8, adj=0)
  time.in.exp = resid.mat[1,ncol(resid.mat)]-resid.mat[1,1]
  abline(v=lag, col="cadet blue", lty=2)
  abline(v=resid.mat[1,ncol(resid.mat)], col="cadet blue", lty=2)
      
  #dx = diff(x)/(t[2]-t[1])
  #par(usr=c(par("usr")[1:2],min(dx)*1.05, max(dx)*1.05))
  #points(t[1:(length(t)-1)],dx, pch=18, type="o", col="dark grey", lty=1)
  #axis(4, col.axis="dark grey", col.ticks="dark grey")
  #mtext("delta(x)/delta(t)", side=4, line=3, col="dark grey", las=3)
    
	#plot
  par(usr=usr.old)
  points(t,x,pch=20)
  abline(lm(x[max[1]:(max[1]+int-1)]~t[max[1]:(max[1]+int-1)]), col="red", lty=2, lwd=2)
  points(t[max[1]:(max[1]+int-1)], x[max[1]:(max[1]+int-1)], col="red")
  mtext(paste("m =",round(max[3],3)), side=3, line=-1, at=0, cex=0.8, adj=0)
  mtext(paste("r2 =",round(max[4],4)), side=3, line=-2, at=0, cex=0.8, adj=0)
  }
  return(c("m"=max[3], "r2"=max[4], "lag"=lag))
}