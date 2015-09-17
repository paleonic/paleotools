#################################
### Plot Minimum Convex Hulls ###
#################################

##DESCRIPTION
# This function plots the minimum convex hulls for a series of specified points.

##ARGUMENTS
#   xvar/yvar - numeric vectors of equal length, representing the x and y values on which to compute the convex hull.
#   col/lty/lwd   - as per the 'lines' function. Specifies the colour and line type for the convex hull.

##VALUES
#No values are returned

MCH<-function(xvar,yvar,col,lty=1,lwd=1){
  mch<-chull(xvar,yvar)
  mch<-c(mch,mch[1])
  lines(xvar[mch],yvar[mch],col=col,lty=lty,lwd=lwd)
}