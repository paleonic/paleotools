#################################
### Plot Minimum Convex Hulls ###
#################################

##DESCRIPTION
# This function plots the minimum convex hulls for a series of specified points.

##ARGUMENTS
#   xvar/yvar - numeric vectors of equal length, representing the x and y values on which to compute the convex hull.
#   type - character specifying the treatment of the shape, whether to “fill” or “outline”
#   col/lty/lwd   - as per the 'lines' function. Specifies the colour and line type for the convex hull.
#   alpha - if type = “fill”, numeric value specifying the level of opacity

##VALUES
#No values are returned

MCH <- function(xvar,yvar,type,col,lty = NULL,alpha = NULL) {
  mch<-chull(xvar,yvar)
  mch<-c(mch,mch[1])
  if(type == "outline") lines(xvar[mch],yvar[mch],col=col,lty=lty)
  if(type == "fill") polygon(xvar[mch],yvar[mch],col = adjustcolor(col,alpha.f = alpha),border = NA)
}