##################################################
### Linear Statistics/Plot for Multiple Groups ###
##################################################

##DESCRIPTION
#   This function returns a summarized table of line-fitting results for the overall sample, as well as for specified groups.
#   A bivariate plot is also returned dividing the sample by specified group. Analyses can be done usin three line-fitting 
#   models: Ordinary Least Squares, Standardized Major Axis (AKA Reduced Major Axis), and Major Axis. Requires 'smatr' library.

#Requires the 'smatr' package
library(smatr)

##ARGUMENTS
#   xvar/yvar - numeric vectors of equal length, representing the x and y data to be used in analysis. NAs are accounted for.
#               Function log-transforms the data.
#   xlab/ylab - x and y labels, as per the 'plot' function.
#   groups    - vector of class 'character' and of equal length to xvar/yvar that specifies the group to which each
#               observation belongs.
#   model     - line-fitting model to be used, model=c('SMA','OLS','MA').

##VALUES
#Returns a kx8 matrix (where 'k' is the number of groups) with:
#   N                       - sample size used in the specific analysis
#   m                       - slope of the line
#   m.lower95CI/m.upper95CI - 95% confidence intervals of the slope
#   b                       - intercept (or elevation) of the line
#   b.lower95CI/b.upper95CI - 95% confidence intervals of the intercept
#   r-squared               - coefficient of determination

linear.stats<-function(xvar,yvar,xlab=NULL,ylab=NULL,groups,model) {
  grps<-levels(groups)
  col.grps<-rainbow(length(grps))
  x<-log10(xvar)
  y<-log10(yvar)
  linear.mat<-matrix(,length(grps)+1,8)
  rownames(linear.mat)<-c("All",grps)
  colnames(linear.mat)<-c('N','m','m.lower95CI','m.upper95CI','b','b.lower95CI','b.upper95CI','r.squared')
  plot(x,y,xlab=xlab,ylab=ylab,type='n',main=model)
  tot.line<-sma(y~x,method=model)
  abline(tot.line)
  linear.mat[1,]<-as.numeric(c(tot.line$n,tot.line$coef[[1]][2,],tot.line$coef[[1]][1,],tot.line$r2))
    for(i in 1:length(grps)) {
      grp<-which(groups==grps[i])
      points(x[grp],y[grp],col=col.grps[i],pch=19)
      line<-sma(y[grp]~x[grp],method=model)
      abline(line,col=col.grps[i])
      linear.mat[i+1,]<-as.numeric(c(line$n,line$coef[[1]][2,],line$coef[[1]][1,],line$r2))
    }
    legend("bottomright",legend=c("All",grps),pch=c(NA,rep(19,length(grps))),lty=1,col=c('black',col.grps),cex=0.75)
  return(linear.mat)
}