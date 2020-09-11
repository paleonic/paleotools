#####################################################
### Intercept Comparisons @ x != 0 Using a t-test ###
#####################################################

##DESCRIPTION
#   This function is meant to make comparisons between line intercepts where comparisons are to made at an x-value other than 0
#   (i.e., not the true intercept). Lines are not adjusted to be equal.

#Requires the 'smatr' package
library(smatr)

##ARGUMENTS
#   xvar/yvar - numeric vectors of equal length, represent x,y data used in regression analyses. NAs are accounted for.
#               Function does not log-transform the datat.
#   grp1/grp2 - numeric vectors specifying which observations belong to groups of interest
#   x.min     - x-value at which comparison should be made.

##VALUES
#Returns a 2x8 matrix with:
#   Int.grp1/Int.grp2           - True intercept values for the two groups being compared
#   Int(adj).grp1/Int(adj).grp2 - Adjusted intercepts based on specified 'x.min'
#   t-stat                      - t-statistic for comparison
#   P-value                     - p-value for intercept comparison
#Rows correspond to comparison made using model II (SMA) and model I (OLS) lin-fitting approaches.


int.comp <- function(xvar,yvar,grp1,grp2,x.min) {
  nas.a <- union(which(is.na(xvar[grp1])),which(is.na(yvar[grp1])))
  nas.b <- union(which(is.na(xvar[grp2])),which(is.na(yvar[grp2])))
  xvar <- xvar[union(grp1,grp2)]
  yvar <- yvar[union(grp1,grp2)]
  nas <- union(which(is.na(xvar)),which(is.na(yvar)))
  if(length(nas)>0) {
    xvar <- xvar[-nas]
    yvar <- yvar[-nas]
  }
  grp1 <- (1:(length(grp1)-length(nas.a)))
  grp2 <- ((max(grp1)+1):length(xvar))
  #group1
  n.a <- length(xvar[grp1])
  mean.x.a <- mean(xvar[grp1])
  sum.X.a <- sum(xvar[grp1])
  sum.X2.a <- sum(xvar[grp1]^2)
  sum.x2.a <- sum.X2.a-(sum.X.a^2)/n.a
  mean.y.a <- mean(yvar[grp1])
  sum.Y.a <- sum(yvar[grp1])
  sum.Y2.a <- sum(yvar[grp1]^2)
  sum.y2.a <- sum.Y2.a-(sum.Y.a^2)/n.a
  sum.XY.a <- sum(xvar[grp1]*yvar[grp1])
  sum.xy.a <- sum.XY.a-(sum.X.a*sum.Y.a)/n.a
  ols.m.a <- sum.xy.a/sum.x2.a
  ols.b.a <- coef(lm(yvar[grp1]~xvar[grp1]))[1]
  sma.m.a <- coef(sma(yvar[grp1]~xvar[grp1]))[2]
  sma.b.a <- coef(sma(yvar[grp1]~xvar[grp1]))[1]
  x.cor.a <- xvar[grp1]-x.min #intercept correction (cor)
  sma.b.cor.a <- coef(sma(yvar[grp1]~x.cor.a))[1]
  ols.b.cor.a <- coef(sma(yvar[grp1]~x.cor.a,method = "OLS"))[1]
  SS.a <- sum.y2.a-(sum.xy.a^2)/sum.x2.a
  SS.sma.a <- sum(residuals(sma(yvar[grp1]~xvar[grp1]))^2)
  SS.ols.a <- sum(residuals(sma(yvar[grp1]~xvar[grp1],method = "OLS"))^2)
  DF.a <- n.a-2
  #group2
  n.b <- length(xvar[grp2])
  mean.x.b <- mean(xvar[grp2])
  sum.X.b <- sum(xvar[grp2])
  sum.X2.b <- sum(xvar[grp2]^2)
  sum.x2.b <- sum.X2.b-(sum.X.b^2)/n.b
  mean.y.b <- mean(yvar[grp2])
  sum.Y.b <- sum(yvar[grp2])
  sum.Y2.b <- sum(yvar[grp2]^2)
  sum.y2.b <- sum.Y2.b-(sum.Y.b^2)/n.b
  sum.XY.b <- sum(xvar[grp2]*yvar[grp2])
  sum.xy.b <- sum.XY.b-(sum.X.b*sum.Y.b)/n.b
  ols.m.b <- sum.xy.b/sum.x2.b
  ols.b.b <- coef(lm(yvar[grp2]~xvar[grp2]))[1]
  sma.m.b <- coef(sma(yvar[grp2]~xvar[grp2]))[2]
  sma.b.b <- coef(sma(yvar[grp2]~xvar[grp2]))[1]
  x.cor.b <- xvar[grp2]-x.min #intercept correction (cor)
  sma.b.cor.b <- coef(sma(yvar[grp2]~x.cor.b))[1]
  ols.b.cor.b <- coef(sma(yvar[grp2]~x.cor.b,method = "OLS"))[1]
  SS.b <- sum.y2.b-(sum.xy.b^2)/sum.x2.b
  SS.sma.b <- sum(residuals(sma(yvar[grp2]~xvar[grp2]))^2)
  SS.ols.b <- sum(residuals(sma(yvar[grp2]~xvar[grp2],method = "OLS"))^2)
  DF.b <- n.b-2
  #pooled regression
  SS.pool <- SS.a+SS.b
  SS.sma.pool <- SS.sma.a+SS.sma.b
  SS.ols.pool <- SS.ols.a+SS.ols.b
  DF.pool <- DF.a+DF.b
  #Common regression
  sum.x2.com <- sum.x2.a+sum.x2.b
  sum.xy.com <- sum.xy.a+sum.xy.b
  sum.y2.com <- sum.y2.a+sum.y2.b
  ols.m.com <- sum.xy.com/sum.x2.com
  SS.com <- sum.y2.com-(sum.xy.com^2)/sum.x2.com
  DF.com <- (n.a+n.b)-2-1
  #total regression
  n.tot <- length(xvar[union(grp1,grp2)])
  sum.X.tot <- sum(xvar[union(grp1,grp2)])
  sum.X2.tot <- sum(xvar[union(grp1,grp2)]^2)
  sum.x2.tot <- sum.X2.tot-(sum.X.tot^2)/n.tot
  sum.Y.tot <- sum(yvar[union(grp1,grp2)])
  sum.Y2.tot <- sum(yvar[union(grp1,grp2)]^2)
  sum.y2.tot <- sum.Y2.tot-(sum.Y.tot^2)/n.tot
  sum.XY.tot <- sum(xvar[union(grp1,grp2)]*yvar[union(grp1,grp2)])
  sum.xy.tot <- sum.XY.tot-(sum.X.tot*sum.Y.tot)/n.tot
  sma.m.tot <- coef(sma(yvar[union(grp1,grp2)]~xvar[union(grp1,grp2)]))[2]
  sma.b.tot <- coef(sma(yvar[union(grp1,grp2)]~xvar[union(grp1,grp2)]))[1]
  SS.tot <- sum.y2.tot-(sum.xy.tot^2)/sum.x2.tot
  SS.sma.tot <- sum(residuals(sma(yvar[union(grp1,grp2)]~xvar[union(grp1,grp2)]))^2)
  SS.ols.tot <- sum(residuals(sma(yvar[union(grp1,grp2)]~xvar[union(grp1,grp2)],method = "OLS"))^2)
  DF.tot <- n.tot-2
  #T-test
  #pooled residual mean square difference
  MS <- (SS.a+SS.b)/(DF.a+DF.b)
  MS.sma <- (SS.sma.a+SS.sma.b)/(DF.a+DF.b)
  MS.ols <- (SS.ols.a+SS.ols.b)/(DF.a+DF.b)
  #standard error of difference
  SE <- sqrt((MS/sum.x2.a)+(MS/sum.x2.b))
  SE.sma <- sqrt(MS.sma*((1/n.a)+(1/n.b)+(((mean.x.a)^2)/sum.x2.a)+((mean.x.b)^2)/sum.x2.b))
  SE.ols <- sqrt(MS.ols*((1/n.a)+(1/n.b)+(((mean.x.a)^2)/sum.x2.a)+((mean.x.b)^2)/sum.x2.b))
  t.sma.slp <- (sma.m.a-sma.m.b)/SE.sma
  t.ols.slp <- (ols.m.a-ols.m.b)/SE.ols
  if(t.sma.slp>0) {
    t.sma.slp <- t.sma.slp*-1
  }
  if(t.ols.slp>0) {
    t.ols.slp <- t.ols.slp*-1
  }
  tp.sma.slp <- 2*pt(t.sma.slp,df = DF.pool)
  tp.ols.slp <- 2*pt(t.ols.slp,df = DF.pool)
  #intercept t-test
  DF.com.t <- n.a+n.b-3
  MS.int <- SS.com/DF.com.t
  t.sma.int <- (sma.b.a-sma.b.b)/SE.sma
  t.ols.int <- (ols.b.a-ols.b.b)/SE.ols
  if(t.sma.int>0) {
    t.sma.int <- t.sma.int*-1
  }
  if(t.ols.int>0) {
    t.ols.int <- t.ols.int*-1
  }
  tp.sma.int <- 2*pt(t.sma.int,df = DF.com.t)
  tp.ols.int <- 2*pt(t.ols.int,df = DF.com.t)
  t.sma.cor.int <- (sma.b.cor.a-sma.b.cor.b)/SE.sma
  t.ols.cor.int <- (ols.b.cor.a-ols.b.cor.b)/SE.ols
  if(t.sma.cor.int>0) {
    t.sma.cor.int <- t.sma.cor.int*-1
  }
  if(t.ols.cor.int>0) {
    t.ols.cor.int <- t.ols.cor.int*-1
  }
  tp.sma.cor.int <- 2*pt(t.sma.cor.int,df = DF.com.t)
  tp.ols.cor.int <- 2*pt(t.ols.cor.int,df = DF.com.t)
  t.mat <- matrix(,2,8)
  colnames(t.mat) <- c("Int.grp1","Int.grp2","t-stat","P-value","Int(adj).grp1","Int(adj).grp2","t-stat","P-value")
  rownames(t.mat) <- c("SMA","OLS")
  t.mat[1,] <- c(sma.b.a,sma.b.b,t.sma.int,tp.sma.int,sma.b.cor.a,sma.b.cor.b,t.sma.cor.int,tp.sma.cor.int)
  t.mat[2,] <- c(ols.b.a,ols.b.b,t.ols.int,tp.ols.int,ols.b.cor.a,ols.b.cor.b,t.ols.cor.int,tp.ols.cor.int)
  return(t.mat)
}