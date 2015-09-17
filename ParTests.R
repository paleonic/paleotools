#############################################################
### Standard Parametric Tests for Multivariate Statistics ###
#############################################################

##DESCRIPTION
#   This function runs standard tests on a table of numeric values meant for parametric approaches (e.g., PCA). Normality is tested using 
#   a shapiro test, heteroscedasticity a F-test, and linearity a Harvey-Collier Test. Requires the packages qdapTools and strucchange.

##ARGUMENTS
#   dataset  - a matrix with numeric values
#   adjust   - specifies the correction method by which p-values are adjusted for multiple comparisons. Default is "fdr", methods are those
#              available in p.adjust.methods

##VALUES
#Returns a list of length=3, including:
#   $normality          - a q x 5 matrix with mean, median, W-statistic, raw p-value, and p-value adjusted for multiple comparisons
#   $heteroscedasticity - q x q matrix with the p-values of the F-test. Upper triangle holds the raw p-values, lower tiangle holds
#                         adjusted p-values
#   $linearity          - q x q matrix with the p-values of the Harvey-Collier Test
#Also returns two level plots using package lattice, that colour-codes the p-values for the heteroscedasticity and linearity tests.

parametric.tests <- function(dataset, adjust = "fdr") {
  require(qdapTools)
  require(strucchange)
  require(lattice)
  q <- ncol(dataset)
  
  #Normality
  norm.res <- matrix(,q,5)
  rownames(norm.res) <- colnames(dataset)
  colnames(norm.res) <- c("mean","median","W.statistic","p.value","adj.p.value")
  norm.res[,1] <- apply(X = dataset,MARGIN = 2,FUN = mean,na.rm = TRUE)
  norm.res[,2] <- apply(X = dataset,MARGIN = 2,FUN = median,na.rm = TRUE)
  norm.test <- apply(X = dataset,MARGIN = 2,FUN = shapiro.test)
  norm.res[,3] <- as.numeric(matrix(unlist(norm.test),ncol = q)[1,])
  p.values <- as.numeric(matrix(unlist(norm.test),ncol = q)[2,])
  norm.res[,4] <- p.values
  norm.res[,5] <- p.adjust(p.values,method = adjust)
  
  #Heteroscedasticity
  heter.res <- v_outer(x = dataset,FUN = function(x,y) var.test(x,y)$p.value)
  diag(heter.res) <- NA
  heter.res[lower.tri(heter.res)] <- p.adjust(heter.res[upper.tri(heter.res)],method = adjust)
  
  #Linearity
  line.res <- v_outer(x = dataset,FUN = function(x,y) t.test(recresid(lm(y~x)),mu = 0)$p.value)
  diag(line.res) <- NA
  line.res[lower.tri(line.res)] <- p.adjust(line.res[upper.tri(line.res)],method = adjust)
  
  final.res <- list(normality = norm.res,heteroscedasticity = heter.res,linearity = line.res)
  rgb.palette <- colorRampPalette(c("blue", "yellow"), space = "rgb")
  print(levelplot(final.res$heteroscedasticity,main = "Heteroscedasticity p-values Plot",xlab="", ylab="",col.regions=rgb.palette(120),at = seq(0,1,0.05)))
  print(levelplot(final.res$linearity,main = "Linearity p-values Plot",xlab="", ylab="",at = seq(0,1,0.05),col.regions=rgb.palette(120)))
  return(final.res)
}
  
