#' @title Summary of lm and Anova
#' @description hdlinear can be used to carry out the summary of fitted linear models and print anova tables with a large amount of data in Rhadoop
#' @usage hdlinear(y, x, data)
#' @param y integer: the colnum number of a dependent variable
#' @param x vector for integer: the colnum numbers of independent variables
#' @param data an optional data frame containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called
#' @return This function returns a quantile of residuals, coefficients, residaul standard error, R.squared, F.statistic, and analysis of variance table.
#' @return t-test and f-test can be done with the results of a function.
#' @export
#' @examples
#' library(rhdfs)
#' hdfs.init()
#' library(rmr2)
#' rmr.options(backend = "hadoop")
#' hdlinear(y=5, x=c(7,11), bigdata)


hdlinear <- function(y, x, data) {

  options(digits = 5)
  dat <- data.frame(data[,y],data[,x])
  colnames(dat) = c('y',paste0(rep('x',length(x)),1:(length(x))))

  xnam <- paste0("x", 1:length(x))
  k = as.formula(paste("y ~ ", paste(xnam, collapse= "+")))

  X = model.matrix(k, dat)
  y = as.matrix(data[,y])
  X.index = to.dfs(cbind(1:nrow(X), X))

  reduce = function(., v) {
    XtX = Reduce("+",v[seq_along(v) %% 3 == 1])
    Xty = Reduce("+",v[seq_along(v) %% 3 == 2])
    yty = Reduce("+",v[seq_along(v) %% 3 == 0])
    keyval(1, list(XtX , Xty, yty))
  }


  ress = values( from.dfs( mapreduce(
    input = X.index,
    map = function(., Xi){
      yi = y[Xi[,1],]
      Xi = Xi[,-1]
      keyval(1, list(t(Xi) %*% Xi, t(Xi) %*% yi, t(yi) %*% yi))
    }, reduce = reduce, combine = T) ) )

  XtX <- ress[[1]]
  Xty <- ress[[2]]
  yty <- ress[[3]]

  beta.hat = solve(XtX, Xty)
  nn = XtX[1,1]
  df = c(ncol(X)-1, nn - ncol(X))
  y.hat <- X %*% beta.hat
  sigma.sq.hat <- sum((y - y.hat)^2) / df[2]
  beta.std.error <- sqrt( sigma.sq.hat * diag(solve(XtX)) )
  t.value <- beta.hat /( beta.std.error )
  p.value <- 2 * pt( abs(t.value), df = nn - 2, lower.tail = F)
  resid = y - X%*%beta.hat
  y.mean <- mean(y)
  SST <- sum((y - y.mean)^2)
  SSE <- sum((y - y.hat)^2)
  SSR <- SST - SSE
  ss <- c(SSR,SSE); ms <- c(SSR/df[1], SSE/df[2])
  res = list(n=nn,df=df, beta.hat = beta.hat, y.hat = y.hat, sigma.sq.hat=sigma.sq.hat, beta.std.error=beta.std.error, t.value=t.value, p.value=p.value, resid = resid, SST= SST, SSE = SSE, SSR = SSR, ss=ss, ms=ms)


  res1 <- quantile(res$resid)
  res2 <- data.frame(res$beta.hat, res$beta.std.error, res$t.value, res$p.value)
  colnames(res2) <- c('Estimate', 'Std. Error', 't value', 'Pr(>t)')
  res3 <- paste('Residaul standard error:',round(sqrt(res$sigma.sq.hat),digits= 4),'on', res$df[2], "degrees of freedom ", seq= ' ')
  res4 <- paste('Multiple R-squared : ', round(res$SSR/res$SST, digits=4),',', 'Adjusted R-squared:', round(1-res$ms[2]/(res$SST/(res$n-1)), digits=4), seq =' ')
  res5 <- paste('F-statistic:', round(res$ms[1]/res$ms[2], digits=4), 'on', res$df[1], 'and', res$df[2],'DF', 'p-value:',round(pf(res$ms[1]/res$ms[2], res$df[1] , res$df[2], lower.tail = FALSE), digits= 3),  seq= ' ')
  res6 <- data.frame(res$df, res$ss, res$ms, res$ms[1]/res$ms[2], pf(res$ms[1]/res$ms[2], res$df[1] , res$df[2], lower.tail = FALSE))
  res6[2,4:5] <- NA
  dimnames(res6) <- list(c('Regression', "Residuals"), c("Df","Sum Sq", "Mean Sq", "F value", "Pr(>F)"))
  list(Residuals = res1, Coefficients = res2, Residuals = res3, R.squared = res4, F.statistic = res5, anova = res6)
}

