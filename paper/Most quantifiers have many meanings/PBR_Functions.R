#Function to simulate logistic curves
## One response
get.resp <- function(p) rbinom(length(p), 1, p)

## Curve simulation
sim.curve <- function(boundary, perc.noise, resp.error, bound.noise){
  x <- seq(0, 1, .01)
  xmat <- matrix(rep(x, 100000), nrow = length(x))
  res <- apply(xmat, 1, resp
               , noise = perc.noise
               , boundary = boundary
               , resp.error = resp.error
               , bound.noise = bound.noise
  )
  return(colMeans(res))
}

## Response rule: respond TRUE if perceived quantity is greater than perceived threshold. After this is decided the response errors kick in. 
resp <- function(x, boundary, noise = 0, resp.error = 0.05, bound.noise = 0.2){
  x.perc <- rnorm(length(x), x, noise)
  bound.perc <- rnorm(length(x), boundary, bound.noise)
  p <- ifelse(x.perc > bound.perc, 1, 0)
  p.noise <- msm::rtnorm(length(p)
                         , p
                         , ifelse(resp.error == 0, .00001, resp.error)
                         , lower = 0, upper = 1)
  get.resp(p.noise)
}

#Functions to run the model
getDat <- function(dat){
  
  D <- 5
  N <- nrow(dat)
  I <- length(unique(dat$workerid))
  
  y <- dat$resp
  cperc <- (dat$percent - 50) / 100
  sub <- as.numeric(factor(dat$workerid))
  few <- ifelse(dat$qq == 2, 1, 0)
  fewer <- ifelse(dat$qq == 3, 1, 0)
  many <- ifelse(dat$qq == 4, 1, 0)
  more <- ifelse(dat$qq == 5, 1, 0)
  most <- ifelse(dat$qq == 6, 1, 0)
  above <- as.numeric(cperc > 0)
  
  
  standat <- list(D = D, N = N, I = I, y = y, cperc = cperc, sub = sub,
                  few = few, fewer = fewer, many = many, more = more, most = most
                  , above = above)  
  return(standat)
}

myRunner <- function(standat, iter = 1000, warmup = 400, mod, nchains = 4, ...){
  inits <- function(...){
    list(delta = runif(standat$D, -.5, .5)
         , sigma2 = runif(standat$D, 1, 1.5)
         , beta = matrix(runif(standat$I * standat$D, -1, 1), nrow = standat$I, ncol = standat$D)
         , nu = runif(standat$D, .1, 1)
         , sigma2alpha = runif(standat$D, 1, 1.5)
         , alpha = matrix(runif(standat$I * standat$D, .1, 1), nrow = standat$I, ncol = standat$D))
  }
  fit <- sampling(mod, verbose=T,
                  data = standat, 
                  iter = iter,
                  warmup = warmup,
                  chains = nchains,
                  cores = nchains,
                  # init = lapply(1:4, inits),
                  pars = c("delta", "nu", "beta", "alpha", "gamma", "sigma2", "sigma2alpha")
                  , ...)
  return(fit)}

#Function to calculate response curve based on parameters
curve.calc <- function(x, p = cperc){
  a <- (p - x[1])/x[2]
  x[3]  + (1 - 2* x[3]) * exp(a)/(exp(a) + 1)
}  

#Plot model fit (for Fig 2)
ind.obs.fig <- function(mat, quant.col, quant.title){
  matplot(t(mat)
          , type = "l"
          , col = adjustcolor(1, .1)
          , lty = 1
          , xaxt = "n"
          , xlim = c(.5, 10.5)
          , ylab = "Proportion 'true'"
          , xlab = "Presented Percentage"
          , frame.plot = F
          , yaxt = "n")
  axis(2, c(0, .5, 1))
  axis(1, at = seq(.5, 10.5, 1)
       , labels = seq(0, 100, 10))
  abline(v = 5.5, col = adjustcolor(1, .4))
  abline(h = .5, col = adjustcolor(1, .4))
  avg <- colMeans(mat, na.rm = T)
  points(1:10, avg, cex = 1.5, col = quant.col, pch = 19)
  title(quant.title, adj = 0, line = 1)
}

#recalculate the mean thresholds
thrPer <- function(threshold){
  threshold*100 + 50
}

#calculate different between thresholds
thrDiff <- function(x){
  x[2] - x[1]
}


meanDiff <- function(baseQ, compQ) {
  mdiff <- baseQ
for(i in 1:ncol(baseQ)){
  mdiff[,i] <- apply(cbind(compQ[,i], baseQ[,i]), 1, thrDiff)
}
  return(mdiff)
}

#Plot correlations between parameters
make.cor.plot <- function(mat, ...){
  M <- cor(mat)
  pval <- corr.test(mat, adjust = "bonferroni") #add for significance labels
  corrplot(M, p.mat = pval$p, insig = 'label_sig', sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9, diag=FALSE, type = "upper", mar = c(0,0,1,0), ...)
  return(pval)
}

#Check the outliers
outliersComp <- function(mat, tit){
  mat2 <- NULL
  m1 <- lm(mat[,1]~mat[,2])
  m2 <- lm(mat[,1]~mat[,3])
  m3 <- lm(mat[,2]~mat[,3])
  mat2 <- cbind(ols_plot_cooksd_bar(m1)$data$color,ols_plot_cooksd_bar(m2)$data$color,ols_plot_cooksd_bar(m3)$data$color)
  mat2 <- as.data.frame(mat2)
  out <-  ifelse(mat2=="outlier", 1, 0)
  out <- as.data.frame(out)
  out$out_sum <- out[,1]+out[,2]+out[,3]
  out$out_fin <- as.factor(ifelse(out$out_sum==0,"normal","outlier"))
  pairs(mat, lower.panel = NULL, pch = 19, col = my_cols[out$out_fin], main = tit)
}