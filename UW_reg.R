UW<-expression(log(
  sigma*log(2)/(y*(-log(mu)))*(log(y)/log(mu))^(sigma-1)*2^(-(log(y)/log(mu))^(sigma))
)
)
m1UW<-D(UW,"mu")
s1UW<-D(UW,"sigma")
ms2UW<-D(m1UW,"sigma")

UW<-function (mu.link = "logit", sigma.link = "log") 
{
  tau<-.5
  mstats <- checklink("mu.link", "UW", substitute(mu.link), 
                      c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  dstats <- checklink("sigma.link", "UW", substitute(sigma.link), 
                      c("inverse", "log", "identity", "own"))
  structure(list(family = c("UW", "UWeibull"), 
                 parameters = list(mu = TRUE, sigma = TRUE), 
                 nopar = 2, 
                 type = "Continuous", 
                 mu.link = as.character(substitute(mu.link)), 
                 sigma.link = as.character(substitute(sigma.link)), 
                 
                 mu.linkfun = mstats$linkfun, 
                 sigma.linkfun = dstats$linkfun, 
                 
                 mu.linkinv = mstats$linkinv, 
                 sigma.linkinv = dstats$linkinv, 
                 
                 mu.dr = mstats$mu.eta, 
                 sigma.dr = dstats$mu.eta, 
                 
                 dldm = function(y, mu, sigma) {#ok
                   tau<-.5
                   dldm <- eval(m1UW)
                   dldm
                 }, 
                 d2ldm2 = function(y,mu, sigma) {
                   tau<-.5
                   dldm <- eval(m1UW)
                   d2ldm2 <- -dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15) 
                   d2ldm2
                 }, 
                 dldd = function(y, mu, sigma) {#ok
                   tau<-.5
                   dldd <- eval(s1UW)
                   dldd
                 },
                 d2ldd2 = function(y,mu, sigma) {
                   tau<-.5
                   dldd <- eval(s1UW)
                   d2ldd2 = -dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
                   d2ldd2
                 },
                 d2ldmdd = function(y,mu, sigma) {
                   tau<-.5
                   dldm <- eval(m1UW)
                   dldd <- eval(s1UW)
                   d2ldmdd = -(dldm * dldd)
                   d2ldmdd<-ifelse(is.na(d2ldmdd)==TRUE,0,d2ldmdd)
                   d2ldmdd  
                 }, 
                 G.dev.incr = function(y, mu, sigma, w, ...) -2 * log(dUW(y=y, mu=mu, sigma=sigma)), 
                 rqres = expression(
                   rqres(pfun = "pUW",  type = "Continuous", y = y, mu = mu, sigma = sigma)
                 ),
                 
                 mu.initial = expression(     mu <- rep(0.5#median(y)
                                                        ,length(y))),   
                 sigma.initial = expression(sigma<- rep(0.5, length(y))),
                 mu.valid = function(mu) all(mu > 0 & mu < 1), 
                 sigma.valid = function(sigma)  all(sigma > 0),
                 y.valid = function(y) all(y > 0 &  y < 1)
  ), 
  class = c("gamlss.family", "family"))
}


# density function
dUW<-function(y, mu = 0.7, sigma = 0.5, log = FALSE)
{
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", "")) 
  if (any(y <= 0) | any(y >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  
  fy1 <- sigma*log(2)/(y)*(-log(mu))^(-1)*(log(y)/log(mu))^(sigma-1)*
    2^(-(log(y)/log(mu))^sigma)
  if(log==FALSE) fy<-fy1 else fy<-log(fy1)
  fy
}
#------------------------------------------------------------------------------------------ #ok
# cumulative distribution function
pUW<-function(q, mu = 0.7, sigma = 0.5, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(q <= 0) | any(q >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  cdf1<-  .5^((log(q)/log(mu))^sigma)
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  cdf
}
#------------------------------------------------------------------------------------------ #ok
# quantile function
qUW<-function(u,mu,sigma)
{
  g<-sigma
  q<- mu^(-log(u)/log(2))^(1/g)
  q
}

# inversion method for randon generation
rUW<-function(n,mu,sigma,a=0,b=1)
{
  u<- runif(n)
  tau<-.5
  y<- qUW(u,mu =mu, sigma =sigma)
  y
}
