lRUBXII<-expression(log(
  sigma*log(1/(1-tau))*(log(1/(1-y)))^(sigma-1)*
    (1+(log(1/(1-y)))^sigma)^(log(1-tau)/log(1+(log(1/(1-mu)))^sigma)-1)/  
    ((1-y)*log(1+(log(1/(1-mu)))^sigma)) 
  )
  )
m1<-D(lRUBXII,"mu")
s1<-D(lRUBXII,"sigma")
ms2<-D(m1,"sigma")


RUBXII<-function (mu.link = "logit", sigma.link = "identity") 
{
  tau<-.5
  mstats <- checklink("mu.link", "RUBXII", substitute(mu.link), 
                      c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  dstats <- checklink("sigma.link", "RUBXII", substitute(sigma.link), 
                      c("inverse", "log", "identity", "own"))
  structure(list(family = c("RUBXII", "RUBXII"), 
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
                   dldm <- eval(m1)
                   dldm
                 }, 
                 d2ldm2 = function(y,mu, sigma) {
                   tau<-.5
                   d2ldm2 <- -eval(m1)^2
                   d2ldm2
                 }, 
                 dldd = function(y, mu, sigma) {#ok
                   tau<-.5
                   dldd <- eval(s1)
                   dldd
                 },
                 d2ldd2 = function(y,mu, sigma) {
                   tau<-.5
                   d2ldd2 <- -eval(s1)^2
                   d2ldd2
                 },
                 d2ldmdd = function(y,mu, sigma) {
                   tau<-.5
                   d2ldmdd <- eval(ms2)
                   d2ldmdd
                 }, 
                 G.dev.incr = function(y, mu, sigma, w, ...) -2 * log(dRUBXII(y=y, mu=mu, sigma=sigma)), 
                 rqres = expression(
                   rqres(pfun = "pRUBXII",  type = "Continuous", y = y, mu = mu, sigma = sigma)
                 ),
                 
                 mu.initial = expression(     mu <- rep(median(y),length(y))),   
                 sigma.initial = expression(sigma<- rep(.5, length(y))),
                 mu.valid = function(mu) all(mu > 0 & mu < 1), 
                 sigma.valid = function(sigma)  all(sigma > 0),
                 y.valid = function(y) all(y > 0 &  y < 1)
  ), 
  class = c("gamlss.family", "family"))
}


# density function
dRUBXII<-function(y, mu = 0.7, sigma = 0.5, log = FALSE)
{
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", "")) 
  if (any(y <= 0) | any(y >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  tau<-.5
  fy1 <-  sigma*log(1/(1-tau))*(log(1/(1-y)))^(sigma-1)*
    (1+(log(1/(1-y)))^sigma)^(log(1-tau)/log(1+(log(1/(1-mu)))^sigma)-1)/  
    ((1-y)*log(1+(log(1/(1-mu)))^sigma)) 
  if(log==FALSE) fy<-fy1 else fy<-log(fy1)
  fy
}
#------------------------------------------------------------------------------------------ #ok
# cumulative distribution function
pRUBXII<-function(q, mu = 0.7, sigma = 0.5, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(q <= 0) | any(q >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  tau<-.5
  cdf1<- 1-(1+(log(1/(1-q)))^sigma)^(log(1-tau)/log(1+(log(1/(1-mu)))^sigma))
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  cdf
}
#------------------------------------------------------------------------------------------ #ok
# quantile function
qRUBXII<-function(u,mu,sigma)
{
  tau<-.5
  q<- 1-exp(-((1-u)^(log(1+(log(1/(1-mu)))^sigma)/log(1-tau))-1)^(1/sigma)) #qf RUBXII qf
  q
}

# inversion method for randon generation
rRUBXII<-function(n,mu,sigma)
{
  u<- runif(n)
  tau<-.5
  y<- qRUBXII(u,mu =mu, sigma =sigma)
  y
}
RUBXII<-expression(log(
  sigma*log(1/(1-tau))*(log(1/(1-y)))^(sigma-1)*
    (1+(log(1/(1-y)))^sigma)^(log(1-tau)/log(1+(log(1/(1-mu)))^sigma)-1)/  
    ((1-y)*log(1+(log(1/(1-mu)))^sigma)) 
)
)
m1RUBXII<-D(RUBXII,"mu")
m2RUBXII<-D(m1RUBXII,"mu")
s1RUBXII<-D(RUBXII,"sigma")
s2RUBXII<-D(s1RUBXII,"sigma")
ms2RUBXII<-D(m1RUBXII,"sigma")


RUBXII<-function (mu.link = "logit", sigma.link = "identity") 
{
  tau<-.5
  mstats <- checklink("mu.link", "RUBXII", substitute(mu.link), 
                      c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  dstats <- checklink("sigma.link", "RUBXII", substitute(sigma.link), 
                      c("inverse", "log", "identity", "own"))
  structure(list(family = c("RUBXII", "RUBXII"), 
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
                   dldm <- eval(m1RUBXII)
                   dldm
                 }, 
                 d2ldm2 = function(y,mu, sigma) {
                   tau<-.5
                   dldm <- eval(m1RUBXII)
                   d2ldm2 <- -dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15) 
                   d2ldm2
                 }, 
                 dldd = function(y, mu, sigma) {#ok
                   tau<-.5
                   dldd <- eval(s1RUBXII)
                   dldd
                 },
                 d2ldd2 = function(y,mu, sigma) {
                   tau<-.5
                   dldd <- eval(s1RUBXII)
                   d2ldd2 = -dldd * dldd
                   d2ldd2 <- ifelse(d2ldd2 < -1e-15, d2ldd2,-1e-15)  
                   d2ldd2
                 },
                 d2ldmdd = function(y,mu, sigma) {
                   tau<-.5
                   dldm <- eval(m1RUBXII)
                   dldd <- eval(s1RUBXII)
                   d2ldmdd = -(dldm * dldd)
                   d2ldmdd<-ifelse(is.na(d2ldmdd)==TRUE,0,d2ldmdd)
                   d2ldmdd  
                 }, 
                 G.dev.incr = function(y, mu, sigma, w, ...) -2 * log(dRUBXII(y=y, mu=mu, sigma=sigma)), 
                 rqres = expression(
                   rqres(pfun = "pRUBXII",  type = "Continuous", y = y, mu = mu, sigma = sigma)
                 ),
                 
                 mu.initial = expression(     mu <- rep(median(y),length(y))),   
                 sigma.initial = expression(sigma<- rep(3, length(y))),
                 mu.valid = function(mu) all(mu > 0 & mu < 1), 
                 sigma.valid = function(sigma)  all(sigma > 0),
                 y.valid = function(y) all(y > 0 &  y < 1)
  ), 
  class = c("gamlss.family", "family"))
}


# density function
dRUBXII<-function(y, mu = 0.7, sigma = 0.5, log = FALSE)
{
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma <= 0)) stop(paste("sigma must be positive", "\n", "")) 
  if (any(y <= 0) | any(y >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  tau<-.5
  fy1 <-  sigma*log(1/(1-tau))*(log(1/(1-y)))^(sigma-1)*
    (1+(log(1/(1-y)))^sigma)^(log(1-tau)/log(1+(log(1/(1-mu)))^sigma)-1)/  
    ((1-y)*log(1+(log(1/(1-mu)))^sigma)) 
  if(log==FALSE) fy<-fy1 else fy<-log(fy1)
  fy
}
#------------------------------------------------------------------------------------------ #ok
# cumulative distribution function
pRUBXII<-function(q, mu = 0.7, sigma = 0.5, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(sigma < 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(q <= 0) | any(q >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  tau<-.5
  cdf1<- 1-(1+(log(1/(1-q)))^sigma)^(log(1-tau)/log(1+(log(1/(1-mu)))^sigma))
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  cdf
}
#------------------------------------------------------------------------------------------ #ok
# quantile function
qRUBXII<-function(u,mu,sigma)
{
  tau<-.5
  q<- 1-exp(-((1-u)^(log(1+(log(1/(1-mu)))^sigma)/log(1-tau))-1)^(1/sigma)) #qf RUBXII qf
  q
}

# inversion method for randon generation
rRUBXII<-function(n,mu,sigma,a=0,b=1)
{
  u<- runif(n)
  tau<-.5
  y<- qRUBXII(u,mu =mu, sigma =sigma)
  y
}

#
# set.seed(100)
# n<-1000
# X<-runif(n)
# ddd<-make.link("logit")
# y<-rRUBXII(n,ddd$linkinv(-1.6+1.25*X),2.3)
# 
# saida<-gamlss(y~X, family="RUBXII")
# sai<-summary(saida)
# round(sai,4)
# 
# 
# 

