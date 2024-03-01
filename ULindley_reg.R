library(lamW)
library(gamlss)

##########################################################################
ULindley<-function (mu.link = "logit")
{
  mstats <- checklink("mu.link", "ULindley", substitute(mu.link),
                      c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  
  structure(list(family = c("ULindley", "Unit-Lindley"),
                 parameters = list(mu = TRUE),
                 nopar = 1,
                 type = "Continuous",
                 mu.link = as.character(substitute(mu.link)),
                 mu.linkfun = mstats$linkfun,
                 mu.linkinv = mstats$linkinv,
                 mu.dr = mstats$mu.eta,
                 dldm = function(y, mu) {
                   dldm <- (y/(1-y)*(1-mu)-mu*(mu+1))/(mu^2*(1-mu))
                   dldm
                 },
                 d2ldm2 = function(y,mu) {
                   dldm <- (y/(1-y)*(1-mu)-mu*(mu+1))/(mu^2*(1-mu))
                   d2ldm2 <- -dldm * dldm
                   d2ldm2 <- ifelse(d2ldm2 < -1e-15, d2ldm2,-1e-15)
                   d2ldm2
                 },
                 G.dev.incr = function(y, mu, w, ...) -2 * log(dULindley(y=y, mu=mu)),
                 rqres = expression(
                   rqres(pfun = "pULindley", type = "Continuous", y = y, mu = mu)
                 ),
                 mu.initial = expression(mu <- rep(mean(y),length(y))),
                 
                 mu.valid = function(mu) all(mu > 0 & mu < 1),
                 
                 y.valid = function(y) all(y > 0 & y < 1)
  ),
  class = c("gamlss.family", "family"))
}
#------------------------------------------------------------------------------------------
# density function
dULindley<-function(y, mu, log = FALSE){
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(y <= 0) | any(y >= 1)) stop(paste("x must be between 0 and 1", "\n", ""))
  fy1 <- ((1-mu)^2/(mu*(1-y)^3))*exp((-y*(1-mu))/(mu*(1-y)))
  if(log==FALSE) fy<-fy1 else fy<-log(fy1)
  fy
}

#------------------------------------------------------------------------------------------
# cumulative distribution function
pULindley<-function(q, mu, lower.tail = TRUE, log.p = FALSE)
{
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(q <= 0) | any(q >= 1)) stop(paste("x must be between 0 and 1", "\n", "")) 
  cdf1 <- 1-(1-(1-mu)*q/(q-1))*exp(-(1-mu)*q/(mu*(1-q)))
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<-1-cdf1
  if(log.p==FALSE) cdf<-cdf else cdf<-log(cdf)
  cdf
}

#------------------------------------------------------------------------------------------
# quantile function
qULindley<-function(u,mu){
  rep <- 1/mu
  q <- (rep+lambertWm1((rep)*(u-1)*exp(-(rep))))/
    (1+lambertWm1((rep)*(u-1)*exp(-(rep))))
  q
}

#------------------------------------------------------------------------------------------
# inversion method for random generation
rULindley<-function(n,mu)
{
  u<- runif(n)
  y<- qULindley(u=u,mu=mu)
  y
}
#------------------------------------------------------------------------------
# 
