## Unit Weibull Distribution for GAMLSS

### Description

The functions `UW()` define the Unit Weibull distribution, a two-parameter distribution, for a `gamlss.family` object to be used in GAMLSS fitting using the function `gamlss()`. `UW()` has the \(\tau \times 100\)-th quantile equal to the parameter `mu` and `sigma` as a shape parameter. In this implementation, \(\tau = 0.5\) is used to model the median. The functions `dUW`, `pUW`, `qUW` and `rUW` define the density, cumulative distribution function, quantile function and random generation for the Unit Weibull distribution.

### Usage
```r
UW(mu.link = "logit", sigma.link = "log")
dUW(x, mu = 0.5, sigma = 1, log = FALSE)
pUW(q, mu = 0.5, sigma = 1, lower.tail = TRUE, log.p = FALSE)
qUW(p, mu = 0.5, sigma = 1, lower.tail = TRUE, log.p = FALSE)
rUW(n, mu = 0.5, sigma = 1)
```

### Details

The Unit Weibull distribution, UW, is given as:

$$
f(y; \mu, \sigma) = \frac{\sigma}{y}\left(\frac{\log \tau}{\log \mu}\right)
\left(\frac{\log y}{\log \mu}\right)^{\sigma - 1} 
\tau^{\left(\frac{\log y}{\log \mu}\right)^\sigma}
$$

for \( 0 < y < 1 \), \( 0 < \mu < 1 \), and \( \sigma > 0 \). In this implementation, \( \tau = 0.5 \); see Mazucheli et al. (2020).

### Author

Renata Rojas Guerra

### Reference

Mazucheli, J., Menezes, A. F. B., Fernandes, L. B., & de Oliveira, R. P. (2020). A quantile-based unit Weibull regression model for continuous bounded data. Statistical Papers, 61, 1549–1572.