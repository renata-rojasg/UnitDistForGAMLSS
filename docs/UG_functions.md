## Unit gamma distribution for GAMLSS

### Description

The function `UG()` defines the unit gamma distribution, a two-parameter distribution, for a `gamlss.family` object to be used in GAMLSS fitting using the function `gamlss()`. `UG()` has mean equal to the parameter `mu` and `sigma` as precision parameter. The functions `dUG`, `pUG`, `qUG` and `rUG` define the density, cumulative distribution function, quantile function and random generation for the unit gamma distribution.

### Usage
```r
UG(mu.link = "logit", sigma.link = "log")
dUG(x, mu = 0.5, sigma = 1, log = FALSE)
pUG(q, mu = 0.5, sigma = 1, lower.tail = TRUE, log.p = FALSE)
qUG(p, mu = 0.5, sigma = 1, lower.tail = TRUE, log.p = FALSE)
rUG(n, mu = 0.5, sigma = 1)
```

### Details

The unit gamma distribution distribution, UG, is given as:

$$
f(y;\mu,\sigma^2)=\left[\frac{\mu^{1/\sigma}}{1-\mu^{1/\sigma}}\right]^\sigma\,\frac{1}{\Gamma(\sigma)}\,y^{\,\mu^{1/\sigma}/(1-\mu^{1/\sigma})-1}\,[-\log(y)]^{\sigma-1},
$$

for \( 0 < y < 1 \), \( 0 < \mu < 1 \), and \( \sigma > 0 \); see Mousa et al. (2016).

### Value

`UG()` returns a `gamlss.family` object which can be used to fit a UG distribution in the `gamlss()` function.

### Author

Renata Rojas Guerra

### Reference

Mousa, A. M., El-Sheikh, A. A., & Abdel-Fattah, M. A. (2016). A gamma regression for bounded continuous variables. Advances and Applications in Statistics, 49, 305–326.