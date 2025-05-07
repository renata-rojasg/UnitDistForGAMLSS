## Unit Lindley Distribution for GAMLSS

### Description

The function `ULindley()` defines the unit Lindley distribution, a one-parameter distribution, for a `gamlss.family` object to be used in GAMLSS fitting using the function `gamlss()`. `ULindley()` has mean equal to the parameter `mu`. The functions `dULindley`, `pULindley`, `qULindley` and `rULindley` define the density, cumulative distribution function, quantile function and random generation for the unit Lindley distribution.

### Usage
```r
ULindley(mu.link = "logit", sigma.link = "log")
dULindley(x, mu = 0.5, log = FALSE)
pULindley(q, mu = 0.5, lower.tail = TRUE, log.p = FALSE)
qULindley(p, mu = 0.5, lower.tail = TRUE, log.p = FALSE)
rULindley(n, mu = 0.5)
```

### Details

The unit Lindley distribution distribution, ULindley, is given as:

$$
f(y;\mu,\sigma^2)=\left[\frac{\mu^{1/\sigma}}{1-\mu^{1/\sigma}}\right]^\sigma\,\frac{1}{\Gamma(\sigma)}\,y^{\,\mu^{1/\sigma}/(1-\mu^{1/\sigma})-1}\,[-\log(y)]^{\sigma-1},
$$

for \( 0 < y < 1 \), and \( 0 < \mu < 1 \); see Mazucheli et al. (2019).

### Value

`ULindley()` returns a `gamlss.family` object which can be used to fit a ULindley distribution in the `gamlss()` function.

### Author

Renata Rojas Guerra

### Reference

Mazucheli, J., Menezes, A. F. B., & Chakraborty, S. (2019). On the one parameter unit-Lindley distribution and its associated regression model for proportion data. Journal of Applied Statistics, 46, 700–714.