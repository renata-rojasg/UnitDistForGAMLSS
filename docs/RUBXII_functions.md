## Reflected unit Burr XII distribution for GAMLSS

### Description

The function `RUBXII()` defines the Reflected Unit Burr XII distribution, a two-parameter distribution, for a `gamlss.family` object to be used in GAMLSS fitting using the function `gamlss()`. `RUBXII()` has the \(\tau \times 100\)-th quantile equal to the parameter `mu` and `sigma` as shape parameter. In this implementation, \(\tau = 0.5\) is used to model the median. The functions `dRUBXII`, `pRUBXII`, `qRUBXII` and `rRUBXII` define the density, cumulative distribution function, quantile function and random generation for the Reflected Unit Burr XII distribution.

### Usage
```r
RUBXII(mu.link = "logit", sigma.link = "log")
dRUBXII(x, mu = 0.5, sigma = 1, log = FALSE)
pRUBXII(q, mu = 0.5, sigma = 1, lower.tail = TRUE, log.p = FALSE)
qRUBXII(p, mu = 0.5, sigma = 1, lower.tail = TRUE, log.p = FALSE)
rRUBXII(n, mu = 0.5, sigma = 1)
```

### Details

The Reflected unit Burr XII distribution, RUBXII, is given as:

$$
f(y; \mu, \sigma) =\frac{\log(1 - \tau)^{-\sigma} \log^{\sigma - 1}(1 - y)^{-1}}
    {(1 - y) \log[1 + \log^\sigma (1 - \mu)^{-1}]}
    [1 + \log^\sigma(1 - y)^{-1}]^{\frac{\log(1 - \tau)}{\log[1 + \log^\sigma(1 - \mu)^{-1}]}-1},
$$

for \(\tau=0.5,\) \( 0 < y < 1 \), \( 0 < \mu < 1 \), and \( \sigma > 0 \); see Ribeiro et al. (2021).

### Value

`RUBXII()` returns a `gamlss.family` object which can be used to fit a Reflected unit Burr XII distribution in the `gamlss()` function.

### Author

Renata Rojas Guerra

### Reference

Ribeiro, T.F.; Cordeiro, G.M.; Peña-Ramírez, F.A.; Guerra, R.R. A New Quantile Regression for the COVID-19 Mortality Rates in the United States. Comput. Appl. Math. 2021, 40, 255.