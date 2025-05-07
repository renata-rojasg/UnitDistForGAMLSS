## Kumaraswamy Distribution for GAMLSS

### Description

The functions `KW()` define the Kumaraswamy distribution, a two-parameter distribution, for a `gamlss.family` object to be used in GAMLSS fitting using the function `gamlss()`. `KW()` has median equal to the parameter `mu` and `sigma` as dispersion parameter. The functions `dKW`, `pKW`, `qKW` and `rKW` define the density, cumulative distribution function, quantile function and random generation for the Kumaraswamy distribution.

### Usage
```r
KW(mu.link = "logit", sigma.link = "log")
dKW(x, mu = 0.5, sigma = 1, log = FALSE)
pKW(q, mu = 0.5, sigma = 1, lower.tail = TRUE, log.p = FALSE)
qKW(p, mu = 0.5, sigma = 1, lower.tail = TRUE, log.p = FALSE)
rKW(n, mu = 0.5, sigma = 1)
```

### Details

The Kumaraswamy distribution, KW, is given as:

$$
f(y; \mu, \sigma) = \frac{\log(0.5)}{\sigma^2 \log(1 - \mu^{1/\sigma^2})}
y^{1/\sigma^2} (1 - y^{1/\sigma^2})^
{{\log(0.5)}/{\log(1 - \mu^{1/\sigma^2})} -1}
$$

for \( 0 < y < 1 \), \( 0 < \mu < 1 \), and \( \sigma > 0 \); see Mitnik (2013).

### Value

`KW()` returns a `gamlss.family` object which can be used to fit a Kumaraswamy distribution in the `gamlss()` function.

### Author

Renata Rojas Guerra

### Reference

 Mitnik, P.A.; Baek, S. The Kumaraswamy Distribution: Median-Dispersion Re-Parameterizations for Regression Modeling and Simulation-Based Estimation. Stat. Pap. 2013, 54, 177–192.