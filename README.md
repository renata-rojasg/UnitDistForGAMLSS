# 📦 Unit Distributions for GAMLSS

This repository provides R implementations of unit distributions to be used within the `gamlss.family` framework, allowing regression modeling for continuous variables bounded in the (0, 1) interval. These distributions are not yet included in the official `gamlss` R package. We focus on one and two-parameter distributions, with parameterizations based on the mean or median (denoted as &mu;) and shape or precision/dispersion parameters (denoted as &sigma;).

## 📚 List of Distributions
---

### Kumaraswamy Distribution

The functions `KW()` define the two-parameter Kumaraswamy distribution for use with `gamlss()`, where `mu` corresponds to the median and `sigma` to a dispersion parameter. Supporting functions include: `dKW`, `pKW`, `qKW`, and `rKW`.

📌 See the R code: [`KW_reg.R`](KW_reg.R)  
📄 See detailed documentation: [docs/KW_functions.md](docs/KW_functions.md)

> **Reference:**  
> Mitnik, P.A.; Baek, S. The Kumaraswamy Distribution: Median-Dispersion Re-Parameterizations for Regression Modeling and Simulation-Based Estimation. Stat. Pap. 2013, 54, 177–192.



### Reflected Burr XII Distribution

The `RUBXII()` function implements a reflected Burr XII distribution, parameterized by median (`mu`) and shape (`sigma`).  Supporting functions include: `dRUBXII`, `pRUBXII`, `qRUBXII`, `rRUBXII`.

📌 See the R code: [`RUBXII_reg.R`](RUBXII_reg.R)  
📄 See detailed documentation: [docs/RUBXII_functions.md](docs/RUBXII_functions.md)

> **Reference:**  
> Ribeiro, T.F.; Cordeiro, G.M.; Peña-Ramírez, F.A.; Guerra, R.R. A New Quantile Regression for the COVID-19 Mortality Rates in the United States. Comput. Appl. Math. 2021, 40, 255.
---

### Unit Gamma Distribution

The `UG()` function defines a unit gamma distribution , parameterized by mean (`mu`) and precision (`sigma`).  Supporting functions include: `dUG`, `pUG`, `qUG`, `rUG`.

📌 See the R code: [`UG_reg.R`](UG_reg.R)  
📄 See detailed documentation: [docs/UG_functions.md](docs/UG_functions.md)

> **Reference:**  
> Mousa, A. M., El-Sheikh, A. A., & Abdel-Fattah, M. A. (2016). A gamma regression for bounded continuous variables. Advances and Applications in Statistics, 49, 305–326.

---

### Unit Lindley Distribution

The `ULindley()` function defines a unit Lindley distribution, , parameterized by mean (`mu`).  Supporting functions include: `dULindley`, `pULindley`, `qULindley`, `rULindley`.

📌 See the R code: [`ULindley_reg.R`](ULindley_reg.R)  
📄 See detailed documentation: [docs/ULindley_functions.md](docs/ULindley_functions.md)

> **Reference:**  
> Mazucheli, J., Menezes, A. F. B., & Chakraborty, S. (2019). On the one parameter unit-Lindley distribution and its associated regression model for proportion data. Journal of Applied Statistics, 46, 700–714.


---

### Unit Weibull Distribution

The `UW()` function provides a unit Weibull distribution , parameterized by median (`mu`) and shape (`sigma`).  Supporting functions include: `dUW`, `pUW`, `qUW`, `rUW`.

📌 See the R code: [`UW_reg.R`](UW_reg.R)  
📄 See detailed documentation: [docs/UW_functions.md](docs/UW_functions.md)

> **Reference:**  
> Mazucheli, J., Menezes, A. F. B., Fernandes, L. B., & de Oliveira, R. P. (2020). A quantile-based unit Weibull regression model for continuous bounded data. Statistical Papers, 61, 1549–1572.

---
