---
layout: default
title: Airsea Gas Exchange
parent: Model Description
nav_order: 9
---


# Airsea Gas Exchange
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---
## Airsea Gas Exchange

OMG uses the method of \citet{Orr2017} to calculate air-sea gas exchange for a given tracer $\text{C}$.

$$
J_{gasex}^{C} = k_w \rho  (\text{C}_{sat} - \text{C}) (1-F_{seaice})
$$

where:

$$
k_w = a \big(\frac{\text{Sc}}{660} \big) ^{-1/2} u^{2}
$$

{:.note} 
Values for $a$ and $Sc$ are taken from \citet{Orr2017} (eqn. 13 and Table 1)

---

## Terms, Units and Relevant Parameters
* $J_{gasex}^{C}$ - flux of gas $\text{C}$ (mol m$^{-2}$ year$^{-1}$)
* $\rho$ - density (kg m$^{-3}$)
* $\text{C}_{sat}$ - saturation concentration of gas $\text{C}$ in equilibrium with atmosphere (mol kg$^{-1}$)
* $\text{C}$ - concentration of gas $\text{C}$ mol kg$^{-1}$)
* $F_{seaice}$ - fraction of seaice-cover (unitless) 
* $k_w$ - instanteous gas transfer velocity (m year$^{-1}$)
* $a$ - empirical constant to scale windspeed (?): `ocn_pars.gastransfer_a`
* $u$ - windspeed (m s$^{-1}$) ???
* $Sc$ - Schmidt number (unitless)
