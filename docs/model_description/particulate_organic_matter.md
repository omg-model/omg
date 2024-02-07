---
layout: default
title: Particulate Organic Matter
parent: Model Description
nav_order: 6
---

# Particulate Organic Matter
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Production

The production of POP is scaled to the total export production:

$$ J_{up}^{POP} = J_{up} (1 - \lambda_{DOP}) $$

---

## Remineralisation

### Implicit Double Exponential 

`bgc_pars.remin_scheme = 'exponential'`

$$ J_{remin}^{POP} (z)  =\underbrace{(1-r) J_{up}^{POP} \cdot exp(\dfrac{(z-z0)}{el_1})}_\text{fraction 1} + \underbrace{rJ_{up}^{POP} \cdot exp(\dfrac{(z-z0)}{el_2})}_\text{fraction 2} $$

### Implicit Power Law (Martin Curve) 

`bgc_pars.remin_scheme = 'martin'`

$$ J_{remin}^{POP} (z) = J_{up}^{POP} \cdot \big(\frac{z}{z0}\big)^{-b} $$

{: .note } 
Any POP remaining at the seafloor is remineralised in the deepest grid-box of the water column.

{: .note } 
Implicait reminerailsation schemes can be made spatially variable by passing a 1D array to `bgc_pars.POC_martin_b` 

---

## Terms, Units and Relevant Parameters

* $J_{up}$ - net export production (mol kg$^{-1}$ year$^{-1}$)
* $J_{up}^{POP}$ - export of POP (mol kg$^{-1}$ year$^{-1}$)
* $J_{remin}^{POP}$ - remineralisation of POP (mol kg$^{-1}$ year$^{-1}$)
* $\lambda_{DOP}$ - fraction of export production to DOP (unitless): `bgc_pars.DOP_frac`
* $z$ - depth (m)
* $r$ - fraction of POM to second exponential decay term (unitless): `bgc_pars.POC_frac_2` 
* $el_1$ - _e_-folding depth of POM in first exponential decay term (m): `bgc_pars.POC_eL1` 
* $el_2$ - _e_-folding depth of POM in second exponential decay term (m): `bgc_pars.POC_eL2` 
* $b$ - exponent controlling attenutation (unitless): `bgc_pars.POC_martin_b` 


