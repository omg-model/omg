---
layout: default
title: Dissolved Organic Matter
parent: Model Description
nav_order: 5
---

# Dissolved Organic Matter
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Production

The production of DOM is scaled to the total organic matter production: 

$$ J_{up}^{DOP} = J_{up} \lambda_{DOP} $$

---

## Remineralisation

The remineralisation of DOM is calculated as:

$$ J_{remin}^{DOP} = \tau_{DOP}\text{DOP} $$

---

$$ Terms, Units and Relevant Parameters

* $J_{up}$ - net export production (mol kg$^{-1}$ year$^{-1}$)
* $J_{up}^{DOP}$ - export production of DOP (mol kg$^{-1}$ year$^{-1}$)
* $\lambda_{DOM}$ - fraction of export production to DOM (unitless): `bgc_pars.DOP_frac`
* $J_{remin}^{DOP}$ - remineralisation of DOP (mol kg$^{-1}$ year$^{-1}$)
* $\text{DOP}$ - concentration of DOP (mol kg$^{-1}$)
* $\tau_{DOM}$ - lifetime of DOM (years): `bgc_pars.DOP_k` 

