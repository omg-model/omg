---
layout: default
title: Calcium Carbonate
parent: Model Description
nav_order: 7
---


# Calcium Carbonate
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Export

The production of calcium carbonate (CaCO$_3$) in the surface ocean is modelled as:

$$
J_{up}^{CaCO_3} =\underbrace{\big( J_{up}^{POP} \cdot R_{C:P} \big)}_{\text{POC export}}\cdot \underbrace{r_{0}^{CaCO_3:POC}}_{\text{PIC:POC ratio}} \cdot \underbrace{\gamma}_{\Omega \text{ modifier}}
$$

where $r_{0}^{CaCO_3:POC}$ is a thermodynamically-based modifier based on the saturation state of CaCO$_3$:

$$
\gamma=(\Omega-1)^\eta,		\Omega>1.0 \\
\gamma=0.0, \Omega\leqslant1.0
$$

{:.note}
To set a globally uniform PIC:POC ratio, set `bgc_pars.red_PIC_POC_mod = 0.0`

---

## Remineralisation

The remineralisation of CaCO$_3$ in the water column in modelled as:

$$
J_{diss}^{CaCO_3} = F^{CaCO_3}_{z}=\underbrace{(1-r) J_{up}^{CaCO3} \cdot exp(\dfrac{(z-z0)}{el_1})}_\text{fraction 1} + \underbrace{rJ_{up}^{CaCO3} \cdot exp(\dfrac{(z-z0)}{el_2})}_\text{fraction 2}
$$

{:.note}
All excess CaCO$_{3}$ not remineralised is remineralised entirely in the deepest grid-box of the corresponding water column.

---

## Terms, Units and Relevant Parameters
* $J_{up}^{CaCO3}$ - production of CaCO$_3$ in surface ocean (mol kg$^{-1}$ year$^{-1}$)
* $J_{up}^{POP}$ - production of organic matter in surface ocean (mol kg$^{-1}$ year$^{-1}$)
* $R_{C:P}$ - stoichiometic ratio of carbon to phosphorus (unitless): `bgc_pars.red_P_to_C`
* $r_{0}^{CaCO_3:POC}$ - PIC:POC ratio (unitless): `bgc_pars.red_PIC_POC`
* $\gamma$ - local saturation state modifier (unitless) 
* $\Omega$ - saturation state of CaCO$_3$ (unitless) 
* $\eta$ - exponent controlling sensitivity to $\Omega$ (unitless): `bgc_pars.red_PIC_POC_mod`
* $J_{diss}^{CaCO3}$ - remineralisation of CaCO$_3$  (mol kg$^{-1}$ year$^{-1}$)
* $z0$ - bottom depth of surface layer (m)
* $z$ - depth (m)
* $r$ - fraction of CaCO$_3$ production to fraction 2 (unitless): `bgc_pars.CaCO3_frac_2` 
* $eL_1$ - \textit{e}-folding depth of first fraction (m): `bgc_pars.CaCO3_eL1`
* $eL_2$ - \textit{e}-folding depth of second fraction (m): `bgc_pars.CaCO3_eL2`


