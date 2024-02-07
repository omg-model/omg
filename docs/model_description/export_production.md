---
layout: default
title: Biological Uptake
parent: Model Description
nav_order: 4
---


# Biological Uptake
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Uptake of PO4

The uptake of PO$_4$ and net export production of organic matter is selected by `bgc_pars.uptake_scheme`.

### Michaelis-Menten / Monod

`bgc_pars.uptake_scheme = 'MM'`


$$ J_{up}=V_{max} \cdot \dfrac{\text{PO}_4}{K_{PO4}+\text{PO}_4} \cdot \dfrac{I}{I_0}  (1-F_{seaice}) $$

### Restoring to Observations

`bgc_pars.uptake_scheme = 'restore'`

$$ J_{up}= \frac{1}{\tau}   \text{max}\big( (\text{PO}_4 - \text{PO}_4^{obs}),0 \big)   (1-F_{seaice}) $$

### Fixed

`bgc_pars.uptake_scheme = ??`

$$
 J_{up}  = \begin{cases}
 	\phi_{prescribed} , &  \text{if}  ([\text{PO}_4] - \phi_{prescribed} )> 0 \\
 	0, & \text{otherwise}
 	\end{cases}
$$

### No Biology

`bgc_pars.uptake_scheme = 'strangelove'`

$$ J_{up}= 0.0 $$

---

## Terms, Units and Relevant Parameters

1. $J_{up}$ - net biological uptake of PO$_4$ (mol kg$^{-1}$ year$^{-1}$)
2. $\text{PO}_{4}$ - PO$_4$ concentration in grid box (mol kg$^{-1}$)
3. $V_{max}$ - maximum uptake rate (mol kg$^{-1}$ year$^{-1}$): `bgc_pars.u0PO4`
4. $K_{PO4}$ - half saturation constant (mol kg$^{-1}$): `bgc_pars.KPO4`
5. $I$ - insolation (units??)
6. $I_0$ - reference insolation
7. $F_{seaice}$ - fraction of seaice cover (unitless)
8. $\tau$ - restoring timescale (days): `bgc_pars.restore_timescale`
9. $\text{PO}_4^{obs}$ - PO$_4$ observations (mol kg$^{-1}$): `bgc_pars.restore_data_file`
10. $F_{seaice}$ - fraction of seaice cover (unitless)

