---
layout: default
title: Forcings
parent: Model Description
nav_order: 11
---


# Forcings
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

All ocean and atmosphere state variables ($C$) can be subjected to an external forcing ($\Delta C$):

$$
J_{force}^{C} = C + \Delta C 
$$

State variables can also be restored to specified values ($C_{restore}$) over a specified timescale ($\tau_{restore}$):

$$
J_{force}^{C} = \frac{( C - C_{restore} )}{\tau_{restore}} 
$$

---

## Relevant Parameters
* $\Delta C$ - forcing of tracer (mol year$^{-1}$)
* $C_{restore}$ - restoring target (mol year$^{-1}$ in the ocean, mol mol$^{-1}$ in the atmosphere)
* $\tau_{restore}$ - restoring timescale (year$^{-1}$)
