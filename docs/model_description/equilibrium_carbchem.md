---
layout: default
title: Equilibrium Carbonate Chemistry
parent: Model Description
nav_order: 8
---


# Equilibrium Carbonate Chemistry
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

OMG uses the method of \citet{Follows2006} to solve the local equilibrium carbonate chemistry from Dissolved Inorganic Carbon and Alkalinity. This is primarily used to find the partial pressure of CO2 in the surface ocean in order to calculate $J_{gasex}^{CO_2}$ (Section REF) and the saturation state of calcite ($\Omega$) to calculate $J_{up}^{CaCO_3}$ (Section REF)

{:.note}
The solver will take 20 steps before giving a fatal warning that will stop the simulation
