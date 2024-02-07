---
layout: default
title: Governing Equations
parent: Model Description
nav_order: 2
---

# Governing Equations
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

The default state variables are phosphate (PO4) and dissolved organic phosphorus (DOP). In the following equaitons, $\mathbf{A}$ is the transport matrix that represents the net ocean transport diagnosed from GENiE (see Section LINK). Biogeochemical source/sink terms ($J$) are described in the following sections.

---

## Phosphate Cycle

$$ \frac{d\text{PO}_4}{dt} = \mathbf{A} \text{PO}_{4} - \underbrace{\big( J_{up}^{POP} + J_{up}^{DOP} \big)}_{\text{Net export production -} J_{up}}+ J_{remin}^{POP} + J_{remin}^{DOP} + J_{force}^{PO4} $$

$$ \frac{d\text{DOP}}{dt} = \mathbf{A} \text{DOP} + J_{up}^{DOP} - J_{remin}^{DOP} + J_{force}^{DOP} $$

---

## Carbon Cycle

A carbon cycle can be selected by setting `bgc_pars.CARBCHEM_select=true`:

$$ \frac{d\text{DIC}}{dt} = \mathbf{A} \text{DIC} - \big( J_{up}^{POP} + J_{up}^{DOP} \big) R_{C:P} - J_{up}^{CaCO_3} +  J_{remin}^{POC} + J_{remin}^{DOC} + J_{diss}^{CaCO_3} + J_{gasex}^{\text{CO}_2} + J_{force}^{DIC} $$

$$ \frac{d\text{ALK}}{dt} = \mathbf{A} \text{ALK} - \big( J_{up}^{POP} + J_{up}^{DOP} \big) (-R_{N:P})  + \big( J_{remin}^{POP} + J_{remin}^{DOP} \big) R_{N:P} -2\big( J_{up}^{CaCO_3} + J_{diss}^{CaCO_3} \big)+ J_{force}^{ALK} $$

$$ \frac{dx\text{CO}_2}{dt} = x\text{CO}_2 - J_{gasex}^{\text{CO}_2} + J_{force}^{xCO_{2}} $$

---

## Oxygen Cycle

An oxygen cycle can be selected by setting `bgc_pars.O2_select=true`:

$$ \frac{d\text{O}_2}{dt} = \mathbf{A} \text{O}_2 - \big( J_{up}^{POP} + J_{up}^{DOP}  \big) R_{O:P} + \big( J_{remin}^{POP} + J_{remin}^{DOP} \big) R_{O:P} + J_{gasex}^{\text{O}_2} + J_{force}^{O_{2}} $$

$$ \frac{dx\text{O}_2}{dt} = x\text{O}_2 - J_{gasex}^{\text{O}_2} + J_{force}^{xO_{2}} $$

---

## Revelant Parameters

1. `bgc_pars.CARBCHEM_select` (logical): turn carbon cycle tracers on (DIC, Alkalinity and atmospheric CO2)
2. `bgc_pars.O2_select` (logical): turn oxyegn cycle tracers on (dissolved oxygen, atmospheric O2)


