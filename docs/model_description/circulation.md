---
layout: default
title: Ocean Circulation
parent: Model Description
nav_order: 3
---


# Ocean Circulation
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Ocean Circulation via the transport matrix

The net transport of dissolved tracers by ocean circulation from timestep $n$ to $n+1$ is represented as a transport matrix \citep{Khatiwala2005,Khatiwala2007} diagnosed from cGENIE:

$$ \mathbf{c^{n+1}} = \mathbb{A} \mathbf{c^{n}} + \mathbf{q} $$

where: 

* $\mathbf{c}$ - a tracer (mol kg$^{-1}$)
* $\mathbf{A}$ - the transport matrix
* $\mathbf{q}$ - sources/sinks (mol kg$^{-1}$ year$^{-1}$)

Each transport matrix is associated with boundary conditions:

* $F_{seaice}$ - the fractional cover of seaice (unitless, from 0 to 1)
* $T$ - temperature (deg C)
* $S$ - salinity (PSU)
* $u$ - windspeed (m s$^{-1}$)
* $I$ - insolation (units?!)

The transport and associated boundary conditions are available as annually or monthly averaged. The matrix and boudnary conditions are linearly interpolated between monthly averages when timestepping. The transport matrix is discrete in time and so is associated with the timestep that the original GENIE experiment was run at. 

---

## Construction of the GENIE Transport Matrix

The Transport Matrix as described in \citet{Khatiwala2005} and \citet{Khatiwala2007} contains the finite difference tendency calculated in an ocean model:

$$ \frac{d\textbf{c}}{dt}=\frac{\textbf{c}^{n+1}-\textbf{c}^{n}}{\Delta t}=\textbf{A}'^{t}\textbf{c}^{n}+\textbf{q}'^{n} $$

where $\textbf{A}'$ is the transport matrix containing the finite difference tendencies, $\textbf{c}$ is a tracer (mol kg$^{-1}$) and $\textbf{q}'$ is a vector of source/sinks (mol kg$^{-1}$ t$^{-1}$).  The superscripts ($n$) refer to the time step index.  Rearranging Equation (\ref{eq:TM_tendency}) for $\textbf{c}^{t+1}$:

$$ \textbf{c}^{n+1}=(\textbf{I}+\textbf{A}'\Delta t)\textbf{c}^{n}+\textbf{q}'^{n}\Delta t $$

The equation above predicts the tracer field after one timestep ($\textbf{c}^{n+1}$) from the effect of ocean circulation ($\textbf{I}+\textbf{A}'\Delta t$) on the tracer at the previous timestep ($\textbf{c}^{n}$) plus any sources or sinks for the tracer over the timestep ($\textbf{q}'^{n}\Delta t$). In GENIE each tracer experiment diagnoses the tracer distribution resulting from the modelled ocean circulation acting on a unit concentration in a single grid-box, i.e., 1 mol kg$^{-1}$, after one time step.  The resulting concentrations form the coefficients in each column. Therefore, the GENIE transport matrix ($\textbf{A}$) described and used in this study is equivalent to ($\textbf{I}+\textbf{A}'\Delta t$) from equation (\ref{eq:TM_steadystate}).

---

