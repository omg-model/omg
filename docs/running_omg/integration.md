---
layout: default
title: Integration
parent: Running OMG
nav_order: 6
---

# Integration
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Forward Integration

OMG can be integrated using the forward Euler method ("forward timestepping"). To do this, set `gen_pars.integrate_scheme` to `fwd`. This is the default option for OMG.

---

## Integration using MATLAB's ODE solvers

{:.warning}
This option works but is not recommended. Forward timestepping or equilibrium solutions are faster! 

OMG can be alternatively integrated using the suite of ODE solvers in MATLAB. To do this set `gen_pars.integrate_scheme` to the name of one of the solvers, e.g., `ode45`. Options for the ode solvers can be set via `gen_pars.odeoptions.<value>`, with default options set in `code/general/gen_params.m`.

{:.note}
Solvers aimed at stiff solutions work better but may need additional information like the Jacobian.

