---
layout: default
title: Equilibrium Solutions
parent: Running OMG
nav_order: 7
---

# Equilibrium Solutions
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---
## Equilibrium Solutions

OMG can (attempt to) find a steady-state or equilibrium solution directly. This approach can be significantly faster than forward timestepping towards steady-state or equilibrium (e.g., it takes $\sim$30 secs for the examples below) but is not guaranteed to find a solution. OMG will report an error saying so if this the case. If OMG converges on a solution, it may be not be meaningful, so it is recommended you check! Finally, equilibrium solutions also require a more nuanced set-up of the model. 

Currently OMG uses two algorithms:

1) Newton's method - only applicable for annual average models
2) Jacobian Free Newton Krylov methods - applicable for both annual and seasonally varying models

Once a solution has been found, OMG will time-step for the specified runtime in the function call and generate output as with forward time-stepping.

---
## Newton's method

OMG uses the solver by \citep[the \textit{nsold} algorithm from][]{Kelley2003}. This algorithm requires that the model set-up uses an annually-averaged transport matrix.  To use this option, set `'gen_pars.integrate_scheme'` parameter to `newton`. 

---
## Jacobian Free Newton Krylov methods

OMG uses the Jacobian Free Newton Krylov solver by \citep[the `nsoli' algorithm in][]{Kelley2003}. This option can be run with an annually-forced model and in some cases seems to find a solution in less time. 

To use this option, set `'gen_pars.integrate_scheme'` parameter to `newton_krylov`. 

---
### Additional Parameters Needed for all Equilibrium Solutions

Mass is not conserved when searching for an equilibrium solution so additional parameters need to set to ensure the solution is associated with the correct tracer inventories. This is achieved by selecting a forcing that restores the average concentration to a defined value at an appropriate timescale. 

The following example spins up the carbon cycle using the newton solver. PO$_4$ and alkalinity are restored globally to an average value on a geologically-relevant timescale \cite{LiPrimeau2008}. Note that DOP is controlled by the PO$_4$ inventory and biological fluxes. Atmospheric CO$_2$ is restored to a pre-industrial value and controls the inventory of dissolved inorganic carbon. 

```matlab
    outputdir = OMG(1...
        ,'ocean_config'                      , 'worbe2'...
        ,'gen_pars.integrate_scheme'         , 'newton'...
        ,'bgc_pars.CARBCHEM_select'          , true...
        ,'gen_pars.save_output_directory'    , 'newton_spinup'...
        ,'gen_pars.forcings_directory'       , 'restore_PO4_restore_ALK_restore_pCO2'...
        ,'bgc_pars.restore_PO4_val'          , 2.15e-6...
        ,'bgc_pars.restore_PO4_timescale'    , 53000...
        ,'bgc_pars.restore_ALK_val'          , 2363e-6...
        ,'bgc_pars.restore_ALK_timescale'    , 53000...
        ,'bgc_pars.restore_pCO2_val'         , 278e-6 ...
        ,'bgc_pars.restore_pCO2_timescale'   , 1...
    );
```

The newton algorithm could also be used with freely-evolving atmospheric CO$_2$, e.g., when varying a biological pump parameter, by removing the corresponding forcing and using the above example as the restart.

### Finer Details: Newton

The model constructs the Jacobian numerically via the native MATLAB \texttt{numjac} function using a sparsity pattern. To make this as independent of the tracers selected as possible, a sparsity pattern is first estimated using the transport matrix and particulate matrix sparsity patterns which is simply repeated for every tracer. The texttt{numjac} function is then used to remove the excess non-zeros.

The default parameters for the \texttt{nsold} functions are used which uses a Newton-Armijo algorithm where the Jacobian is updated after 1000 nonlinear iterations or whenever the ratio of successive norms is greater than 0.5. The absolute and relative tolerances are set to 1e-10. These parameters are found in the `solve_newton` function in `code/OMG/OMG_functions.m` if you want to alter them.

### Finer Details: Newton Krylov

The model applies a left preconditioner using the preconditioner definition outlined by \citet{Khatiwala2008}. This is calculated using the interpolated transport matrix at the mid-point of the model period (1 year) and the same `numjac` routine used above. The preconditioner is calculated in the first function evaluation and then stored. Whilst \citet{Khatiwala2008} suggests updating the preconditioner with each nonlinear iteration, in practice, we find an initial calculation is sufficient to find solutions. 

The default parameters for the `nsoli` functions are used which sets an upper limit on nonlinear iterations and linear iterations per nonlinerar iteration as well specifying the GMRES algorithm. The absolute and relative tolerances are set to 1e-10. These parameters are found in the `solve_newton_krylov` function in `code/OMG/OMG_functions.m` if you want to alter them.
