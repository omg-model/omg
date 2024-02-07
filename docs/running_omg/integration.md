---
layout: default
title: Integration and Equilibrium Solutions
parent: Running OMG
nav_order: 6
---

# Integration and Equilibrium Solutions
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
So far this option does not seem as useful as forward integration or equilibrium solutions. 

OMG can be alternatively integrated using the suite of ODE solvers in MATLAB. To do this set `gen_pars.integrate_scheme` to the name of one of the solvers, e.g., `ode45`. Options for the ode solvers can be set via `gen_pars.odeoptions.<value>`, with default options set in `code/general/gen_params.m`.

{:.note}
Solvers aimed at stiff solutions work better but may need additional information like the Jacobian.

---

## Equilibrium solutions for annual-average circulations

OMG can (attempt to) find a steady-state solution for the state-variables directly using Newton's method \citep[the \textit{nsold} algorithm from][]{Kelley2003}. This algorithm requires that the model set-up uses an annually-averaged transport matrix. This option can be significantly faster to find a steady-state solution (e.g., it takes $\sim$30 secs for the example below) but the newton algorithm is not guaranteed to converge on a solution! In this case, OMG will report an error saying so. If the newton algorithm does converge on a solution, it may be not be meaningful, so it is recommended you check!

To use this option, set the integrate scheme option to `newton`. The run will also need a corresponding set of restoring forcings to ensure that the solution has the correct inventory as this is not conserved in this approach. The following example spins up the carbon cycle using the newton scheme. PO$_4$ and alkalinity are restored globally to an average value on a geologically-relevant timescale \cite{LiPrimeau2008}. Note that DOP is controlled by the PO$_4$ inventory and biological fluxes. Atmospheric CO$_2$ is restored to a pre-industrial value and controls the inventory of dissolved inorganic carbon. In the example the runtime is set to 1 year. OMG will find a steady-state solution and then time-step the model for runtime to save output (The runtime can be set to longer).

```matlab
    outputdir = OMG(1...
        ,'ocean_config'                      , 'worbe2_annual'...
        ,'gen_pars.integrate_scheme'         , 'newton'...
        ,'bgc_pars.CARBCHEM_select'          , true...
        ,'gen_pars.save_timeslice_output'    ,[99.5:100:899.5 999.5:1000:9999.5]...
        ,'gen_pars.save_timeseries_output'   ,[0.5:1:99.5 100.5:100:999.5 1000.5:1000:9999.5]...
        ,'gen_pars.save_output_directory'    , 'newton_spinup'...
        ,'gen_pars.dt_ratio'                 , 1...
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

### Finer Details

The model constructs the Jacobian numerically via the native MATLAB \texttt{numjac} function using a sparsity pattern. To make this as independent of the tracers selected as possible, a sparsity pattern is first estimated using the transport matrix and particulate matrix sparsity patterns which is simply repeated for every tracer. The texttt{numjac} function is then used to remove the excess non-zeros.

The default parameters for the \texttt{nsold} functions are used which uses a Newton-Armijo algorithm where the Jacobian is updated after 1000 nonlinear iterations or whenever the ratio of successive norms is greater than 0.5. The absolute and relative tolerances are set to 1e-10. These parameters are found in the `solve_newton` function in `code/OMG/OMG_functions.m` if you want to alter them.

---

## Equilibrium solutions for seasonally-varying circulations

OMG can also (attempt to) find a steady-state solution for the state-variables directly in a seasonally-forced model, i.e., a periodic solution. This uses a Jacobian Free Newton Krylov solver  \citep[the `nsoli' algorithm in][]{Kelley2003}. This option can be run with an annually-forced model and in some cases seems to find a solution in less time.

To use this option, set the integrate scheme option to `newton_krylov`. As with the Newton algorithm above, the same forcings are required and OMG will time-step for the specified runtime and generate output.

### Finer Details

The model applies a left preconditioner using the preconditioner definition outlined by \citet{Khatiwala2008}. This is calculated using the interpolated transport matrix at the mid-point of the model period (1 year) and the same `numjac` routine used above. The preconditioner is calculated in the first function evaluation and then stored. Whilst \citet{Khatiwala2008} suggests updating the preconditioner with each nonlinear iteration, in practice, we find an initial calculation is sufficient to find solutions. 

The default parameters for the `nsoli` functions are used which sets an upper limit on nonlinear iterations and linear iterations per nonlinerar iteration as well specifying the GMRES algorithm. The absolute and relative tolerances are set to 1e-10. These parameters are found in the `solve_newton_krylov` function in `code/OMG/OMG_functions.m` if you want to alter them.
