---
layout: default
title: Forcings
parent: Running OMG
nav_order: 5
---

# Forcings
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

The ocean and atmosphere tracers in OMG can be perturbed with external forcings, _e.g.,_ a source of carbon to the atmosphere or nutrients to the ocean, or restored to a specified value _e.g.,_ atmospheric CO$_{2}$ is restored to future trajectory as defined in IPCC RCP scenarios. 

---

## The forcing directory

Forcings are defined in directories within the `OMG/forcings` directory, _e.g.,_ `OMG/forcings/my_forcing`. A single directory can contain information for multiple tracers. Each tracer to be forced has a set of files with the following extensions:

1. a `.dat` file that contains information about the type of forcing
2. a `.sig` file that contains information about how the forcing varies with time
3. (optional) a `.sur` file that contains information for spatially-varying forcings.

Each file should be named with the name of the tracer to be forced. The tracer is identifed by its name in OMG (see Section REF and the `setup_array_indices` in `OMG/code/general/general_functions.m`). For example, an atmospheric CO$_2$ forcing will have files: `pCO2.dat` and `pCO2.sig`. Forcings of multiple tracers can be defined within the same directory, each with their corresponding set of .dat, .sig and .sur files.

### The .dat file

This script contains options for the forcing: whether it is a perturbation or restoring forcing, and whether the forcing is applied uniformly or with a spatial pattern. A `1` in the corresponding column selects the option. Note: you cannot select restoring and forcing at the same time. 

```matlab
% col 1 - restoring?
% col 2 - force?
% col 3 - global forcing?
% col 4 - explicit 2D field? 
0 1 1 0 
```

### The .sig file

Ths file contains information about how the forcing varies in time. It contains the time-points in the lefthand column and corresponding scaling values in the righthand column. OMG will linearly interpolate between these points during the run. 

```matlab
0.0  1.0
9.0  1.0
10.0  0.0
999999.0  0.0
```

OMG will automatically add values for the very start of the run with a scaling of 0.0 if you don't specify this. 

### The .sur file

This file contains information about how the forcing varies spatially.

---

## Applying a forcing

The forcings are designed to allow flexibility in how they are applied. Each tracer has a corresponding parameter to specify a value at runtime, _e.g.,_ `bgc_pars.force_pCO2_val` that is combined with the information in the specified forcings directory. The forcing is calculated as:

$$
\text{Forcing} = \text{Spatial} * \text{Time} * \text{Forcing Value}
$$

where Forcing Value is the user-defined value, _e.g.,_ `bgc_pars.force_pCO2_val`, Time is the time-dependent scaling in the .sig file and Spatial is either defined explicity by the user in the .sur file or by choosing the global option in the .dat file.  Each component can have a scaling value, _e.g.,_ 0, 1, or a unit conversion value. At least one component must have a numerically meaningful value, _e.g.,_ 10 $\mu$mol year$^{-1}$. 

{:.note} 
Only the ocean can have a spatially-varying forcing.

---

## Example: Restoring atmospheric CO2

Create a directory containing the forcing information - for exampled called `restoring_atm_CO2`.

Within the directory, create a file called `pCO2.dat`, select options for a restoring of atmospheric CO$_2$ and global forcing:
```matlab
% col 1 - restoring?
% col 2 - force?
% col 3 - global forcing?
% col 4 - explicit 2D field? 
1 0 1 0 
```

Create a file called `pCO2.sig`, set the beginning and end timepoints with a scaling of 1.0 to indicate a constant time signal. The end timepoint can extend past the end of the experiment - here a large value ensures atmospheric CO$_{2}$ is always restored whatever the length of the run:
```maltab
0.0  1.0
999999.0  1.0
```

Run OMG specifying the forcing directory and restoring value:

```matlab
[output]=OMG(3000...
	,'bgc_pars.forcings_directory','restoring_atm_CO2'...
	,'bgc_pars.restore_pCO2',278e-6 ...
	);
```

---

## Example: Historical and future atmospheric CO2 using IPCC RCP scenarios

As in the above example, create a directory and the `.dat` file to select a global restoring of atmospheric CO$_{2}$. This example differs from that above in that the `.sig` file contains the actual values whilst the user-defined value `bgc_pars.restore_pCO2` keeps a default value of 1.0.
```matlab
1765.0  276e-6
1766.0  278e-6
.
.
.
```

Run OMG specifying the forcings directory. You can specify the restoring value but by default this is 1.0 so not neccessary. Note, the model is run with a starting year of 1765 corresponding with the forcings.
```matlab
[output]=OMG(3000...
	,'bgc_pars.forcings_directory','RCP_atm_CO2'...
	,'bgc_pars.restore_pCO2',1.0 ...
	);
```

---

## Example: Idealised atmospheric CO2 emission scenario

As above, create a directory and the same `.dat` file but now select a global forcing of atmospheric CO$_{2}$ in the `.dat` file:
```matlab
% col 1 - restoring?
% col 2 - force?
% col 3 - global forcing?
% col 4 - explicit 2D field? 
0 1 1 0 
```

An emission of CO$_2$ is applied over 1 year. The `.sig` file contains a scaling values. The scenario here is a top-hat function with a constant emission starting in year 1 of the experiment:
```matlab
0.0  0.0
1.0  1.0
2.0  1.0
2.0  0.0
9999.0  0.0
```

Run OMG specifying the forcings directory and the size of the emissions. Here we specify a release of 10 Pg C year$^{-1}$ and rescale the value to mol year$^{-1}$. 
```maltab
[output]=OMG(3000...
	,'bgc_pars.forcings_directory','Emissions_atm_CO2'...
	,'bgc_pars.force_pCO2',10*1e15/12.01 ...
	);
```

{:.note}
The rescaling can also be hidden in the `.sig` file by replacing `1.0` with `8.3e13` (a conversion factor of 1e15*12.01).

Alternatively, we can ramp up (and down) the emissions:
```matlab
0.5  1.0
1.0  1.0
1.5  1.0
2.0  0.0
9999.0  0.0
```

---

## Example: Adding a spatially-explicit flux of PO4 to the ocean surface

Similarly to the spatially-varying example, `.sig` file can contain 0s and 1s to simply scale a value given in `bgc_pars.force_*_val`. 
