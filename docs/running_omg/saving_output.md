---
layout: default
title: Saving Output
parent: Running OMG
nav_order: 3
---

# Saving Output
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

OMG can save output in multiple ways that are suitable for different tasks. Output can be saved as:

1) text format (.res) for timeseries of averaged values such as inventories or global means.
2) netCDF (.nc) for 2D/3D fields, which are more suitable for sharing and use by anyone less familiar with MATLAB.
3) matlab variables in the workspace, for immediate use

All output will be written to the _OMG/experiments/_ directory. You can specify a name for the experiment output directory. If not specified, a default name containing the date and time will be generated. 

---

## Specifying when OMG saves output

OMG will automatically save the last year of the run to avoid you accidentally missing it! Output years for T=text format (timeseries) and netcdf (timeslices) are specified seperarely as netcdf saving can slow a run down.

You can select output years by specifying an array using the `[ ]` with the years you want saved, e.g., specifying 2.0 will save the 2nd year of the simulation. You can specify a custom choice of years or use the `:` operator or functions like `linspace` to easily specify an orderded choice.  This array is then passed as the input for either the timeseries or timeslice output parameters:

```matlab
% setup arrays with years to save output from
netcdf_output=[4999.0];
timeseries_output=[0:1:99 100:200:999 1000:1000:49999];

% run experiment 
OMG(3000...
	,'gen_pars.save_output_directory','basic_run'...         % name of output directory
	,'gen_pars.save_timeslice_output',netcdf_output...       % timeslice years
	,'gen_pars.save_timeseries_output',timeseries_output...  % timeseries years
	);
```

---

## Intra-annual outputs

By default. OMG will save an annually averaged output but each save-point can be further split into smaller averaging periods such as seasonal or monthly averaged periods.

```matlab
% setup arrays with years to save output from
netcdf_output=[4999.0];
timeseries_output=[0:1:99 100:200:999 1000:1000:49999];

% run experiment 
OMG(3000...
	,'gen_pars.save_output_directory','basic_run'...         % name of output directory
	,'gen_pars.save_timeslice_output',netcdf_output...       % timeslice years
	,'gen_pars.save_timeseries_output',timeseries_output...  % timeseries years
	.`gen_pars.save_intra',12...				 % save monthly outputs for each year selected
	);
```

---

## Saving output directly to the MATLAB workspace

Lastly, output can be loaded directly into the MATLAB workspace like a normal function. This will create a structure array in your workspace. This option allows the collation of specific ouputs within an ensemble 

```matlab
[output]=OMG(3000,...
	'gen_pars.save_output_directory','basic_run',... % name of experiment
	'gen_pars.save_timeslice_output',netcdf_output,... % timeslice years
	'gen_pars.save_timeseries_output',timeseries_output,... % timeseries years
	);
```

---

## Relevant Parameters

1. `gen_pars.save_output_directory` (string): name of output directory to write output to
2. `gen_pars.save_intra` (integer): number of intra-annual steps to save (1 = annual, 4 = seasonal, 12 = monthly)
3. `gen_pars.save_timeseries_output` (numerical array): year(s) for text output, e.g., [0:1:10])
4. `gen_pars.save_timeslice_output` (numerical array): year(s) for 2D/3D outputs,  e.g., [10]
5. `gen_pars.save_matlab_output` (logical): save as MATLAB format?

