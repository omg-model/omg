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

OMG can save output in multiple ways that are suitable for different tasks. Output can be saved as text format for timeseries and netCDF for 2D/3D fields, which are more suitable for sharing and use by anyone less familiar with MATLAB.  Different save points for timeseries and timeslices can be saved to avoid excess storage/computational time. Output can also be saved in native MATLAB format files for analysis in MATLAB. 

---

## Specifying when OMG saves output

OMG will automatically save the last year of the run to avoid you accidentally missing it! You can select other years by specifying an array using the `[ ]` with the mid-point of years you want saved, e.g., specifying 1.5 will save the 2nd year of the simulation. You can specify a custom choice of years or use the `:` operator or functions like `linspace` to easily specify an orderded choice.  

By default. OMG will save an annually averaged output but each save-point can be further split into smaller averaging periods such as seasonal or monthly averaged periods.

All output will be written to the _OMG/experiments/_ directory. You can specify a name for the experiment. If not specified, a default name containing the date and time will be generated. 

```matlab
% setup arrays with save-points
netcdf_output=[9999.5];
timeseries_output=[0.5:10:99.5 100.5:200:999.5 1000.5:1000:49999.5];

% run experiment
OMG(3000...
	,'gen_pars.save_output_directory','basic_run'...         % name of output directory
	,'gen_pars.save_timeslice_output',netcdf_output...       % timeslice years
	,'gen_pars.save_timeseries_output',timeseries_output...  % timeseries years
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
3. `gen_pars.save_timeseries_output` (numerical array): year mid-points for text output, e.g., [0.5:1:9.5])
4. `gen_pars.save_timeslice_output` (numerical array): year mid-points for 2D/3D outputs (year mid-points, e.g., [9.5])
5. `gen_pars.save_matlab_output` (logical): sae as MATLAB format?

