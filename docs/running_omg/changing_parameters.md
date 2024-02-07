---
layout: default
title: Changing Parameters
parent: Running OMG
nav_order: 2
---

# Changing Parameters
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

The default parameter set can overriden in two ways which are not mutually exclusive: 

1) passing additional parameters directly in the OMG function

2) via a pre-defined script.

---

## Changing parameters in the OMG function 

Parameters can be changed from their default value when calling the OMG function. The advantage of this is that parameters can be altered interactively, defined within loops, and you can explicitly see parameter values when running the model. Parameters in the OMG function call as a pair of additional inputs:

1) the parameter name as a string
   
2) the parameter value, which may be a string, logical, float etc...

For example, you can add a parameter name/value pair to define the name of the output directory to _basic_run_

```matlab
OMG(3000,'gen_pars.save_output_directory','basic_run'); 
```

You can change any number of parameters this way: 

```matlab
OMG(3000...
    	,'gen_pars.save_output_directory','basic_run'... 	% name of experiment
	,'bgc_pars.uptake_scheme','restoring'... 		% restore to PO4 observations
	,'bgc_pars.restore_timescale',30... 			% PO4 restoring timescale (days)
	,'bgc_pars.restore_data_file','PO4_Obs.mat'...		% PO4 observations
	);
```
{: .note }
To make things easier to read, comment and edit, use the line-continuation (three dots: `...`) to separate out each parameter pair onto different lines. 

{: .note }
Each parameter name/value pair begins with a comma and ends with `...` This makes it easy to quickly delete a line or copy/paste a new line in.

{: .note }
Adding a comment on each parameter line helps keep track of what you're doing.

---

## Changing parameters via a script

You can alternatively define parameters in a script that is read in when OMG runs. This is useful if you want to set a number of parameter values that won't be changed for a set of experiements, saving you to repeat the same lines of code in the OMG function call. Parameter scripts are kept in the _omg/config_files_ directory. The script contains the same parameter name/value pairs but are assignd like normal MATLAB variables. The script name is then used as the value for the special _ocean\_config_ parameter:

For example, you can include all the parameters from the example above in a script called _new_setup.m_:

```matlab
gen_pars.save_output_directory			='basic_run';  	% name of experiment
bgc_pars.uptake_scheme				='restoring';  	% restore to PO4 observations
bgc_pars.restore_timescale			=30;  	       	% restoring timescale (days)
bgc_pars.restore_data_file			='PO4_Obs.mat'; % observations
```

You then point OMG to _new_setup.m_ in the function call:

```matlab
% Run OMG with default parameters with a number of different parameter values
OMG(3000...
	,'ocean_config','new_setup'...
	);
```

---

## Changing parameters via both a script and the function call

Both options for setting parameters can be used simultaneously. For example, you can set a large number of parameters that are fixed for an ensemble of experiments and then change specific parameters for each experiment within in the function call:

```matlab
% Run OMG with default parameters with a number of different parameter values
OMG(3000...
	,'ocean_config','new_setup'...			% parameter script
	,'bgc_pars.restore_timescale',60 ... 		% a weaker restoring timescale for biological PO4 uptake (days)
	);
```

---

## Relevant Parameters
1. `ocean_config` (string) - the name of a script containing parameter values in `OMG/config_files` 


