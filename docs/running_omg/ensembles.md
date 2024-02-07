---
layout: default
title: Running Ensembles
parent: Running OMG
nav_order: 7
---

# Running Ensembles 
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

OMG can be used to run an ensemble of experiments easily. All that needs to be done is to create an array of parameter values, place the OMG function call within a loop and index the array to get the parameter value. 

---

## An example ensemble

In this example, we setup an array of _e_-folding depths for POC remineralisation. We first run a spinup as a starting point for each experiment. We then create a loop within which we create a unique experiment name, call OMG and input the _e_-folding depth from the array of values:

```matlab
% define perturbed parameter array
efold_depth=linspace(50,3000,10);

% run a spinup
OMG(10000...
	,'gen_pars.save_output_directory','spinup'...
	);

% run the ensemble
for n=1:numel(efold_depth)

	filename=strcat('efold_exp_',num2str(n)); % unique experiment name
	
	OMG(3000...
    	,'gen_pars.save_output_directory',filename ... % exp. name
      ,'gen_pars.OMG_restart_select',true... % restart
	    ,'gen_pars.OMG_restart_file','spinup'... % restart experiment name
		  ,'bgc_pars.POC_eL1',efold_depth(n,1)... % parameter value from array
		);
end	
```

---

##  A parallel ensemble varying 2 parameters

{:.warning}
This requires the parallel toolbox

You can create an ensemble varying multiple parameters simultaneously. In this case the number of experiments required to explore the model behaviour across the parameter input space can get big very quickly. If you have a computer with multiple cores and the parallel toolbox installed, you can easily save time by invoking a parallel for loop by changing `for` to `parfor`.

The following code shows an example ensemble varying the biological uptake rate of PO4 and the remineralisation length scale of POM and quantifying the impact on atmospheric CO2. To view the key output quickly the CO2 output is saved directly to the workspace and copied to an array which can be plotted against the inputs. 

```matlab
% define perturbed parameter arrays
POC_efold_ref=550.5195;
POC_uptake_ref=1.9582242E-06;
POC_efold=linspace(POC_efold_ref/2,POC_efold_ref*2,10);
POC_uptake=linspace(POC_uptake_ref/2,POC_uptake_ref*2,10);

% combine and reshape to get parameter list
n_params=numel(POC_efold)*numel(POC_uptake);
params=zeros(n_params,2); 
[X,Y] = meshgrid(POC_efold,POC_uptake);
params(:,1)=reshape(X,n_params,1);
params(:,2)=reshape(Y,n_params,1);

% initialise an output array for CO2
CO2=zeros(n_params,1); 

% run parallel ensemble
parfor n=1:numel(efold_depth)
	[output]=OMG(3000,...
		'bgc_pars.POC_eL1',params(n,1),...
		'bgc_pars.u0PO4',params(n,2),...
		);
	
	CO2(n,1)=output.CO2;	
end	

% plot output
contourf(X,Y,CO2)
```
