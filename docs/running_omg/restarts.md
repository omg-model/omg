---
layout: default
title: Restarts
parent: Running OMG
nav_order: 4
---

# Restarts
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

### Initial conditions

The default option in OMG is to initialise tracers with uniform distributions defined by global mean concentrations set in `OMG/code/biogeochemistry/initialise_bgc.m`. 

## Initialising OMG from a previous run

To alternatively initialise tracers from a previous OMG run, set `gen_pars.OMG_restart_file` to the name of the relevant output directory:

```matlab
[output]=OMG(3000...
	,'gen_pars.OMG_restart_file','spinup' ...
	);
```

---

## Relevant Parameters

* `gen_pars.netcdf_restart`: restart from GENIE netCDF files? (true/false)
* `gen_pars.OMG_restart_file`: name of experiment to restart from (string)
* `bgc_pars.PO4_init`: initial concentration of phosphate ($\mu$mol kg$^{-1}$)
* `bgc_pars.DOP_init`: initial concentration of dissolved organic phosphorus ($\mu$mol kg$^{-1}$)
* `bgc_pars.O2_init`: initial concentration of oxygen ($\mu$mol kg$^{-1}$)
* `bgc_pars.DIC_init`: initial concentration of dissolved inorganic carbon ($\mu$mol kg$^{-1}$)
* `bgc_pars.ALK_init`: initial concentration of alkalinity ($\mu$mol kg$^{-1}$)
* `bgc_pars.pCO2_init`: initial concentration ofatmospheric CO$_2$ ($\mu$mol kg$^{-1}$)
