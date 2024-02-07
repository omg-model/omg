---
layout: default
title: A Basic Run
parent: Running OMG
nav_order: 1
---

OMG is designed to be run simply by calling as a function in MATLAB. 

The function has one mandatory input: the experiment length in years. OMG will run with a 8-depth-level annual-averaged preindustrial set-up corresponding to [Ridgwell et al., (2007)](https://bg.copernicus.org/articles/4/87/2007/bg-4-87-2007.html). The model will be time-stepped forward in time with output from the final year saved.

```matlab
OMG(3000);
```

The output will be written to OMG/output/ with an automatically generated experiment name consisting of the date and a unique number for that date. 
