---
layout: default
title: A Basic Run
parent: Running OMG
nav_order: 1
---

OMG is designed to be run simply by calling as a function in MATLAB. The function has one mandatory input: the experiment length in years. OMG will run with a 16-level seasonally-averaged preindustrial set-up corresponding to Cao et al., (2009)

```matlab
OMG(3000)
```

The output will be written to OMG/output/ with an automatic experiment name consisting of the date and a unique number for that date. 
