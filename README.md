<<<<<<< HEAD
# OMG

Offline Matrix GENIE (OMG) is MATLAB code to run a basic Transport Matrix version of GENIE based on work developed in my PhD thesis. 

Created by Jamie D Wilson (jamie.wilson@bristol.ac.uk) The code is **unpublished** and very much **unfinished!**

To run the model: 

  (1) Requires GENIE transport matrix and associated files (contact me: jamie.wilson@bristol.ac.uk)
  
  (2) A config file in OMG/config_files contains default parameters specific to the run in GENIE you want to use.
  
  (3) The function OMG() is called with two parameters; runtime in years and the name of the config file. Default parameters can be overridden by passing them as a pair of extra arguments consisting of the parameter name and value. OMG_profiling.m contains an example of running the code. The config file contains definitions of most parameters.
  
  (4) The function OMG() returns output as one structure variable arranged as timeseries and timeslices, e.g. output.TimeSlice.PO4 gives you the PO4 fields at the times indicated to save.
  
  
  Useful Stuffs:
  
  Various functions for reorganising data from vectors to model grids are contained in general_functions.m; this is a work-around to get MATLAB to have multiple functions in a single file like FORTRAN subroutines or oject-oriented code. Calling gen_fcns=general_functions; will create a structure (gen_fcns or the name of your choice) which can be used to call "sub-functions" within like: gen_fcns.v2f(....)
  
  
  
=======
# omg
omg!
>>>>>>> 5d337c325a8cb524f7f55e42344e6cf3767b716f
