%% GENERAL FUNCTIONS
gen_fcns                            = general_functions;

%% DEFAULT GENERAL PARAMETERS

gen_pars.TM_path                    = '../../TM_data/worbe2';              % default configuration
gen_pars.TM_correction              = false;                               % apply correction to TMs (for MMM)

gen_pars.n_dt                       = 96;                                  % number timesteps per year


% gen_pars.dt_ratio needs  to be made robust so that it is not used if Matlab ode
% solvers are used (not compatible with variable timestep)
% This is currently done in OMG.m (parameters.gen_pars.dt_ratio=1;)
gen_pars.dt_ratio                   = 1;                                   % biogeochemistry to circulation timestep ratio

gen_pars.start_year                 = 1;                                   % start year of run

gen_pars.forcings_directory          ='';                                  % name of directory containing forcings

gen_pars.progress_reporting         = false;                               % print %run completed to screen?

% move below to initialise_gen.m as don't need to be user-definable?
gen_pars.conv_d_s                   = 24*3600;   
gen_pars.conv_d_yr                  = 365.25;
gen_pars.conv_yr_s                  = 365.25*24.0*3600.0;

gen_pars.save2netcdf                = true;                                % flag to save netcdf files
gen_pars.save2text                  = true;                                % flag to save text timeseries files
gen_pars.save2matobj                = false;                               % flag to save to MatFile objects


gen_pars.netcdf_restart              = false;                              % flag to restart from netcdf files
gen_pars.OMG_restart_file           = '';                                  % previous OMG run to restart from
gen_pars.OMG_restart_year           = -1;                                  % previous OMG run year to restart from (-1 is latest available year)
gen_pars.tseries_index              = [];                                  % no timeseries sites by default

gen_pars.integrate_scheme            = 'fwd';                             % integration scheme ('fwd','ode_45','ode23') 

% ODE solver info
% set solver options: positive definite, tolerances, min/max tsteps
gen_pars.odeoptions.AbsTol          = 1e-9;
gen_pars.odeoptions.RelTol          = 1e-6;
gen_pars.odeoptions.NormControl     = 'off';
gen_pars.odeoptions.MaxStep         = 360;
gen_pars.odeoptions.InitialStep     = 1; % day
gen_pars.odeoptions.Stats           = 'off';

% I/O
gen_pars.save_output_directory      = '';                                   % name of output directory - blank ('') = automatic name generated
gen_pars.save_timeslice_output =[1];                                       % timeslice output years 
gen_pars.save_timeseries_output=[1];                                       % timeseries output years 
gen_pars.save_matlab_output=false;                                          % additional save to matlab variables?

gen_pars.save_intra_annual_n        = 1;                                   % number of intra-annual saving intervals
gen_pars.n_sub_tsteps               = 1;                                    % number of extra time steps per standard timestep

gen_pars.save_carbchem              = false;                                % save carbonate chemistry output?
% CURRENTLY NOT SET UP FOR INTRA-ANNUAL OUTPUT
                                                                            %   1 = yearly
                                                                            %   4 = seasonal
                                                                            %  12 = monthly
                                                                            %  52 = weekly
                                                                            % 365 = daily
