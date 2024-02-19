%% 8 level annual preindustrial 

% Transport Matrix
gen_pars.TM_path='../../TM_data/worbe2';                                  % path to TM files

% Model Setup
gen_pars.TM_dt=1/96;                                                        % TM timestep (year)
gen_pars.dt_ratio = 1;                                                      % biogeochemistry to circulation timestep ratio

% Biogeochemistry
bgc_pars.CARBCHEM_select=false;                                             % flag to calculate carbonate chemistry and DIC/ALK tracers
bgc_pars.O2_select=false;                                                   % flag to calculate O2

bgc_pars.uptake_scheme='MM';                                                % Biological uptake scheme, default: Michaelis Menton
bgc_pars.remin_scheme='matrix';                                             % POM remineralisation scheme
bgc_pars.remin_function='exponential';                                      % POM remineralisation function
