%% 16 level seasonal preindustrial 

% Transport Matrix
gen_pars.TM_path='../../TM_data/worjh2_seasonal';                           % path to TM files

% Model Setup
gen_pars.TM_dt=1/96;                                                        % TM timestep (year)
gen_pars.dt_ratio = 1;                                                      % biogeochemistry to circulation timestep ratio

% Biogeochemistry
bgc_pars.CARBCHEM_select=false;                                             % flag to calculate carbonate chemistry and DIC/ALK tracers
bgc_pars.O2_select=false;                                                   % flag to calculate O2

bgc_pars.uptake_scheme='MM';                                                % Biological uptake scheme, default: Michaelis Menton
bgc_pars.remin_scheme='matrix';                                             % POM remineralisation scheme
bgc_pars.remin_function='exponential';                                      % POM remineralisation function

bgc_pars.u0PO4=8.9876e-006;                                                 % Michaelis-Menton maximum PO4 uptake (mol kg-1 yr-1)
bgc_pars.KPO4=8.9369e-007;                                                  % Michaelis-Menton half saturation concentration (mol kg-1);

bgc_pars.DOP_frac=0.66;                                                     % fraction of export for DOP  
bgc_pars.DOP_k=0.5;                                                         % DOP decay rate (yr))

bgc_pars.POC_eL1=589.9451;                                                  % efolding depth of labile POP (m)    
bgc_pars.POC_eL2=1000000.0;                                                 % efolding depth of refractory POP (m)    
bgc_pars.POC_frac2=0.0557;                                                  % initial fraction of 'refractory' POP

bgc_pars.red_POC_CaCO3=0.0485;                                              % underlying ratio of CaCO3 to POC
bgc_pars.POC_CaCO3_pP=0.7440;                                               % exponent for modifier of CaCO3:POC ratio

bgc_pars.CaCO3_frac2=0.45;                                                  % initial fraction of "refractory" CaCO3 
bgc_pars.CaCO3_eL1=1.8905e+003;                                             % e-folding depth of CaCO3
bgc_pars.CaCO3_eL2=1000000.0;    




