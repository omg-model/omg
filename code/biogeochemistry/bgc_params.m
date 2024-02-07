%% BIOGEOCHEMISTRY FUNCTIONS
bgc_fcns                        =biogeochemistry_functions;

%% DEFAULT BIOGEOCHEMISTRY PARAMETERS

bgc_pars.uptake_scheme          = 'MM';                                   % Biological uptake scheme, default: Michaelis Menton
bgc_pars.remin_scheme           = 'matrix';                                % POM remineralisation scheme
bgc_pars.remin_function         = 'exponential';                           % remineralisation function

bgc_pars.u0PO4                  = 1.9582242E-06;                           % Michaelis-Menton maximum PO4 uptake (mol kg-1 yr-1)
bgc_pars.KPO4                   = 2.1989611E-07;                           % Michaelis-Menton half saturation concentration (mol kg-1);
bgc_pars.PO4_restore_timescale  = 30;                                      % nutrient restoring timescale (days)
bgc_pars.PO4_restore_data_file  = '../experiments/PO4_OBS.mat';            % location of PO4 observations to restore to

bgc_pars.CARBCHEM_select        = false;                                   % flag to calculate carbonate chemistry and DIC/ALK tracers
bgc_pars.O2_select              = false;                                   % flag to calculate O2

bgc_pars.parfrac                = 0.43;                                    % PAR fraction

bgc_pars.DOP_frac               = 0.66;                                    % fraction of export for DOP  
bgc_pars.DOP_k                  = 0.5;                                     % DOP decay rate (yr))

bgc_pars.POC_eL1                = 550.5195;                                % efolding depth of labile POP (m)    
bgc_pars.POC_eL2                = 1000000.0;                               % efolding depth of refractory POP (m)    
bgc_pars.POC_frac2              = 6.4591110E-02;                           % initial fraction of 'refractory' POP

bgc_pars.PIC_POC                = 0.044372;                                % underlying ratio of CaCO3 to POC
bgc_pars.PIC_POC_omega_mod      = 0.8053406;                               % exponent for modifier of CaCO3:POC ratio

bgc_pars.CaCO3_frac2            = 0.4325;                                  % initial fraction of "refractory" CaCO3 
bgc_pars.CaCO3_eL1              = 1083.486;                                % e-folding depth of CaCO3
bgc_pars.CaCO3_eL2              = 1000000.0; 

bgc_pars.C_to_P                 = 106.0;                                   % Stoichiometric ratio of C:P 
bgc_pars.N_to_P                 = 16.0;                                    % Stoichiometric ratio of N:P 
bgc_pars.O_to_P                 = -170.0;                                  % Stoichiometric ratio of O:P 
bgc_pars.ALK_to_P               = -16.0;                                   % Stoichiometric ratio of alkalinity:P (-N:P)

bgc_pars.gastransfer_a          = 0.31;                                    % gas transfer

bgc_pars.restore_ocnatm_Cinv    = false;                                   % restore ocean atmosphere inventory of carbon to initial value of run?

bgc_pars.restore_pCO2_val=1;                                               % global modifier for restoring atm pCO2
bgc_pars.restore_pO2_val=1;                                                % global modifier for restoring atm pO2
bgc_pars.restore_PO4_val=1;                                                % global modifier for restoring ocn PO4
bgc_pars.restore_DOP_val=1;                                                % global modifier for restoring ocn DOP
bgc_pars.restore_DIC_val=1;                                                % global modifier for restoring ocn DIC
bgc_pars.restore_ALK_val=1;                                                % global modifier for restoring ocn ALK
bgc_pars.restore_O2_val=1;                                                 % global modifier for restoring ocn O2

bgc_pars.restore_pCO2_timescale=1;                                         % timescale (year) for restoring atm pCO2
bgc_pars.restore_pO2_timescale=1;                                          % timescale (year) for restoring atm pO2
bgc_pars.restore_PO4_timescale=1;                                          % timescale (year) for restoring ocn PO4
bgc_pars.restore_DOP_timescale=1;                                          % timescale (year) for restoring ocn DOP
bgc_pars.restore_DIC_timescale=1;                                          % timescale (year) for restoring ocn DIC
bgc_pars.restore_ALK_timescale=1;                                          % timescale (year) for restoring ocn ALK
bgc_pars.restore_O2_timescale=1;                                           % timescale (year) for restoring ocn O2

bgc_pars.force_pCO2_val=1;                                                 % global modifier for forcing atm pCO2
bgc_pars.force_pO2_val=1;                                                  % global modifier for forcing atm pO2
bgc_pars.force_PO4_val=1;                                                  % global modifier for forcing ocn PO4
bgc_pars.force_DOP_val=1;                                                  % global modifier for forcing ocn DOP
bgc_pars.force_DIC_val=1;                                                  % global modifier for forcing ocn DIC
bgc_pars.force_ALK_val=1;                                                  % global modifier for forcing ocn ALK
bgc_pars.force_O2_val=1;                                                   % global modifier for forcing ocn O2

bgc_pars.PO4_init=2.159/1e6;                                               % initial PO4 (umol kg-1)    (Ridgwell et al. 2007)
bgc_pars.DOP_init=0.0235/1e6;                                              % inital DOP (umol kg-1)      (Najjar et al. 2007)
bgc_pars.O2_init=169.6/1e6;
bgc_pars.DIC_init=2244/1e6;
bgc_pars.ALK_init=2363/1e6;
bgc_pars.pCO2_init=278/1e6;
bgc_pars.pO2_init=0.2095;