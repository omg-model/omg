%% BIOGEOCHEMISTRY FUNCTIONS
bgc_fcns                        =biogeochemistry_functions;

%% DEFAULT BIOGEOCHEMISTRY PARAMETERS

bgc_pars.uptake_scheme          = 'MM';                                   % Biological uptake scheme, default: Michaelis Menton
bgc_pars.remin_scheme           = 'matrix';                                % POM remineralisation scheme
bgc_pars.remin_function         = 'exponential';                           % remineralisation function

bgc_pars.u0PO4                  = 1.9582242E-06;                           % Michaelis-Menton maximum [PO4] uptake (mol kg-1 yr-1)
bgc_pars.KPO4                   = 2.1989611E-07;                           % Michaelis-Menton [PO4] half saturation concentration (mol kg-1);
bgc_pars.bio_tau                = 95.6337;                                 % biological uptake time-scale (days)
bgc_pars.bio_kT0                = 0.59;                                    % biological uptake temperature dependence scaling constant (dimensionless)
bgc_pars.bio_keT                = 15.8;                                    % biological uptake temperature dependence e-folding temperature (deg C)
bgc_pars.KFe                    = 0.1e-9;                                  % [Fe] half saturation (mol kg-1)
bgc_pars.PO4_restore_timescale  = 30;                                      % nutrient restoring timescale (days)
bgc_pars.PO4_restore_data_file  = '../experiments/PO4_OBS.mat';            % location of PO4 observations to restore to

bgc_pars.CARBCHEM_select        = false;                                   % flag to calculate carbonate chemistry and DIC/ALK tracers
bgc_pars.O2_select              = false;                                   % flag to calculate O2
bgc_pars.Fe_cycle               = false;                                   % flag to include Iron cycle

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

bgc_pars.det_Fe_sol             = 0.002014275;                             % aeolian Fe solubility
bgc_pars.det_Fe_sol_exp         = 0.500;                                   % aeolian Fe solubility exponent (use 1.0 for uniform solubility)
bgc_pars.scav_Fe_sf_POC         = 1.338130;                                % modifier of scavenging rate of dissolved Fe
bgc_pars.scav_Fe_sf_CaCO3       = 0.0;                                     % modifier of scavenging rate of dissolved Fe
bgc_pars.scav_fremin            = 0.0;                                     % scavenged regeneration
bgc_pars.scav_Fe_k0             = 0.079;                                   % initial scavenging rate (d-1; Parekh et al., 2005)
bgc_pars.scav_Fe_exp            = 0.58;                                    % (Parekh et al., 2005)
bgc_pars.K_FeL_pP               = 11.0;                                    % adjust pK' (FeL)
bgc_pars.FetoC_pP               = -0.4225;                                 % [FeT] dependent Fe:C ratio [Ridgwell, 2001] -- power
bgc_pars.FetoC_K                = 103684.0;                                % [FeT] dependent Fe:C ratio [Ridgwell, 2001] -- scaling
bgc_pars.FetoC_C                = 0.0;                                     % [FeT] dependent Fe:C ratio [Ridgwell, 2001] -- constant
bgc_pars.no_fsedFe              = false;                                   % return POFe?
bgc_pars.fixed_FetoC            = false;                                   % Variable Fe:C?
bgc_pars.conv_Fe_g_mol          = 17.86e-3;                                % convert g Fe to mol Fe
bgc_pars.det_Fe_frac            = 0.035;                                   % mass fraction of Fe in aeolian dust deposition
bgc_pars.conv_det_mol_g         = 60;                                      % convert mol det to g det 

bgc_pars.C_to_P                 = 106.0;                                   % Stoichiometric ratio of C:P 
bgc_pars.N_to_P                 = 16.0;                                    % Stoichiometric ratio of N:P 
bgc_pars.O_to_P                 = -170.0;                                  % Stoichiometric ratio of O:P 
bgc_pars.ALK_to_P               = -16.0;                                   % Stoichiometric ratio of alkalinity:P (-N:P)
bgc_pars.Fe_to_C                = 250000.0;                                % max Fe:C organic matter ratio for variable Fe:C

bgc_pars.gastransfer_a          = 0.31;                                    % gas transfer

bgc_pars.restore_ocnatm_Cinv    = false;                                   % restore ocean atmosphere inventory of carbon to initial value of run?

% bgc_pars.restore_pCO2_val=1;                                               % global modifier for restoring atm pCO2
% bgc_pars.restore_pO2_val=1;                                                % global modifier for restoring atm pO2
% bgc_pars.restore_PO4_val=1;                                                % global modifier for restoring ocn PO4
% bgc_pars.restore_DOP_val=1;                                                % global modifier for restoring ocn DOP
% bgc_pars.restore_DIC_val=1;                                                % global modifier for restoring ocn DIC
% bgc_pars.restore_ALK_val=1;                                                % global modifier for restoring ocn ALK
% bgc_pars.restore_O2_val=1;                                                 % global modifier for restoring ocn O2
% bgc_pars.restore_CaCO3_val=1; 
% bgc_pars.restore_TDFe_val=1;                                                % global modifier for restoring ocn Fe
% bgc_pars.restore_TL_val=1;                                                 % global modifier for restoring ocn TL
% bgc_pars.restore_Det_val=1;
% bgc_pars.restore_POP_val=1;
% bgc_pars.restore_POC_val=1;
% bgc_pars.restore_POFe_val=1;
% bgc_pars.restore_DOC_val=1;
% bgc_pars.restore_DOFe_val=1;
% 
% bgc_pars.restore_pCO2_timescale=1;                                         % timescale (year) for restoring atm pCO2
% bgc_pars.restore_pO2_timescale=1;                                          % timescale (year) for restoring atm pO2
% bgc_pars.restore_PO4_timescale=1;                                          % timescale (year) for restoring ocn PO4
% bgc_pars.restore_DOP_timescale=1;                                          % timescale (year) for restoring ocn DOP
% bgc_pars.restore_DIC_timescale=1;                                          % timescale (year) for restoring ocn DIC
% bgc_pars.restore_ALK_timescale=1;                                          % timescale (year) for restoring ocn ALK
% bgc_pars.restore_CaCO3_timescale=1; 
% bgc_pars.restore_O2_timescale=1;                                           % timescale (year) for restoring ocn O2
% bgc_pars.restore_TDFe_timescale=1;                                           % timescale (year) for restoring ocn TFe
% bgc_pars.restore_TL_timescale=1;                                           % timescale (year) for restoring ocn TL
% bgc_pars.restore_Det_timescale=1;
% bgc_pars.restore_POP_timescale=1;
% bgc_pars.restore_POC_timescale=1;
% bgc_pars.restore_POFe_timescale=1;
% bgc_pars.restore_DOC_timescale=1;
% bgc_pars.restore_DOFe_timescale=1;
% 
% bgc_pars.force_pCO2_val=1;                                                 % global modifier for forcing atm pCO2
% bgc_pars.force_pO2_val=1;                                                  % global modifier for forcing atm pO2
% bgc_pars.force_PO4_val=1;                                                  % global modifier for forcing ocn PO4
% bgc_pars.force_DOP_val=1;                                                  % global modifier for forcing ocn DOP
% bgc_pars.force_DIC_val=1;                                                  % global modifier for forcing ocn DIC
% bgc_pars.force_ALK_val=1;                                                  % global modifier for forcing ocn ALK
% bgc_pars.force_CaCO3_val=1;                             
% bgc_pars.force_O2_val=1;                                                   % global modifier for forcing ocn O2
% bgc_pars.force_TDFe_val=1;                                                  % global modifier for forcing ocn TFe
% bgc_pars.force_TL_val=1;                                                   % global modifier for forcing ocn TL
% bgc_pars.force_Det_val=1;
% bgc_pars.force_POP_val=1;
% bgc_pars.force_POC_val=1;
% bgc_pars.force_POFe_val=1;
% bgc_pars.force_DOC_val=1;
% bgc_pars.force_DOFe_val=1;


bgc_pars.PO4_init=2.159/1e6;                                               % initial PO4 (umol kg-1)    (Ridgwell et al. 2007)
bgc_pars.DOP_init=0.0235/1e6;                                              % inital DOP (umol kg-1)      (Najjar et al. 2007)
bgc_pars.O2_init=169.6/1e6;
bgc_pars.DIC_init=2244/1e6;
bgc_pars.ALK_init=2363/1e6;
bgc_pars.pCO2_init=278/1e6;
bgc_pars.pO2_init=0.2095;
bgc_pars.DOC_init=0.0235*106/1e6; 
bgc_pars.TDFe_init=0.650E-09;
bgc_pars.TL_init=1.000E-09;
bgc_pars.DOFe_init=1.00E-11;