%% ECOSYSTEM & EVOLUTION FUNCTIONS
eco_fcns                = ecoevo_functions;

%% DEFAULT ECOSYSTEM & EVOLUTION PARAMETERS

eco_pars.plankton_advect = '2D'; % plankton advection scheme: 3D (advected out of surface layer)
                % or 2D (horizontal transport only)
% Trait dimensions
eco_pars.nsize      = 10; % total number of plankton size classes
eco_pars.nTopt      = 10; % total number of plankton thermal optima
eco_pars.ntroph     =  1; % total number of plankton trophic classes

% thermal traits
eco_pars.Twidth =  6.00;
eco_pars.Tslope =  0.05;
eco_pars.Tref   = 20.00; 

eco_pars.alpha_a    = +3981;     % m^3/mmol P/d                            % Edwards et al. 2012
eco_pars.alpha_b    = -0.55;    % 0.63;     %    

eco_pars.linearmort = +0.025;	% /d                                       % Maranon et al EcoLett 2013

% Light and Photosynthesis
eco_pars.kw         = 0.04;     % light attenuation by water
eco_pars.kchl       = 0.03;     % light attenuation by chlorophyll
eco_pars.kPAR       = 40;       % PAR 1/2 saturation (actually "1/sqrt(2) saturation")

% grazing
eco_pars.clearance_a   = +103.68;	% m3 / mmol P / d                         % coefficient of clearance rate
eco_pars.clearance_b   = -0.16;                                               % exponent of clearance rate

eco_pars.gmax_a        = +21.9;                                               % coefficient of maximum grazing rate
eco_pars.gmax_b        = -0.16;                                               % exponent of maximum grazing rate

eco_pars.refuge        = -100;                                                % grazing refuge parameter (negative)
                                                                           
eco_pars.pp_opt     = +1000;                                               % optimum pred:prey ratio
eco_pars.pp_sig     = +2;%log(10);                                         % width grazing kernel
eco_pars.lambda     = +0.7;                                                 % grazing efficiency

% trophic trade off
eco_pars.tau        = 1.0;      % trophic trade-off parameter (1 = linear)

%% initialisation
eco_pars.PHY_init   = 0.0;                                                 % initial biomass *everywhere*
eco_pars.seed_loc   = [];                                                  % spatial indices for seeding
eco_pars.seed_PHY   = [];                                                  % plankton indices for seeding
eco_pars.seed_val   = 1e-12;                                               % biomass for seeding

% interference function
eco_pars.sigma_interference_size    = log10(2); % interference distance in size dimension 
eco_pars.sigma_interference_troph   = 0.01;     % interference distance in trophic dimension 
eco_pars.interference_strength      = 0;        % interference strength

% evolution (trophic strategy, size)
eco_pars.mutrat     = [1e-15 1e-15]; % mutation rate relative to size of trait space (will be normalised)

% Extinction thresholds
eco_pars.functional_extinction   = 1e-15;

% plankton advection scaling
eco_pars.trscale=1;

eco_pars.eco_tstep              = 0.25;             

eco_pars.ngenes  = 0;
eco_pars.nrgb    = 0;

%% quota (used to calculate mumax)
eco_pars.muinf_a    = +4.70;	% /d                                        % Ward et al. AmNat 2017
eco_pars.muinf_b    = -0.26; 	%                                           % Ward et al. AmNat 2017

eco_pars.Vmax_a     = +0.024;   % pgN/cell/d                               % Ward et al. AmNat 2017
eco_pars.Vmax_b     = +1.10;	%   

eco_pars.Qmin_a     = +0.032;   % pgN/cell                                 % Ward et al. AmNat 2017
eco_pars.Qmin_b     = 0.76;  	%                                          % Ward et al. AmNat 2017


