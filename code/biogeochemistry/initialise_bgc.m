
bgc_pars.PO4_restore_data_file = gen_fcns.GetFullPath(bgc_pars.PO4_restore_data_file); % Convert relative to absolute path
%% Biogeochemistry Parameters

bgc_pars.rDOP_frac=1-bgc_pars.DOP_frac;                                     % fraction of POM exported    
if bgc_pars.DOP_k>0
    bgc_pars.DOP_k=(1/bgc_pars.DOP_k);                                          % reciprical and year
    bgc_pars.DOP_k=bgc_pars.DOP_k/360;
else
    bgc_pars.DOP_k=0;
end

% calculate POM remineralisation matrix if needed
if strcmp(bgc_pars.remin_scheme,'matrix')
    [POM_matrix , CaCO3_matrix]=bgc_fcns.create_remin_matrices(bgc_pars,ocn_pars);
end
bgc_pars.POM_matrix   = POM_matrix;
bgc_pars.CaCO3_matrix = CaCO3_matrix;

% stoichiometry - JDW: move this to tracer definitions as defines
% relationships between tracers?
bgc_pars.stoichiometry = zeros(ocn_pars.nb,gen_pars.n_bgc_tracers);

bgc_pars.stoichiometry(:,I.PO4)=1.0;
bgc_pars.stoichiometry(:,I.DOP)=0.0;
if(bgc_pars.CARBCHEM_select)
    bgc_pars.stoichiometry(:,I.DIC)=bgc_pars.C_to_P;
    bgc_pars.stoichiometry(:,I.ALK)=bgc_pars.ALK_to_P;
end
if(bgc_pars.O2_select)
    bgc_pars.stoichiometry(:,I.O2)=bgc_pars.O_to_P;
end
if(bgc_pars.Fe_cycle)
    bgc_pars.stoichiometry(:,I.TDFe)=bgc_pars.C_to_Fe;
    bgc_pars.stoichiometry(:,I.TL)=0;
end


% air-sea gas exchange constants
[bgc_pars.Sc_constants,bgc_pars.sol_constants]=bgc_fcns.initialise_gasex_constants(gen_pars,I);

% precalculate part of piston velocity
OCEAN.wspeed(:,:)=OCEAN.wspeed(:,:).^2 ...
        *bgc_pars.gastransfer_a ...
        .*OCEAN.seaice(:,:) ...                                              
        *0.01 ... % cm -> m
        *24; % hr -> day
    
OCEAN.PO4_obs = zeros(ocn_pars.nb,ocn_pars.n_A);                                              % BODGE !?!?!?!

% JDW - converting to model units of day-1
bgc_pars.u0PO4=bgc_pars.u0PO4/360;


%% Load forcing data

forcings=struct();
forcings=gen_fcns.load_forcing_data(gen_pars,bgc_pars,forcings,I,ocn_pars);


%%
