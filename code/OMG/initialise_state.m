%% Set Initial Ocean State (from OMG restart or from GEnIE)

%% Initialise Tracers Array

TRACERS=zeros(ocn_pars.nb,gen_pars.n_tracers);     % tracer array at previous time step
PARTICLES=zeros(ocn_pars.nb,gen_pars.n_particles); % particle array
CARBCHEM=zeros(ocn_pars.nb,gen_pars.n_carbchem);   % carbonate chemistry array
ATM=zeros(1,gen_pars.n_atm);                       % atmosphere array

if strcmp(bgc_pars.uptake_scheme,'eco')
%     age    = sparse(zeros(size(TRACERS(ocn_pars.Ib,I.PHY))));
    bioinf.GENOME = zeros(ocn_pars.ni,eco_pars.jpmax,eco_pars.ngenes);
    bioinf.RGB    = zeros(ocn_pars.ni,eco_pars.jpmax,eco_pars.nrgb);

%     bioinf.RGB(:,round(eco_pars.jpmax/2),1) = double(ocn_pars.i(ocn_pars.Ib))./18.*1e6;
%     bioinf.RGB(:,round(eco_pars.jpmax/2),2) = double(ocn_pars.j(ocn_pars.Ib))./18.*1e6;
%     bioinf.RGB(:,round(eco_pars.jpmax/2),3) = 1e6;
else
%     parameters.eco_pars.nrgb   = 0;
%     parameters.eco_pars.ngenes = 0;
    bioinf = [];
end

% restart from OMG run
if ~isempty(gen_pars.OMG_restart_file)
    
    disp(['Attempting to restart from:  ''' gen_pars.OMG_restart_file ''''])
    
    load(['../../output/' gen_pars.OMG_restart_file '/Restart.mat']);

    disp(['Restart data loaded from year ' num2str(yr) ' of ' gen_pars.OMG_restart_file])
    
    % Check restart parameters match current parameters
    gen_fcns.check_restart_setup(varargin,gen_pars.OMG_restart_file);
    


elseif gen_pars.netcdf_restart % restart from netcdf files 
    
    disp(['Restarting from matrix netCDF files'])
    
    [dum_PO4,dum_DOP,dum_O2,dum_DIC,dum_ALK]=gen_fcns.read_genie_netcdf([gen_pars.TM_path '/fields_biogem_3d.nc'],1,v_index,'ocn_PO4','ocn_DOM_P','ocn_O2','ocn_DIC','ocn_ALK');
    TRACERS(:,I.PO4)=dum_PO4;
    TRACERS(:,I.DOP)=dum_DOP;
    
    if(bgc_pars.O2_select)
        TRACERS(:,I.O2)=dum_O2;
    end
    if(bgc_pars.CARBCHEM_select)
        TRACERS(:,I.DIC)=dum_DIC;
        TRACERS(:,I.ALK)=dum_ALK;
        tmp=gen_fcns.read_genie_netcdf([gen_pars.TM_path '/fields_biogem_2d.nc'],1,v_index,'atm_pCO2');
        ATM(1,I.pCO2)=tmp(ocn_pars.Ib);
        CARBCHEM(:,I.H)=10e-8; % temporary! JDW: remove?
    end
    
    if strcmp(bgc_pars.uptake_scheme,'eco')
        TRACERS(ocn_pars.Ib,I.PHY)=eco_pars.PHY_init;
        
        % seed locations with biomass if set in config file
        if ~isempty(eco_pars.seed_PHY) % specify plankton with initial seed population
            seed_PHY=I.PHY(eco_pars.seed_PHY);
        else % all plankton have initial seed population
            seed_PHY=I.PHY;
        end
        if ~isempty(eco_pars.seed_loc) % if seeding locations specified
            TRACERS(eco_pars.seed_loc,seed_PHY)=eco_pars.seed_val;
        else % else seed everywhere
            TRACERS(ocn_pars.Ib,seed_PHY)=eco_pars.seed_val;
        end
    end
     

else % no restart selected, start from initial conditions
    
    disp(['Initialising OMG from initial conditions'])
    
    TRACERS(:,I.PO4)=bgc_pars.PO4_init;                                     % initial PO4 (umol kg-1 -> mol kg-1) (Ridgwell et al. 2007)
    TRACERS(:,I.DOP)=bgc_pars.DOP_init;                                     % inital DOP umol kg-1 -> mol kg-1) (Najjar et al. 2007)
    if bgc_pars.O2_select
        TRACERS(:,I.O2)=bgc_pars.O2_init;
        ATM(1,I.pO2)=bgc_pars.pO2_init;
    end
    if bgc_pars.CARBCHEM_select
        TRACERS(:,I.DIC)=bgc_pars.DIC_init;
        TRACERS(:,I.ALK)=bgc_pars.ALK_init;
        TRACERS(:,I.DOC)=bgc_pars.DOC_init;
        ATM(I.pCO2)=bgc_pars.pCO2_init;
        CARBCHEM(:,I.H)=10e-8; % JDW: remove?
    end
    if strcmp(bgc_pars.uptake_scheme,'eco')
        TRACERS(ocn_pars.Ib,I.PHY)=eco_pars.PHY_init;                            % initial plankton concentrations
        if isempty(eco_pars.seed_val)
            eco_pars.seed_val=0.0;                                              % dead ocean!
        end
        % seed locations with biomass if set in config file
        if ~isempty(eco_pars.seed_PHY) % specify plankton with initial seed population
            seed_PHY=I.PHY(eco_pars.seed_PHY);
        else % all plankton have initial seed population
            seed_PHY=I.PHY;
        end
        if ~isempty(eco_pars.seed_loc) % if seeding locations specified
            TRACERS(ocn_pars.Ib(eco_pars.seed_loc),seed_PHY)=eco_pars.seed_val;
        else % else seed everywhere
            TRACERS(ocn_pars.Ib,seed_PHY)=eco_pars.seed_val;
        end
    end
    
end


% initial prediction of carbonate chemistry state to speed up convergence
parameters.ocn_pars.T=functions.OMG_fcns.interpolate_var_at_t(360,parameters.ocn_pars,OCEAN.T);
parameters.ocn_pars.S=functions.OMG_fcns.interpolate_var_at_t(360,parameters.ocn_pars,OCEAN.S);
ECC = functions.gchem_fcns.calc_carbonate_constants(ECC,parameters,functions,true,360);
ECC = functions.gchem_fcns.solve_carbonate_system(ECC,TRACERS,parameters,true);

% get initial ocean and atmosphere carbon inventory
if parameters.bgc_pars.restore_ocnatm_Cinv
    parameters.bgc_pars.ocnatm_Cinv_target=sum(TRACERS(:,parameters.ind_pars.DIC).*parameters.ocn_pars.M) + ATM(parameters.ind_pars.pCO2)*parameters.ocn_pars.atm_mol;
end
