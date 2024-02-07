%% Initialise geochemistry

ECC=struct;

if bgc_pars.CARBCHEM_select
    
    ECC.constants=zeros(ocn_pars.nb,13);
    ECC.constants_interp=gchem_fcns.interpolate_constants(OCEAN,I,ocn_pars,gchem_pars);
    
    ECC.state=zeros(ocn_pars.nb,gen_pars.n_carbchem);
    ECC.state(:,I.H)=10e-8; % initial guess of H+
    
    % populate copy of H+
    gchem_pars.H=ECC.state(:,I.H);
    
    
    
end