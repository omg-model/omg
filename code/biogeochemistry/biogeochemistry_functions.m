function [ bgc_fcns ] = biogeochemistry_functions ( )
    
    bgc_fcns.transportmatrix = @transportmatrix;
    bgc_fcns.SurfaceProd=@SurfaceProd;
    bgc_fcns.remin_DOM=@remin_DOM; 
    bgc_fcns.remin_POM=@remin_POM;
    bgc_fcns.remin_CaCO3=@remin_CaCO3;
    bgc_fcns.create_remin_matrices=@create_remin_matrices;
    bgc_fcns.airsea_gas_exchange=@airsea_gas_exchange;
    bgc_fcns.calc_gasexchange_constants=@calc_gasexchange_constants;
    bgc_fcns.atm_forcings=@atm_forcings;
    bgc_fcns.ocn_forcings=@ocn_forcings;
    bgc_fcns.sed_forcings=@sed_forcings;
    bgc_fcns.restore_ocnatm_Cinv=@restore_ocnatm_Cinv;
    bgc_fcns.initialise_gasex_constants=@initialise_gasex_constants;
    bgc_fcns.aeolian_Fe=@aeolian_Fe;
    bgc_fcns.Fe_speciation=@Fe_speciation;
    bgc_fcns.scavenge_Fe=@scavenge_Fe;
    bgc_fcns.calc_variable_Fe_to_P=@calc_variable_Fe_to_P;

end

%%----- SUBROUTINES -----%%
%%
function [dCdt,POM_prodn] = SurfaceProd(dCdt , SURFACE , parameters)
    % BAW changed i/o format to pass dCdt and Particle_Export (implicit particulate production)
    gen_pars=parameters.gen_pars;
    bgc_pars=parameters.bgc_pars;
    ocn_pars=parameters.ocn_pars;
    I=parameters.ind_pars;
    
    Ib       = ocn_pars.Ib;
        
    solfor=parameters.ocn_pars.solfor(Ib);    
    seaice=parameters.ocn_pars.seaice(Ib);
    
    % calculate uptake of PO4
    switch bgc_pars.uptake_scheme
        case 'MM' 
            % Ridgwell et al., (2007) Biogeosciences
            uptake=(seaice.*solfor.*(bgc_pars.u0PO4.*(SURFACE(:,I.PO4)./(SURFACE(:,I.PO4)+bgc_pars.KPO4))));
        case 'MM_Fe'
            % van de Velde et al., (2021) GMD
            FT=bgc_pars.bio_kT0.*exp(parameters.ocn_pars.T(Ib)./bgc_pars.bio_keT);
            FN=min([SURFACE(:,I.PO4)./(SURFACE(:,I.PO4)+bgc_pars.KPO4),SURFACE(:,I.TDFe)./(SURFACE(:,I.TDFe)+bgc_pars.KFe)],[],2);
            Nmin=min([SURFACE(:,I.PO4),bgc_pars.C_to_P.*SURFACE(:,I.TDFe)./bgc_pars.C_to_Fe],[],2);
            uptake=solfor.*seaice.*FT.*FN.*Nmin./bgc_pars.bio_tau;
        case 'restore'
             error('Restore scheme broken by Ben! Need to add PO4_obs to parameter structure.')
%              uptake=(seaice(Ib)).*((PO4(Ib)-PO4_obs(Ib,dt_yr))*bgc_pars.PO4_restore_timescale);
%              uptake(uptake<0)=0;
        case 'strangelove'
            uptake=0.0;
        otherwise
            error('No PO4 uptake routine selected')
    end

    % calculate variable Fe:C as f(TDFe)
    if bgc_pars.Fe_cycle
        bgc_pars.uptake_stoichiometry(Ib,I.TDFe)=calc_variable_Fe_to_P( SURFACE , parameters );
    end
    
    uptake=uptake.*bgc_pars.uptake_stoichiometry(Ib,:);
    dCdt(Ib,1:gen_pars.n_bgc_tracers) = dCdt(Ib,1:gen_pars.n_bgc_tracers)  - uptake;
    dCdt(Ib,1:gen_pars.n_bgc_tracers) = dCdt(Ib,1:gen_pars.n_bgc_tracers)  + uptake.*bgc_pars.DOP_frac*bgc_pars.mapOCN_DOM;
    POM_prodn = uptake.*bgc_pars.rDOP_frac*bgc_pars.mapOCN_POM; 
    
    % Calculate uptake of all BGC tracers from PO4 uptake, using bgc_pars.stoichiometry
    %dCdt(:,1:gen_pars.n_bgc_tracers) = -uptake .* bgc_pars.uptake_stoichiometry(Ib,:);
    % Add dissolved organic matter production to dCdt(:,I.DOP) 
    %dCdt(:,I.DOP) = dCdt(:,I.DOP) + uptake.*bgc_pars.DOP_frac; 
    % Output implicit particulate organic matter production as dPOPdt
    %POM_prodn        = uptake.*bgc_pars.rDOP_frac; 
    
%         
%         rainratio_exponent=(CARBCHEM(Ib,I.SAT_CA)-1).^bgc_pars.red_PIC_POC_mod;
%         rainratio=bgc_pars.red_PIC_POC.*rainratio_exponent;
%         CaCO3=(uptake*bgc_pars.rDOP_frac*bgc_pars.red_P_to_C).*rainratio;
%         
%         PARTICLES(Ib,I.CaCO3)=PARTICLES(Ib,I.CaCO3)+CaCO3;
%         Jpump(Ib,I.DIC)=Jpump(Ib,I.DIC)-CaCO3;
%         Jpump(Ib,I.ALK)=Jpump(Ib,I.ALK)-2*CaCO3; % adjust ALK and DIC for formation of CaCO3 in 2:1 ratio
%     end
%     
%     % update oxygen
%     if bgc_pars.O2_select        
%         Jpump(Ib,I.O2)=Jpump(Ib,I.O2)+uptake*bgc_pars.red_P_to_O;        
%     end
 
end

%%
function [dCdt] = remin_DOM(TRACERS , dCdt , parameters)

    gen_pars = parameters.gen_pars;
    bgc_pars = parameters.bgc_pars;
    I        = parameters.ind_pars;
        
    % remineralisation of DOM
    %DOM_remin=TRACERS(:,I.DOP) .* bgc_pars.DOP_k;
    % first calculate remineralisation for *everything*, then map remin of DOM
    % tracers back to normal tracers
    remin = TRACERS .* bgc_pars.DOP_k * bgc_pars.mapOCN_DOM';
    
    % add effects of organic matter remineralisation stoichiometry
    remin = remin + remin(:,I.PO4) .* bgc_pars.remin_stoichiometry; 

    % update tendencies 
    dCdt = dCdt + remin - remin*bgc_pars.mapOCN_DOM; 
    
    % subtract remin from DOP column
    %dCdt(:,I.DOP) = dCdt(:,I.DOP) - DOM_remin;
    % add remin to BGC columns according to bgc_pars.stoichiometry
    %dCdt(:,1:gen_pars.n_bgc_tracers) = dCdt(:,1:gen_pars.n_bgc_tracers) + DOM_remin.*bgc_pars.stoichiometry; 

end


%%

function [dCdt , PARTICLES] = remin_POM ( dCdt , POM_prodn , PARTICLES , parameters )

    gen_pars = parameters.gen_pars;
    bgc_pars = parameters.bgc_pars;
    ocn_pars = parameters.ocn_pars;
    I        = parameters.ind_pars;

    i_remin=1:gen_pars.n_bgc_tracers;

    % Take POM production from the surface layer and redistribute across water column
    POM_remin=bgc_pars.POM_matrix(:,ocn_pars.Ib)*POM_prodn;
    
    % get flux hitting seafloor
    POM_total_export = sum(POM_prodn.*ocn_pars.M(ocn_pars.Ib)); % mol
    for n=1:size(POM_remin,2)
        benthic_remin(:,n)=(POM_total_export(n)-accumarray(ocn_pars.wc,POM_remin(:,n).*ocn_pars.M)).*ocn_pars.rM(ocn_pars.Iben);
    end
    
    % reflective boundary at seafloor
    if bgc_pars.Fe_cycle
        benthic_flux(:,I.POFe)=0; % except for particulate-bound Fe
    end
    POM_remin(ocn_pars.Iben) = POM_remin(ocn_pars.Iben) + benthic_remin;

    % 
    dCdt(:,i_remin) = dCdt(:,i_remin) + benthic_remin .* bgc_pars.stoichiometry;





    % add remin to dCdt
    %dCdt(:,i_remin) = dCdt(:,i_remin) + POM_remin .* bgc_pars.stoichiometry;

    % get flux hitting seafloor
    benthic_remin=(1-accumarray(ocn_pars.wc,POM_remin.*ocn_pars.M).*ocn_pars.rM(ocn_pars.Ib));
    benthic_remin=benthic_remin.*ocn_pars.M(ocn_pars.Ib).*ocn_pars.rM(ocn_pars.Iben);

    % reflective boundary at seafloor
    if bgc_pars.Fe_cycle
        benthic_flux(:,I.TDFe)=0; % except for particulate-bound Fe
    end
    dCdt(:,i_remin) = dCdt(:,i_remin) + benthic_remin .* bgc_pars.stoichiometry;

    % Add to particle export matrix
    PARTICLES(:,I.POM) = POM_remin;
    
end

%%
function [dCdt , PARTICLES] = remin_CaCO3(dCdt , POM_prodn , ECC , PARTICLES , parameters)

if parameters.bgc_pars.CARBCHEM_select

    bgc_pars = parameters.bgc_pars;
    ocn_pars = parameters.ocn_pars;
    I        = parameters.ind_pars;
    Ib       = parameters.ocn_pars.Ib;
    
    % Rain-ratio
    if bgc_pars.PIC_POC_omega_mod>0.0
        rainratio_exponent=(ECC.state(Ib,I.SAT_CA)-1).^bgc_pars.PIC_POC_omega_mod;
        rainratio_exponent(ECC.state(Ib,I.SAT_CA)<=1)=0.0;
    else 
        rainratio_exponent=0.0;
    end
    rainratio=bgc_pars.PIC_POC.*rainratio_exponent;
    
    % Calculate CaCO3 Production from POM production and rain ratio
    CaCO3_prodn = POM_prodn.*bgc_pars.C_to_P.*rainratio;
    
    % DIC and ALK uptake from CaCO3 production
    dCdt(Ib,I.DIC) = dCdt(Ib,I.DIC) - CaCO3_prodn;
    dCdt(Ib,I.ALK) = dCdt(Ib,I.ALK) - 2*CaCO3_prodn; % adjust ALK for formation of CaCO3 in 2:1 ratio
    
    % Take CaCO3 production from the surface layer and redistribute across other layers
    CaCO3_remin=bgc_pars.CaCO3_matrix(:,ocn_pars.Ib)*CaCO3_prodn;
    
    % DIC and ALK production from CaCO3 dissolution
    dCdt(:,I.DIC) = dCdt(:,I.DIC) + CaCO3_remin;
    dCdt(:,I.ALK) = dCdt(:,I.ALK) + 2*CaCO3_remin;
    
    % Add to particle matrix
    PARTICLES(:,I.CaCO3) = CaCO3_remin;

end

end

%% ]
function [ POM_matrix , CaCO3_matrix ] = create_remin_matrices ( bgc_pars , ocn_pars )

dzt=ocn_pars.zt_edges;
production_max_k=1;


% n.b. this not correct if >1 depth layer
% can add larger number and code below removes unfilled entries
spcol_ind=zeros(ocn_pars.nb,1);
sprow_ind=zeros(ocn_pars.nb,1);
spval=zeros(ocn_pars.nb,1);

count=1; 


for n=1:max(ocn_pars.wc)
    % for each water-column
    ind=find(ocn_pars.wc==n);
    k_btm=numel(ind);
    % if within export boxes
    if k_btm<=production_max_k
        k_start=1;
        col_ind=ind(k_start:k_btm);
        remin=zeros(numel(col_ind),1);
        remin(end)=1.0; % remineralise in bottom-most gridbox
    % otherwise calculate full remineralisation profile    
    else
        k_start=production_max_k;
        col_ind=ind(1:production_max_k);

        % flux curves
        switch(bgc_pars.remin_function)
            case('exponential')
                flux=(1-bgc_pars.POC_frac2)*exp((dzt(production_max_k+1)-dzt)/bgc_pars.POC_eL1) + ...
                    bgc_pars.POC_frac2 *exp((dzt(production_max_k+1)-dzt)/bgc_pars.POC_eL2);
                 
            case('powerlaw')
                flux=(dzt./dzt(production_max_k+1)).^bgc_pars.powerlaw_b;
        end

        remin=-diff(flux);
        remin(1:production_max_k)=0.0;
        %remin(k_btm)=remin(k_btm)+(1-sum(remin(production_max_k:k_btm)));

    end

    % build sparse i,j,s arrays
    for r=k_start:k_btm
        for c=1:numel(col_ind)

            sprow_ind(count,1)=ind(r);
            spcol_ind(count,1)=col_ind(c);
            spval(count,1)=1.0*remin(r)*ocn_pars.M(col_ind(c))/ocn_pars.M(ind(r)); % volume weighted
            count=count+1;

        end
    end


end

spcol_ind(count:end)=[];
sprow_ind(count:end)=[];
spval(count:end)=[];

POM_matrix=sparse(sprow_ind,spcol_ind,spval,ocn_pars.nb,ocn_pars.nb);


if bgc_pars.CARBCHEM_select

    dzt=ocn_pars.zt_edges;
production_max_k=1;


% n.b. this not correct if >1 depth layer
% can add larger number and code below removes unfilled entries
spcol_ind=zeros(ocn_pars.nb,1);
sprow_ind=zeros(ocn_pars.nb,1);
spval=zeros(ocn_pars.nb,1);

count=1; 


for n=1:max(ocn_pars.wc)
    % for each water-column
    ind=find(ocn_pars.wc==n);
    k_btm=numel(ind);
    % if within export boxes
    if k_btm<=production_max_k
        k_start=1;
        col_ind=ind(k_start:k_btm);
        remin=zeros(numel(col_ind),1);
        remin(end)=1.0; % remineralise in bottom-most gridbox
    % otherwise calculate full remineralisation profile    
    else
        k_start=production_max_k;
        col_ind=ind(1:production_max_k);

        % flux curves
        flux=(1-bgc_pars.CaCO3_frac2)*exp((dzt(production_max_k+1)-dzt)/bgc_pars.CaCO3_eL1) + ...
            bgc_pars.CaCO3_frac2 *exp((dzt(production_max_k+1)-dzt)/bgc_pars.CaCO3_eL2);
                 
        remin=-diff(flux);
        remin(1:production_max_k)=0.0;
        remin(k_btm)=remin(k_btm)+(1-sum(remin(production_max_k:k_btm)));

    end

    % build sparse i,j,s arrays
    for r=k_start:k_btm
        for c=1:numel(col_ind)

            sprow_ind(count,1)=ind(r);
            spcol_ind(count,1)=col_ind(c);
            spval(count,1)=1.0*remin(r)*ocn_pars.M(col_ind(c))/ocn_pars.M(ind(r)); % volume weighted
            count=count+1;

        end
    end

end

spcol_ind(count:end)=[];
sprow_ind(count:end)=[];
spval(count:end)=[];

CaCO3_matrix=sparse(sprow_ind,spcol_ind,spval,ocn_pars.nb,ocn_pars.nb);

else
    CaCO3_matrix=1;
end


end


%%
function [ dCdt , dATMdt ] = airsea_gas_exchange ( dCdt , dATMdt , TRACERS , ATM , ECC , parameters )

    if parameters.gen_pars.n_atm>0

        ocn_pars=parameters.ocn_pars;
        bgc_pars=parameters.bgc_pars;
        gen_pars=parameters.gen_pars;
        I=parameters.ind_pars;
        Ib       =ocn_pars.Ib;
        
        % m2*kg / m3
        loc_A_rho=ocn_pars.A.*ocn_pars.rho(Ib);
       
        % loop over atmospheric tracers
        for n=1:gen_pars.n_atm
    
            % n = atm index, find corresponding OCN index
            ind_ocn=I.ATM_OCN{n}; 
            
            % calcaulte constants
            [Sc,Solubility]=calc_gasexchange_constants(parameters);
        
            % piston velocity cm hr-1 -> m day-1
            piston=ocn_pars.wspeed(Ib).*Sc(:,n);
    
            % fluxes (kg yr-1 -> mol day-1)
            % handle special case of CO2*
            if strcmp(I.ATM_names{n},'pCO2')
                ocn=piston.*loc_A_rho.*ECC.state(Ib,I.CO2);
            else
                ocn=piston.*loc_A_rho.*TRACERS(Ib,ind_ocn);
            end
    
            atm=piston.*loc_A_rho.*Solubility(:,n).*ATM(n);
    
            % net flux (mol day-1), +ve gives net sea to air flux
            Flux=(ocn-atm);
            dCdt(Ib,ind_ocn) = dCdt(Ib,ind_ocn) - (Flux.*ocn_pars.rM(Ib));
    
            % update atmosphere tracer
            dATMdt(n) = dATMdt(n) + sum(Flux)/ocn_pars.atm_mol;
    
        end
    end

end
%%
function [Sc,Solubility]=calc_gasexchange_constants( parameters )
        
        T=parameters.ocn_pars.T(parameters.ocn_pars.Ib,:);
        S=parameters.ocn_pars.S(parameters.ocn_pars.Ib,:);
        rho=parameters.ocn_pars.rho(parameters.ocn_pars.Ib,:);
        
        Sc=zeros(numel(parameters.ocn_pars.Ib),parameters.gen_pars.n_atm);
        Solubility=zeros(numel(parameters.ocn_pars.Ib),parameters.gen_pars.n_atm);

        % Schmidt Numbers
        % N.B. this is an empirical fit between windspeed (m s-1) and piston velocity (cm hr-1)
        T(T<0)=0; 
        T(T>30)=30;
        T2=T.*T;
        T3=T2.*T;
        
        for n=1:parameters.gen_pars.n_atm
            
            Sc(:,n) = parameters.bgc_pars.Sc_constants(n,1)-...
                parameters.bgc_pars.Sc_constants(n,2)*T+...
                parameters.bgc_pars.Sc_constants(n,3)*T2-...
                parameters.bgc_pars.Sc_constants(n,4)*T3;
            
        end
        
        % precalculate
        Sc=(Sc*1.515e-3).^-0.5; % (Sc/660)^1/2
        
        
        % Solubility
        T(T<2)=2; 
        T(T>35)=35; 
        T=T+273.15;
        S(S<26)=26; 
        S(S>43)=43;
        Tr100=T./100;
        
        for n=1:parameters.gen_pars.n_atm
            
            Solubility(:,n) = exp(...
                parameters.bgc_pars.sol_constants(n,1)...
                +parameters.bgc_pars.sol_constants(n,2)*(100./T)...
                +parameters.bgc_pars.sol_constants(n,3)*log(Tr100)...
                +S...
                .*(parameters.bgc_pars.sol_constants(n,4)...
                +parameters.bgc_pars.sol_constants(n,5)*Tr100...
                +parameters.bgc_pars.sol_constants(n,6)*Tr100.*Tr100));
            
            % if not CO2
            if ~strcmp(parameters.ind_pars.ATM_names{n},'pCO2')
                Solubility(:,n)=Solubility(:,n)./(rho.*0.022414); % 0.022414 is Standard Molar Volume
            end
            
        end
        
            
end

%%
function [ dATMdt ] = atm_forcings ( t , ATM , dATMdt , parameters , forcings )
    
    for n=1:size(forcings.atm.meta,1)
        
        % restoring (mol mol-1)
        if forcings.atm.meta(n,1)
            sig=forcings.atm.interp{n}(t);
            dATMdt(:,n) = dATMdt(:,n) + (sig .* forcings.atm.restore_scale(1,n) - ATM(1,n)) * forcings.atm.restore_timescale(1,n);   
        end
        % forcing (mol mol-1 day-1)
        if forcings.atm.meta(n,2)
            sig=forcings.atm.interp{n}(t);
            dATMdt(:,n) = dATMdt(:,n) + (sig .* forcings.atm.data(:,n) * forcings.atm.force_scale(1,n) ) ./ parameters.ocn_pars.atm_mol ./ parameters.gen_pars.conv_d_yr;
        end
    end
 
end
    


function [ dCdt ] = ocn_forcings ( t , TRACERS , dCdt , parameters , forcings)
    
    for n=1:size(forcings.ocn.meta,1)
        
        % restoring (mol kg-1)
        if forcings.ocn.meta(n,1)
            sig=forcings.ocn.interp{n}(t);
            %dCdt(:,n) = dCdt(:,n) + ((sig .* forcings.ocn.data(:,n) * forcings.ocn.restore_scale(1,n)) - TRACERS(:,n)) * parameters.bgc_pars.restore_timescale;
            %dCdt(:,n) = dCdt(:,n) + (((sig .* sum(forcings.ocn.restore_scale(1,n)*sum(parameters.ocn_pars.M))) - sum(TRACERS(:,n).*parameters.ocn_pars.M)) * parameters.bgc_pars.restore_timescale)*forcings.ocn.data(:,n)./parameters.ocn_pars.M;  
            dCdt(:,n) = dCdt(:,n) + ((sig .* forcings.ocn.restore_scale(1,n)) - TRACERS(:,n)) * forcings.ocn.restore_timescale(1,n);

        end
        % forcing (mol kg-1 day-1)
        if forcings.ocn.meta(n,2)
            sig=forcings.ocn.interp{n}(t);
            dCdt(:,n) = dCdt(:,n)+(sig .* forcings.ocn.data(:,n) * forcings.ocn.force_scale(1,n) ) .* parameters.ocn_pars.rM ./ parameters.gen_pars.conv_d_yr;
        end
    end
      
end

function [ PARTICLES ] = sed_forcings ( t , PARTICLES , dCdt , parameters , forcings)
    
    for n=1:size(forcings.sed.meta,1)
        
        % restoring (mol kg-1)
        if forcings.sed.meta(n,1)
            sig=forcings.sed.interp{n}(t);
            %dCdt(:,n) = dCdt(:,n) + ((sig .* forcings.ocn.data(:,n) * forcings.ocn.restore_scale(1,n)) - TRACERS(:,n)) * parameters.bgc_pars.restore_timescale;
            %dCdt(:,n) = dCdt(:,n) + (((sig .* sum(forcings.ocn.restore_scale(1,n)*sum(parameters.ocn_pars.M))) - sum(TRACERS(:,n).*parameters.ocn_pars.M)) * parameters.bgc_pars.restore_timescale)*forcings.ocn.data(:,n)./parameters.ocn_pars.M;  
            PARTICLES(:,n) = PARTICLES(:,n) + ((sig .* forcings.sed.restore_scale(1,n)) - PARTICLES(:,n)) * forcings.sed.restore_timescale(1,n);

        end
        % forcing (mol kg-1 day-1)
        if forcings.sed.meta(n,2)
            sig=forcings.sed.interp{n}(t);
            PARTICLES(:,n) = dCdt(:,n)+(sig .* forcings.sed.data(:,n) * forcings.sed.force_scale(1,n) ) .* parameters.ocn_pars.rM ./ parameters.gen_pars.conv_d_yr;
        end
    end
      
end


function [dCdt] = restore_ocnatm_Cinv ( dCdt , TRACERS , ATM , parameters)
    
    if parameters.bgc_pars.restore_ocnatm_Cinv

        ocnatm_Cinv=sum(TRACERS(:,parameters.ind_pars.DIC).*parameters.ocn_pars.M) + ATM(parameters.ind_pars.pCO2)*parameters.ocn_pars.atm_mol; 

        grid_weighting=parameters.ocn_pars.M./sum(parameters.ocn_pars.M);

        % add a small flux of DIC across the ocean weighted by grid size
        dCdt(:,parameters.ind_pars.DIC)=dCdt(:,parameters.ind_pars.DIC) + (parameters.bgc_pars.ocnatm_Cinv_target-ocnatm_Cinv) .*grid_weighting ./5300 ./parameters.ocn_pars.M;

    end

end

function [Sc_constants , sol_constants] = initialise_gasex_constants ( gen_pars , ind_pars )

    % indices of coefficients
    coeff_index={'pCO2','pO2'};

    % Coefficients
    Sc_coeffs=[2073.1 125.62 3.6276 0.043219 ; ... % CO2
               1953.4 128.00 3.9918 0.050091];     % O2

    sol_coeffs=[-60.2409 93.4517 23.3585 0.023517 -0.023656 0.0047036 ; ... % CO2
                -58.3877 85.8079 23.8439 -0.034892 0.015568 -0.0019387];    % O2
            
    % output coefficients needed by OMG
    Sc_constants=zeros(gen_pars.n_atm,4);
    sol_constants=zeros(gen_pars.n_atm,6);
    
    for n=1:gen_pars.n_atm
        for nn=1:numel(coeff_index)
            if strcmp(ind_pars.ATM_names{n},coeff_index{nn})
                Sc_constants(n,:)=Sc_coeffs(nn,:);
                sol_constants(n,:)=sol_coeffs(nn,:);
            end
        end
    end
                
            
end

function [ dCdt ] = aeolian_Fe ( dCdt , PARTICLES , parameters )

    if parameters.bgc_pars.Fe_cycle

    I=parameters.ind_pars;

    % get total Fe input from aeolian deposition of dust
    % using mass fraction, so convert det from mol kg-1 
    dust_Fe = parameters.bgc_pars.conv_Fe_g_mol .* ...
       parameters.bgc_pars.det_Fe_frac .* ...
       parameters.bgc_pars.conv_det_mol_g .* ...
       PARTICLES(:,I.Det).*parameters.ocn_pars.M;

    % solubility
    det_flux=PARTICLES(:,I.Det).*parameters.gen_pars.conv_d_yr; % mol year-1
    Fe_sol = det_flux.^(parameters.bgc_pars.det_Fe_sol_exp-1.0);
    Fe_sol(isinf(Fe_sol))=0.0;

    % scale solubility to global value
    Fe_sol_scale = sum(Fe_sol.*det_flux) ./ sum(det_flux);
    Fe_sol    = Fe_sol .* (parameters.bgc_pars.det_Fe_sol ./ Fe_sol_scale);

    % apply solubility (mol kg-1 day-1)
    dCdt(:,I.TDFe) = dCdt(:,I.TDFe) + (Fe_sol .* dust_Fe)./parameters.ocn_pars.M;

    end

end

function [ Fe , FeL , L ] = Fe_speciation ( TDFe , TL , K_FeL )

    p = [ 1 , -(TL+TDFe+1/K_FeL) , TDFe*TL ];
    r = roots(p);

    if min(r)<1e-15
        FeL=max(r);
    else
        FeL=min(r);
    end

    Fe = TDFe - FeL;
    L  = TL - FeL;

end


function [ scav_Fe ] = scavenge_Fe ( Fe , parameters , POC )

    % convert particle concentration from mol kg-1 to mg l-1
    Cp = POC * 1e6 * 12.01 * (1000/1027.649);

    % scavenging rate (d-1)
    scav_Fe_k = parameters.bgc_pars.scav_Fe_sf_POC ...
        * parameters.bgc_pars.scav_Fe_k0 ...
        (Cp.^parameters.bgc_pars.scav_Fe_exp);
    
    % scavenged Fe (mol kg-1 d-1)
    scav_Fe = scav_Fe_k .* Fe;

    % cap at maximum available Fe
    ind=scav_Fe>Fe;
    scav_Fe(ind)=Fe(ind);

end


function [ Fe_to_P ] = calc_variable_Fe_to_P ( TRACERS , parameters )

    % note, C:FE max is pre-set

    Ib=parameters.ocn_pars.Ib;
    ind=TRACERS(Ib,parameters.ind_pars.TDFe)>0.125e-9;

    Fe_to_C=parameters.bgc_pars.stoichiometry(Ib,parameters.ind_pars.TDFe);
    Fe_to_C(ind) = min([Fe_to_C(ind),...
        (parameters.bgc_pars.FetoC_C + parameters.bgc_pars.FetoC_K .* ...
        (1e9.*TRACERS(Ib,parameters.ind_pars.TDFe) .^ parameters.bgc_pars.FetoC_pP))],[],2);
    Fe_to_C=1./Fe_to_C;

    Fe_to_P=Fe_to_C.*parameters.bgc_pars.C_to_P;

end






