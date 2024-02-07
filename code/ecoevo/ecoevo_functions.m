function [ eco_fcns ] = ecoevo_functions ( )
    eco_fcns.environmental_limits   = @environmental_limits;
    eco_fcns.ecosystem              = @ecosystem;
    eco_fcns.reaper                 = @reaper;
    eco_fcns.mol_clock              = @mol_clock;
    eco_fcns.gene_mutate            = @gene_mutate;
    eco_fcns.maxrows                = @maxrows;
    eco_fcns.maxcols                = @maxcols;
%     eco_fcns.MaxJacobian            = @MaxJacobian;
end

%%
% ----- SUBROUTINES ----- %


%%
function [gamma_I gamma_T]=environmental_limits(P_Biomass,parameters);
    ocn_pars = parameters.ocn_pars;
    eco_pars = parameters.eco_pars;
    I        = parameters.ind_pars;
    
    Ib        = ocn_pars.Ib';
    SLD       = ocn_pars.z0; 
    
    PAR0 = ocn_pars.PAR0(Ib)'; % PAR in W/m^2
    TEMP = ocn_pars.T(Ib)';
    
    ktot = eco_pars.kw + eco_pars.kchl .* sum(P_Biomass,1).*1.59;
    Irr  = PAR0 ./ ktot ./ SLD .* (1-exp(-ktot.*SLD)); % Mean PAR in surface layer

    gamma_I  = Irr./sqrt(Irr.^2 + eco_pars.kPAR^2);  % light-limitation term

    gamma_T  = exp(eco_pars.Tslope*(TEMP-eco_pars.Tref'));

    if eco_pars.nTopt > 1

        %% make thermal optima a range between midpoints of Topt vector

% THIS IS MIGHT BE QUITE COMPUTATIONALLY EXPENSIVE???
    % Seems on a par with other parts of the function
       
        % set default distance matrix to zero
        D=zeros(length(eco_pars.Topt),length(TEMP));

        % distances between ambient temp and shifted phenotypes
        D_upper = TEMP-eco_pars.Topt_hi;
        D_lower = TEMP-eco_pars.Topt_lo;

        % find where they apply ...
        i_upper = (D_upper>0);
        i_lower = (D_lower<0);

        % ... and put them there
        D(i_upper)=D_upper(i_upper);
        D(i_lower)=D_lower(i_lower);
   
        gamma_T = gamma_T .* (1-(D./(eco_pars.Twidth/2)).^2);
        % gamma_T = gamma_T .* exp(-(D./eco_pars.Twidth).^2);
        
    end
    % set < zero to zero
    gamma_T(gamma_T<0)=0;
end


%%
function [dCdt dPOPdt invfit ggr_out]=ecosystem(SURFACE,parameters);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    gen_pars=parameters.gen_pars;
    bgc_pars=parameters.bgc_pars;
    ocn_pars=parameters.ocn_pars;
    eco_pars=parameters.eco_pars;

    I=parameters.ind_pars;
    I.PHY = gen_pars.n_bgc_tracers+1:gen_pars.n_tracers;

    SURFACE = SURFACE.*1e6;

    % Extract state variables 
    % (using non-negative values to calculate rates of change)
    PO4   = max(SURFACE(:,I.PO4),0)'; % [1  x ns]
    D     = max(SURFACE(:,I.DOP),0)'; % [1  x ns]
    P     = max(SURFACE(:,I.PHY)-eco_pars.functional_extinction,0)'; % [jp x ns]

    clear SURFACE

    mumax  = eco_pars.mumax;
    gmax   = eco_pars.gmax;
    v_ones = ones(1,ocn_pars.ni); % vector for padding search_rate array
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Light- and temperature-dependent limiting factors (gamma_I and gamma_T)
    [gamma_I, gamma_T] = environmental_limits(P,parameters);
    
    % NUTRIENT UPTAKE (i.e. AUTOTROPHIC GROWTH)
    alphaPO4   = eco_pars.alpha .* PO4;   % Implicit Expansion of row and column - expands input vectors for matrix of all possible combinations [jp x ns]
    PO4limit   = gamma_I .* alphaPO4 ./ (mumax + alphaPO4);
    PO4uptk    = gamma_T .* mumax .* PO4limit;
    PO4uptk(mumax + alphaPO4<=0) = 0; % Check for possible divide by zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GRAZING
    % total prey available to each predator
    % [mmol Prey] per m^3
    pry_P       = eco_pars.prdpry * P; % [jpred x jprey] *[jprey x ns] = [jpred x ns]
    
    % amount of water that is searched (if no food) per pred per unit time
    % (includes "refuge")
    % m^3 per [mmol Pred] per day
%     search_rate    = eco_pars.clearance .* spfun(@(x) 1-exp(eco_pars.refuge.*x),pry_P); 
    if isinf(eco_pars.refuge)
        search_rate = eco_pars.clearance*v_ones;
    else
        search_rate = eco_pars.clearance .* (1-exp(eco_pars.refuge.*pry_P)); 
    end
    % [jpred x ns] .* [1] = [jpred x ns]
    
    % biomass of all prey encounters for each predator
    % m^3 per [mmol Pred] per d * [mmol Prey] per m^3 
    % = [mmol Prey] per [mmol Pred] per day
    encounter_rate  = search_rate .* pry_P;

    % [jpred x ns] .* [jprey x ns] = [jp x ns]
    
    % amount of water that can be effectively cleared per pred per unit time (accounting for satiation) 
    % [mmol Prey] per [mmol Pred] per day * m^3 per [mmol Pred] per day / [mmol Prey] per [mmol Pred] per day 
    % = m^3 per [mmol Pred] per day
    grazing_limit = encounter_rate ./ (gmax + encounter_rate);
    grazing_rate = gmax  .* gamma_T .* grazing_limit;
    % [jpred x ns   ].*[[jpred x ns]+[jprey x ns]] = [jp x ns]
    
    % Check for possible divide by zero
    grazing_rate(eco_pars.gmax + encounter_rate<=0) = 0;                   
    
% GAINS AS PREDATOR (specific to predator biomass)
    % m^3 per [mmol Pred] per day * [mmol Prey] per m^3 
    % = [mmol Prey] per [mmol Pred] per day
    GGain       = grazing_rate .* pry_P; 
    % [jp x ns   ].*[jprey x ns] = [jp x ns]  
    
% LOSSES AS PREY (specific to predator biomass)
    % m^3 per [mmol Pred] per day * [mmol Pred] per m^3 
    % = per day
    GLoss       = eco_pars.pryprd * (grazing_rate .* P); 
    % [jprey x jpred] *[jpred x ns] = [jp x ns] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GGain and GLoss should sum to equal value
    if nnz(GLoss)~=0 && abs((sum(reshape(GGain.*P,1,[]))-sum(reshape(GLoss.*P,1,[])))/sum(reshape(GLoss.*P,1,[])))>1e-10
        warning('Grazing gains and losses not balanced!!!!')
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LINEAR MORTALITY
    BaseMort = eco_pars.linearmort;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INTERFERENCE
    if eco_pars.interference_strength>0
        Interference = eco_pars.interference_strength .* eco_pars.interference_matrix*P;
    else
        Interference = 0;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GROSS GROWTH RATE
    ggr = PO4uptk + GGain.*eco_pars.lambda;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % REPRODUCTIVE FLUX (including 'mutations')
    mutflux = eco_pars.Pmut * (ggr.*P); % GGRs subject to mutation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ODEs
    dPO4dt = - sum(PO4uptk.*P,1);

    dDdt = + sum((BaseMort + Interference + GGain .* (1-eco_pars.lambda)).*P ,1);  % Mortality + sloppy feeding;
          
    dPdt = + mutflux   ...           % ggr including mutational flux
         - ( BaseMort + Interference + GLoss).*P;    % linear mortality + grazing losses

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    invfit = (+ ggr   ...              % gross growth rate
              - BaseMort ...           % linear mortality
              - GLoss);                % Grazing losses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Output
    dCdt = [dPO4dt;dDdt.*bgc_pars.DOP_frac;dPdt]'./1e6;  
    dPOPdt = (dDdt.*bgc_pars.rDOP_frac)'./1e6;

    ggr_out = ggr';
end

%%
function y = reaper(y,threshold,I)
    i_extinct = y(:,I.PHY)<threshold; % Index all plankon     below extinction threshold
    i_extant  = 1 - i_extinct;        % Index all plankon not below extinction threshold
    
    % pass extincted stuff to POM
    y(:,I.DOP) = y(:,I.DOP) + sum( y(:,I.PHY).*i_extinct ,2 );
    % apply extinction to PHY biomass below threshold
    y(:,I.PHY) =                   y(:,I.PHY).*i_extant;
end


%%
function [genome,rgb] = gene_mutate(genome,rgb,newbio,parameters)

    % get probabilistic mutation rates
    pmut=parameters.eco_pars.pmut;

    % get locations of non-zero biomass
    kk=find(reshape(newbio,[],1)); 

    ngenes=size(genome,2); % number of 64 bit genes
    nrgb =size(rgb   ,2); % number of rgb tracers
    
    %% BINARY GENOME
    if ngenes>0
        % extract genes for extant populations
        extant=full(genome(kk,:));
        
        npopn=numel(kk);      % number of extant populations
        nbit = 53; % log2(FlIntMax) = 53 --> maximum integer precision of floating point values

        % identify genes in which to flip single base (1 in each population)
        igene = find(rand([npopn ngenes])<pmut);

        % identify base to flip in each gene that does mutate
        ibase=randi(nbit,size(igene)); % uniform probability

        % if a mutation occurs
        if numel(igene>0)
            % get original value of mutated bases (bits)
            mut=bitget(extant(igene),ibase);
            % flip relevant bits
            extant(igene)=bitset(extant(igene),ibase,1-mut);
        end

        % place back in full genome array
        genome(kk,:)=extant;

    end
    %% RGB TRACER
    if nrgb>0
        % extract tracers for extant populations
        extant=full(rgb(kk,:));
        % add random mutation
        extant = extant+randn(size(extant));
        % place back in full rgb array
        rgb(kk,:)=extant;
    end
    
end

%%
function [genome,rgb] = mol_clock(genome,rgb,biomass,TM,ggr,Pmut)
    
%% Spatial 'coalescence': Identify maximum (randomised) source to each location (for each population)
    imax=maxrows(TM+eye(size(TM)),biomass.*rand(size(biomass)));
    
    % imax is index of maximum product in intermediate dimension (of matrix multiplication)
    % (i.e. location [matrix row] that dominates supply to each location)
    ii=sub2ind(size(imax),imax,repmat(1:size(imax,2),[size(imax,1),1])); % Convert this to linear index
    
    % reshape index arrays
    ii=reshape(ii,1,[]);
    
    %% Phenotype 'coalescence': Identify maximum source to each population (for each location)
    jmax=maxcols(ggr,Pmut);
    
    % imax is index of maximum product in intermediate dimension (of matrix multiplication)
    % (i.e. population [matrix column] that dominates biomass flux (gross growth and mutation) into each population)
    jj=sub2ind(size(imax),repmat((1:size(jmax,1))',[1,size(jmax,2)]),jmax); % Convert this to linear index
    
    % reshape index arrays
    jj=reshape(jj,1,[]);

    %% Rearrange bioinf arrays
    genome = genome(jj,:); % overcome populations
    rgb    = rgb   (jj,:); % overcome populations
    genome = genome(ii,:); % overcome locations
    rgb    = rgb   (ii,:); % overcome locations
end

%%
function K = maxrows(A,B)
    AT = A.'; % gets the product data contiguous in memory
    K = zeros(size(B));
    zi= 1:size(B,1);
    for j=1:size(B,2) % for each population, find dominant source location at each location
        [v,k] = max(AT.*B(:,j),[],1); % Uses implicit array expansion (all influxes to each location)
        k(v==0) = zi(v==0); % Set index to local point where maximum incoming flux == 0
        K(:,j) = k(:);
    end
end

%%
function K = maxcols(B,A)
    AT = A.'; % gets the product data contiguous in memory
    K = zeros(size(B));
    zi= 1:size(B,2);
    for j=1:size(B,1) % for each location, find dominant source population for each population
        [v,k] = max(B(j,:).*AT,[],2); % Uses implicit array expansion
        k(v==0) = zi(v==0); % Set index to local point where maximum incoming flux == 0
        K(j,:) = k(:);
    end
end
%%

% function [parameters] = MaxJacobian(parameters)
%     % estimate sparsity pattern of Jacobian and add to ode_options
% 
%     gen_pars=parameters.gen_pars;
%     ocn_pars=parameters.ocn_pars;
%     eco_pars=parameters.eco_pars;
% 
% % Maximum direct connectivity of local state variables
%     % Generate empty sparse matrix
%     J_Eco = sparse([],[],[],gen_pars.n_tracers,gen_pars.n_tracers,0);  
%     % Assume bgc variables locally connected
%     J_Eco(1:gen_pars.n_bgc_tracers,1:gen_pars.n_bgc_tracers) = 1; 
%     % Assume plankton see all bgc variables
%     J_Eco(gen_pars.n_bgc_tracers+1:end,1:gen_pars.n_bgc_tracers) = 1; 
%     % Assume plankton affect all bgc variables
%     J_Eco(1:gen_pars.n_bgc_tracers,gen_pars.n_bgc_tracers+1:end) = 1; 
%     % Then add plankton interactions from mutations and grazing
%     J_Eco(gen_pars.n_bgc_tracers+1:end,gen_pars.n_bgc_tracers+1:end) = eco_pars.prdpry | eco_pars.pryprd | eco_pars.Pmut; 
%     
% % Maximum direct connectivity of tracers in physical space
%     J_TM = speye(ocn_pars.nb);
%     for i=1:size(ocn_pars.TMA,1) % collate non-zeros of time-resolved matrices
%         J_TM = (J_TM + (ocn_pars.TMA{i})~=0)~=0;
%     end
%     
% % Tile spatial pattern within local tracer pattern using Kronecker product 
%     J_Pattern = kron(J_Eco,eye(size(J_TM))) | kron(eye(size(J_Eco)),J_TM);
% 
% % Add result to odeoptions
%     parameters.gen_pars.odeoptions=odeset('JPattern',J_Pattern);
% 
%     if false
%         % visualise matrices
%         subplot(311)
%         spy(J_Eco)
%         subplot(312)
%         spy(J_TM)
%         subplot(313)
%         spy(J_Pattern)
%         sparsity_fraction = nnz(J_Pattern)./numel(J_Pattern)
%     end
% end


