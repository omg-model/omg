
bgc_pars.PO4_restore_data_file = gen_fcns.GetFullPath(bgc_pars.PO4_restore_data_file); % Convert relative to absolute path
%% Biogeochemistry Parameters

bgc_pars.rDOP_frac=1-bgc_pars.DOP_frac;                                     % fraction of POM exported    
if bgc_pars.DOP_k>0
    bgc_pars.DOP_k=(1/bgc_pars.DOP_k);                                          % reciprical and year
    bgc_pars.DOP_k=bgc_pars.DOP_k/360;
else
    bgc_pars.DOP_k=0;
end

bgc_pars.K_FeL=10.^bgc_pars.K_FeL_pP; 

% calculate POM remineralisation matrix if needed
if strcmp(bgc_pars.remin_scheme,'matrix')
    [POM_matrix , CaCO3_matrix]=bgc_fcns.create_remin_matrices(bgc_pars,ocn_pars);
end
bgc_pars.POM_matrix   = POM_matrix;
bgc_pars.CaCO3_matrix = CaCO3_matrix;

% uptake stoichiometry - JDW: move this to tracer definitions as defines
% relationships between tracers?
bgc_pars.uptake_stoichiometry = zeros(ocn_pars.nb,gen_pars.n_bgc_tracers);

bgc_pars.uptake_stoichiometry(:,I.PO4)=1.0;
bgc_pars.uptake_stoichiometry(:,I.DOP)=0.0;
if(bgc_pars.CARBCHEM_select)
    bgc_pars.uptake_stoichiometry(:,I.DIC)=bgc_pars.C_to_P;
    bgc_pars.uptake_stoichiometry(:,I.ALK)=bgc_pars.ALK_to_P;
    %bgc_pars.uptake_stoichiometry(:,I.DOC)=-bgc_pars.C_to_P;
end
if(bgc_pars.O2_select)
    bgc_pars.uptake_stoichiometry(:,I.O2)=bgc_pars.O_to_P;
end
if(bgc_pars.Fe_cycle)
    bgc_pars.uptake_stoichiometry(:,I.TDFe)=bgc_pars.Fe_to_C; % n.b. will be converted to Fe_to_P
    bgc_pars.uptake_stoichiometry(:,I.TL)=0;
end


% remin stoichiometry - JDW: move this to tracer definitions as defines
% relationships between tracers?
bgc_pars.remin_stoichiometry = zeros(ocn_pars.nb,gen_pars.n_bgc_tracers);

bgc_pars.remin_stoichiometry(:,I.PO4)=0.0;
bgc_pars.remin_stoichiometry(:,I.DOP)=0.0;
if(bgc_pars.CARBCHEM_select)
    bgc_pars.remin_stoichiometry(:,I.DIC)=0.0;
    bgc_pars.remin_stoichiometry(:,I.ALK)=bgc_pars.ALK_to_P;
    bgc_pars.remin_stoichiometry(:,I.DOC)=0.0;
end
if(bgc_pars.O2_select)
    bgc_pars.remin_stoichiometry(:,I.O2)=bgc_pars.O_to_P;
end
if(bgc_pars.Fe_cycle)
    bgc_pars.remin_stoichiometry(:,I.TDFe)=0.0;
    bgc_pars.remin_stoichiometry(:,I.TL)=0.0;
end

% map OCN to DOM and Particles
row=find(~cellfun(@isempty,I.OCN_to_DOM));
col=0; count=1;
for n=1:numel(I.OCN_to_DOM)
    tmp=strfind(I.OCN_names,I.OCN_to_DOM{n});
        if ~isempty(find(~cellfun(@isempty,tmp)))
        col(count)=find(~cellfun(@isempty,tmp));
        count=count+1;
    end
end    
bgc_pars.mapOCN_DOM=sparse(row,col,ones(size(row)),gen_pars.n_bgc_tracers,gen_pars.n_bgc_tracers);


row=find(~cellfun(@isempty,I.OCN_to_POM));
col=0; count=1;
for n=1:numel(I.OCN_to_POM)
    tmp=strfind(I.SED_names,I.OCN_to_POM{n});
        if ~isempty(find(~cellfun(@isempty,tmp)))
        col(count)=find(~cellfun(@isempty,tmp));
        count=count+1;
    end
end    
bgc_pars.mapOCN_POM=sparse(row,col,ones(size(row)),gen_pars.n_particles,gen_pars.n_particles);

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

% make water column cumulative sum matrix
count=1; 
tmp=0;
for n=1:934
    tmp=tmp+numel(find(ocn_pars.wc==n))^2;
end

spcol_ind=zeros(tmp,1);
sprow_ind=zeros(tmp,1);
spval=zeros(tmp,1);


for wc=1:ocn_pars.nb
    n_wc=numel(find(ocn_pars.wc==wc));
    wc_start=find(ocn_pars.wc==wc,1,'first');
    wc_end=find(ocn_pars.wc==wc,1,'last');
    %tmp_row=reshape(triu(ones(n_wc,n_wc)),[],1);
    tmp_row=reshape(triu(repmat(find(ocn_pars.wc==wc),1,n_wc)),[],1);
    tmp_col=reshape(repmat([wc_start:1:wc_end],n_wc,1),[],1);
    spcol_ind(count:count+numel(tmp_col)-1,1)=tmp_col;
    sprow_ind(count:count+numel(tmp_row)-1,1)=tmp_row;
    count=count+numel(tmp_row);
end

ind=sprow_ind==0;
sprow_ind(ind)=[];
spcol_ind(ind)=[];

ocn_pars.wc_cumsum=sparse(spcol_ind,sprow_ind,ones(size(sprow_ind,1),1),ocn_pars.nb,ocn_pars.nb);

% make water column sum matrix
ocn_pars.wc_sum=sparse(ocn_pars.wc,[1:1:ocn_pars.nb],ones(ocn_pars.nb,1),ocn_pars.nb,ocn_pars.nb);


%% Load forcing data

forcings=struct();
forcings=gen_fcns.load_forcing_data(gen_pars,bgc_pars,forcings,I,ocn_pars);


%%
