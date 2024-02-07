%% -------- Initialise Ecosystem -------- %%
if eco_pars.use_parallel
    
    eco_pars.n_parallel = feature('numcores');
    if eco_pars.n_parallel==40
        eco_pars.n_parallel = 20;
    end
    
    ii=round(linspace(0,ocn_pars.ni,eco_pars.n_parallel+1)); % divide ocean surface among n proc
    idiv = sort(diff(ii),'descend');       % get size of blocks on each proc
    ii=[0 cumsum(idiv)];                   % sort ii as for idiv
    
    for ipar=1:eco_pars.n_parallel
        % get surface index of spatially divided domains
        indx_surf  = (ii(ipar)+1) : (ii(ipar+1));
        indomain   = [ocn_pars.i(ocn_pars.Ib(indx_surf)) ocn_pars.j(ocn_pars.Ib(indx_surf))];
        % find full depth index of spatially divided domains
        indx_full  = find(ismember([ocn_pars.i ocn_pars.j],indomain,'rows'));
        % find interior index of spatially divided domains
        indx_interior  = find(ismember([ocn_pars.i(ocn_pars.Ii) ocn_pars.j(ocn_pars.Ii)],indomain,'rows'));
        
        % subdivide ocn_pars
        ocn_pars_proc{ipar}            = ocn_pars;
        
        ocn_pars_proc{ipar}.indx_surf      = indx_surf;
        ocn_pars_proc{ipar}.indx_full      = indx_full;
        ocn_pars_proc{ipar}.indx_interior  = indx_interior;
        
        ocn_pars_proc{ipar}.i            = ocn_pars.i(indx_full);
        ocn_pars_proc{ipar}.j            = ocn_pars.j(indx_full);
        ocn_pars_proc{ipar}.k            = ocn_pars.k(indx_full);
        ocn_pars_proc{ipar}.rk           = ocn_pars.rk(indx_full);
        ocn_pars_proc{ipar}.nb           = size(ocn_pars_proc{ipar}.i,1);
        ocn_pars_proc{ipar}.ni           = size(ocn_pars_proc{ipar}.Ib,1);
        ocn_pars_proc{ipar}.T            = ocn_pars.T(indx_full,:);
        ocn_pars_proc{ipar}.S            = ocn_pars.S(indx_full,:);
        ocn_pars_proc{ipar}.depth        = ocn_pars.depth(indx_full,:);
        ocn_pars_proc{ipar}.M            = ocn_pars.M(indx_full,:);
        ocn_pars_proc{ipar}.rM           = ocn_pars.rM(indx_full,:);
%         ocn_pars_proc{ipar}.A            = ocn_pars.A(indx_full,:);
        ocn_pars_proc{ipar}.rho          = ocn_pars.rho(indx_full,:);
        ocn_pars_proc{ipar}.solfor       = ocn_pars.solfor(indx_full,:);
        ocn_pars_proc{ipar}.MLD          = ocn_pars.MLD(indx_full,:);
        ocn_pars_proc{ipar}.wspeed       = ocn_pars.wspeed(indx_full,:);
%         ocn_pars_proc{ipar}.atm_A        = ocn_pars.atm_A(indx_full,:);
        ocn_pars_proc{ipar}.atm_temp     = ocn_pars.atm_temp(indx_full,:);
        ocn_pars_proc{ipar}.PAR0         = ocn_pars.PAR0(indx_full,:);
%         ocn_pars_proc{ipar}.atm_V        = ocn_pars.atm_V(indx_full);
        ocn_pars_proc{ipar}.conv_atm_mol = ocn_pars.conv_atm_mol(indx_full,:);
        ocn_pars_proc{ipar}.conv_mol_atm = ocn_pars.conv_mol_atm(indx_full,:);
        ocn_pars_proc{ipar}.PO4_obs      = ocn_pars.PO4_obs(indx_full);
        
        ocn_pars_proc{ipar}.Ib           = find(ocn_pars_proc{ipar}.k==max(ocn_pars_proc{ipar}.k));
        ocn_pars_proc{ipar}.Ii           = find(ocn_pars_proc{ipar}.k~=max(ocn_pars_proc{ipar}.k));

%         for j=1:length(ocn_pars_proc{ipar}.TMA)
%             [~,jj,~]=find(ocn_pars.TMA{j}(indx_full,:));
%             jj=unique(jj);
%             ocn_pars_proc{ipar}.TMA{j} = ocn_pars.TMA{j}(indx_full,jj);
%             
%             [~,kk,~]=find(ocn_pars.TMB{j}(indx_surf,:));
%             kk=unique(kk);
%             ocn_pars_proc{ipar}.TMB{j} = ocn_pars.TMB{j}(indx_surf,kk);
%         end
%         ocn_pars_proc{ipar}.upstream_full = jj;
%         ocn_pars_proc{ipar}.upstream_surf = kk;
        
        % subdivide bgc_pars
        bgc_pars_proc{ipar}            = bgc_pars;
        bgc_pars_proc{ipar}.POM_matrix = bgc_pars.POM_matrix(indx_full,indx_full);
        
        % subdivide parameters
        parameters_proc{ipar}.gen_pars = gen_pars;
        parameters_proc{ipar}.bgc_pars = bgc_pars_proc{ipar};
        parameters_proc{ipar}.eco_pars = eco_pars;
        parameters_proc{ipar}.ocn_pars = ocn_pars_proc{ipar};
        parameters_proc{ipar}.ind_pars = I;
        
    end
    
    
    % Initialise parallel pool
    if isempty(p)
        parpool(eco_pars.n_parallel);
    elseif p.NumWorkers~=eco_pars.n_parallel
        delete(gcp);
        parpool(eco_pars.n_parallel);
    end
    
else
    if ~isempty(p)
        delete(gcp);
    end
    eco_pars.n_parallel = 1;
end

% JDW - check for variables hanging around workspace












