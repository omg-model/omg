
%% -------- Initialise Ecosystem -------- %%
if strcmp(bgc_pars.uptake_scheme,'eco')
    

    

        
    %% Size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Need to make sure size classes 'snap' to optimal pp interval
    n_OoM=min(eco_pars.nsize-1,4); % Minimum number of orders of magnitude to span in length dimension
    if eco_pars.nsize==1
        minsize=6; % only size class
    else
        minsize=0.6; % size of smallest class
    end
    l_ratio=log10(nthroot(eco_pars.pp_opt,3)); % find log10 of optimal length ratio (as cube root of volume ratio)
    proceed=1; % n size classes in optimal pp interval
    d=1;
    while proceed % If smaller spacing can still achieve required length span...
        f=l_ratio/d; % calculate log spacing
        if (l_ratio/(d+1)) >= (n_OoM/(eco_pars.nsize-1))
            d=d+1; % add one more size class in optimal pp interval
        else
            proceed=0;
        end % stop iteration if spacing too close
    end
    ESD = minsize.*10.^((0:(eco_pars.nsize-1)).*f)'; % calculate ESDs of all plankton classes

    %     ESD  = logspace(log10(0.6),log10(6000),eco_pars.nsize)'; % equivalent spherical diameter
    ESR  = 0.5.*ESD;     % equivalent spherical radius
    Vol  = 4/3*pi.*ESR.^3; % cell volume
    eco_pars.unq_ESD = ESD;
    
    %% thermal optima %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % define midpoints of thermal tolerance curves
    Topt = linspace(min(OCEAN.T(:)),max(OCEAN.T(:)),eco_pars.nTopt);
    % define boundaries between thermal tolerance ranges
    Topt_int = [-inf conv(Topt, [0.5 0.5], 'valid') inf]';

    % Transpose Topt
    Topt    = Topt';
    % get high and low boundaries realted to each midpoint
    Topt_lo = Topt_int(1:end-1);
    Topt_hi = Topt_int(2:end);
    
    %% trophic strategy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    trophic = fliplr(linspace(0,1,eco_pars.ntroph))';
    
    nrgb     = eco_pars.nrgb;
    ngenes   = eco_pars.ngenes;

    %% grid traits and reshape to one dimension
    
    [iTopt,itrophic,iVol] = ndgrid(1:eco_pars.nTopt,1:eco_pars.ntroph,1:eco_pars.nsize);
    
    eco_pars.Topt      = Topt(iTopt(:));
    eco_pars.Topt_lo   = Topt_lo(iTopt(:));
    eco_pars.Topt_hi   = Topt_hi(iTopt(:));
    
    eco_pars.trophic   = trophic(itrophic(:));

    eco_pars.V         = Vol(iVol(:));
    eco_pars.ESD       = nthroot(6.*eco_pars.V./pi,3);
    
    %% Trait dimension indices

    ndims  = [eco_pars.nTopt eco_pars.ntroph eco_pars.nsize];
    itraits = find(ndims>1);
    % Two (and only two) trait dimensions are be evaluated at one time
    if size(itraits,2)>2
        error('Cannot have variability in more than two trait dimensions')
    elseif size(itraits,2)==1
        itraits = [1 max(itraits,2)];
    elseif isempty(itraits)
        itraits = 1:2;
    end

    tnames    = {'Thermal','Trophic','Size'};
    traits    = {unique(eco_pars.Topt),unique(eco_pars.trophic),eco_pars.unq_ESD};
    alltraits = {eco_pars.Topt,eco_pars.trophic,log10(eco_pars.ESD)};

    eco_pars.nT1      =  ndims(itraits(1));
    eco_pars.nT2      =  ndims(itraits(2));
    eco_pars.T1       = traits{itraits(1)};
    eco_pars.T2       = traits{itraits(2)};
    eco_pars.Tname1   = tnames{itraits(1)};
    eco_pars.Tname2   = tnames{itraits(2)};
    
    %% Trait diffusion/mutation matrix

    [xn,yn]=ndgrid(1:eco_pars.nT1,1:eco_pars.nT2);
    xl=reshape(xn,numel(xn),1);
    yl=reshape(yn,numel(yn),1);
    tvect=[xl yl];
    
    nn=eco_pars.jpmax;
    tm=zeros(nn);

    % if 'Everything' selected for P diversity
    if isempty(eco_pars.seed_PHY)
        % make mutation matrix that produces mutants of *all* phenotypes
        tm=ones(nn).*eco_pars.mutrat(1);
        tm(eye(nn)==1)=1-sum(tm,2);
    else
    
        T1 = alltraits{itraits(1)};
        T2 = alltraits{itraits(2)};

        traitdiff = [(max(OCEAN.T(:))-min(OCEAN.T(:)))./(eco_pars.nTopt-1) ...
                     1/(eco_pars.ntroph-1) ...
                     (eco_pars.nsize-1).*f/(eco_pars.nsize-1)];
        traitdiff(~isfinite(traitdiff))=1;
        deltaT1 = traitdiff(itraits(1)); 
        deltaT2 = traitdiff(itraits(2));
       
        % differences between trait values in discrete grid
        p1 = range(T1)^2 .* eco_pars.mutrat(1)./deltaT1.^2; % probabilty of mutation in T1
        p2 = range(T2)^2 .* eco_pars.mutrat(2)./deltaT2.^2; % probabilty of mutation in T2
    
        for i=1:nn
            for j=1:nn
                adjvect=abs(tvect(i,:)-tvect(j,:));
                % populations adjacent in first dimension only
                if adjvect(1)==1 && adjvect(2)==0
                    % normalise flux to range of trait space
                    tm(i,j) = p1-2*p1*p2;
                end
                % populations adjacent in second dimension only
                if adjvect(1)==0 && adjvect(2)==1
                    % normalise flux to range of trait space
                    tm(i,j) = p2-2*p1*p2;
                end
                % populations adjacent in both dimensions
                if adjvect(1)==1 && adjvect(2)==1
                    % normalise flux to range of trait space
                    tm(i,j) = p1*p2;
                end
                % clonal reproduction
                if i==j
                    % normalise flux to range of trait space
                    tm(i,j) = (1-2*p1)*(1-2*p2);
                end
            end
        end
    end

    eco_pars.Pmut=sparse(tm);
    
    % adjust threshold by number of bins 
        % more bins will have less biomass in each, producing less flux
        % decreasing threshold makes it easier 
            % for smaller populations to yield successful mutants
    eco_pars.functional_extinction = eco_pars.functional_extinction ./ eco_pars.jpmax;
    % set numerical threshold a million times lower than functional
%     eco_pars.numerical_extinction = eco_pars.functional_extinction .* 1e-6; % NOT CURRENTLY USED
    
    clear trophic Vol xl yl xn yn ESD ESR
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    eco_pars.gamma   = [];
    
    eco_pars.muinf      = eco_pars.muinf_a .* eco_pars.V .^ eco_pars.muinf_b;
    eco_pars.Qmin       = eco_pars.Qmin_a  .* eco_pars.V .^ eco_pars.Qmin_b;
    eco_pars.Vmax       = eco_pars.Vmax_a  .* eco_pars.V .^ eco_pars.Vmax_b;
    eco_pars.alpha      = eco_pars.alpha_a      .* eco_pars.V .^ eco_pars.alpha_b;
    eco_pars.clearance  = eco_pars.clearance_a  .* eco_pars.V .^ eco_pars.clearance_b;
    eco_pars.gmax       = eco_pars.gmax_a       .* eco_pars.V .^ eco_pars.gmax_b;
    
    eco_pars.allometric.muinf       =[eco_pars.muinf_a      eco_pars.muinf_b ];
    eco_pars.allometric.Qmin        =[eco_pars.Qmin_a       eco_pars.Qmin_b  ];
    eco_pars.allometric.Vmax        =[eco_pars.Vmax_a       eco_pars.Vmax_b  ];
    eco_pars.allometric.alpha       =[eco_pars.alpha_a      eco_pars.alpha_b ];
    eco_pars.allometric.clearance   =[eco_pars.clearance_a  eco_pars.clearance_b];
    eco_pars.allometric.gmax        =[eco_pars.gmax_a       eco_pars.gmax_b  ];
    eco_pars=rmfield(eco_pars,{'muinf_a','muinf_b','Qmin_a','Qmin_b','Vmax_a','Vmax_b','alpha_a','alpha_b','clearance_a','clearance_b','gmax_a','gmax_b'});
    
    eco_pars.mumax=eco_pars.muinf .* eco_pars.Vmax ...
        ./ (eco_pars.muinf .* eco_pars.Qmin + eco_pars.Vmax ); % d^{-1}
    
    % Grazing matrix
    [Vpry,Vprd]= meshgrid(eco_pars.V,eco_pars.V);
    prd_pry=Vprd./Vpry;
    herbmat=exp(-(log(prd_pry./eco_pars.pp_opt)).^2 ./ (2.*eco_pars.pp_sig.^2));
    herbmat(herbmat<1e-3)=0;
    eco_pars.prdpry=herbmat;
    eco_pars.pryprd=eco_pars.prdpry';
    
    clear herbmat Vpry Vprd prd_pry
    
    if eco_pars.ntroph>1
        % apply trophic trade-offs
        eco_pars.mumax     = eco_pars.mumax     .*    eco_pars.trophic .^eco_pars.tau;
        eco_pars.gmax      = eco_pars.gmax      .* (1-eco_pars.trophic).^eco_pars.tau;
        eco_pars.alpha     = eco_pars.alpha     .*    eco_pars.trophic .^eco_pars.tau;
        eco_pars.clearance = eco_pars.clearance .* (1-eco_pars.trophic).^eco_pars.tau;
    end
    
    
    
    if eco_pars.trscale~=1         % scale plankton transport terms
        Ptrans=ones(eco_pars.jpmax,1).*eco_pars.trscale;
        eco_pars.trscale=spdiags([1;1;Ptrans],0,gen_pars.n_tracers,gen_pars.n_tracers);
    end
    
    % make matrix for interference competition
    eco_pars = interference(eco_pars);
    
    % set lowest rate of binary gene mutation based on length of simulation
    low_prob = floor(log10(1./(runtime.*gen_pars.n_dt)));
    eco_pars.pmut = logspace(0,low_prob,eco_pars.ngenes); % weight genes to have different mutation rates



    
end





%%
function eco_pars = interference(eco_pars)
    
    % square of trait differences
	delta_x2 = (log10(eco_pars.V)-log10(eco_pars.V)').^2;
	delta_y2 = (eco_pars.trophic -eco_pars.trophic' ).^2;
    
    % variance of intererference kernel
    sigma_x2 = eco_pars.sigma_interference_size.^2;
    sigma_y2 = eco_pars.sigma_interference_troph.^2;
    
    % combined uncorrelated distance
    z = delta_x2 ./ sigma_x2 + delta_y2 ./ sigma_y2;
    
    % interference_matrix
    intmat = exp(-z./2);
    intmat(intmat<1e-6) = 0;
    eco_pars.interference_matrix = sparse(intmat);
    
end





