function [D] = genome_distances(biomass,genome,parameters,isextant)


    options = optimoptions('lsqcurvefit',...
                           'Display','off',...
                           'UseParallel',false,...
                           'Algorithm','levenberg-marquardt',...
                           'SpecifyObjectiveGradient',true);
    
    % get different mutation rates
    pmut = parameters.eco_pars.pmut; 
    
    % get n genes
    ngenes = parameters.eco_pars.ngenes;

    % number of bits in each gene
    nbit = 53; % log2(flintmax) = 53 --> maximum integer precision of floating point values
    
    % Combine all known parameters into single vector p
    p    = -4.*pmut./nbit;

    % initial guess of t
    t0   = 1;
    
    % t must be greater  than 0 and less than 1/min(pmut) 
    tmin = 0;
    tmax = 1/min(pmut); % ~maximum molecular clock can record in slowest gene
        
    % and genome to matrix (1 column per gene)
    genome_vect  = reshape(genome,numel(biomass),ngenes);

    % extract genes from non-zero locations
    genes  = genome_vect(isextant,:);

    % get number of non zero populations
    npopn = numel(isextant);

    % number of pairwise distances
    n_pairs = nchoosek(npopn,2);

    % convert to unsigned 64 bit integers
    gen64   = uint64(genes);    % [npopulations x ngenes]
    gen64_T = transpose(gen64); % [ngenes x npopulations]

    % convert to binary strings and transpose so ...
        % bits occupy columns
        % columns iterate through genes
        % and then populations
    bingene = transpose(dec2bin(gen64_T,nbit)); % reads through genes then populations

    % reshape so each column contains all bits in one population
    bingene = reshape(bingene,nbit,ngenes,npopn);

    % create empty matrix to fill with Hamming distances
    d_ham=zeros(ngenes,n_pairs);
    % Calculate Hamming distances
    for i_gene=1:ngenes
        d_ham(i_gene,:) = pdist(squeeze(bingene(:,i_gene,:))'-'0','hamming');
    end

    % create empty matrix to fill with estimated divergences
    tstep_est = zeros(1,n_pairs);

    jf=java.text.DecimalFormat; % comma for thousands, three decimal places
    numOut= char(jf.format(n_pairs));

    disp(['Fitting curves to ' numOut ' sets of paired gene Hamming distances (which could take a while...)'])

    WaitMessage = parfor_wait(n_pairs, 'Waitbar', true); 
    tic;
    parfor i_pair=1:n_pairs
        WaitMessage.Send;
        % estimate generations from genome        
        tstep_est(i_pair) = lsqcurvefit(@JukesCantorFun,t0,p',d_ham(:,i_pair),tmin,tmax,options); 
    end
    WaitMessage.Destroy
    t1= seconds(toc);

    t1.Format = 'hh:mm:ss';
    disp(['Time taken = ' char(t1) ' (hh:mm:ss)'])
    
    D=tstep_est./parameters.gen_pars.n_dt;



end