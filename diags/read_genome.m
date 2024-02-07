function [binary_gene] = read_genome(genome,biomass)

    % get locations of non-zero biomass
    kk=find(reshape(biomass,[],1)); 

    ngenes=size(genome,2); % number of 64 bit genes
    
    %% BINARY GENOME
    if ngenes>0
        % extract genes for extant populations
        extant=full(genome(kk,:));
        
        npopn=numel(kk);      % number of extant populations
        nbit = 53; % log2(FlIntMax) = 53 --> maximum integer precision of floating point values
        
        % set probabilistic mutation rates
        % each gene in a genome has a different probability of a bit flip
        % first gene is most likely, last gene is 1000 times less likely
        % so first gene will saturate first, while last gene very unlikely to saturate
        ww=fliplr(logspace(0,3,ngenes));
        w=reshape(repmat(ww,nbit,1),1,[]); % weight genes to have different mutation rates
        w=w./sum(w); % normalise weights
        
        igene = randsample(ngenes, npopn, true, ww)'; % identify gene in which to flip single base
        ibase=randi(nbit-1,npopn,1);    % identify base to flip in that gene (excluding left-most bit, which is needed to occupy a space in the sparse genome matrix)
        
        % get index of genes to mutate within overall extant population
        genind=sub2ind([npopn ngenes],1:npopn,igene);
        
        % extract fp_values of mutant genes of extant populations
        fp_values=extant(genind);
        
        % get original value of mutated bases (bits)
        mut=bitget(fp_values,ibase');
        
        % get new genome after flipping a single bit
        extant(genind)=bitset(fp_values,ibase',1-mut);
        
        % set left-most bit of all genes in extant populations to one
        extant=bitset(extant,nbit); % (this is so all populations are assigned a full row in the sparse genome)
        
        % place back in full genome array
        genome(kk,:)=extant;
    end
    
end