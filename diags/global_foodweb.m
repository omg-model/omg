function global_foodweb(output_dir)

disp(['Processing output from ' output_dir])

% load saved data
% get matObj for matfile data
matObj  =matfile([output_dir '/OUTPUT.mat'])
eco_pars=matObj.eco_pars;
gen_pars=matObj.gen_pars;

load([output_dir '/matrix_vars.mat'])
vi=v_index.i(Ib);
vj=v_index.j(Ib);
vk=v_index.k(Ib);

gen_fcns=general_functions;

%% Prepare figure axes
figure(1)
set(gcf,'units','normalized','outerposition',[0 0 0.5 2/3])
drawnow
gap = [0.01 0.005];
marg_h = 0.030;
marg_w = 0.02;

% color axis limits
maxPO4=0.1;
minPHY=1e-6;
maxPHY=1;

crng=log10([minPHY maxPHY]);
crng(2)=crng(1)+range(crng)*64/62;
caxis(crng);
cmap=colormap(parula(84));
cmap=cmap(21:end,:);
cmap(1,:)=[1 1 1].*0.75;
colormap(cmap);


nphy=eco_pars.nsize*eco_pars.ntroph;
%%
% threshold=eco_pars.extnct;
threshold=0;

PO4=matObj.PO4;
i_timeslice=find(~cellfun(@isempty,PO4))'; % find years with output data

disp(['Processing output from ' num2str(numel(i_timeslice)) ' time-slices'])
if ~exist([output_dir '/Figures/'],'dir')
    wd=cd(output_dir);
    !mkdir Figures
    cd(wd)
end
    
% for nyr=i_timeslice
x=1-eco_pars.trophic;
y=log10(eco_pars.V);
for nyr=i_timeslice
    PHY=cell2mat(matObj.PHY(nyr,1));
    disp(['Processing ' num2str(nyr) ' of '  num2str(numel(i_timeslice)) ', with ' num2str(nnz(PHY)) ' non-zero local populations'])
    
    GENOME=full(cell2mat(matObj.GENOME(nyr,1)));
    nbit=53; % log2(FlIntMax) = 53
    str='';
    for i=1:eco_pars.ngenes
        bin_gene = dec2bin(GENOME(usedata,i),nbit); % convert double value to binary string
        str=[str bin_gene];
    end
    % convert to numeric values (1/0)
    binstr=str-'0';
    % get Hamming distance as weighted number of base differences
    % D = size(binstr,2) * pdist(binstr,'hamming');
    ww=nthroot(2,10).^(0:eco_pars.ngenes-1);
    w=reshape(repmat(ww,nbit,1),1,[]); % weight genes to have different mutation rates
    w=w./mean(w);
    binstr=binstr./w; % multiply bits by inverse mutation weights (less likely worth more)
    D = pdist(binstr,'CityBlock')./size(binstr,2);
    % hackish compensation for slowdown of binary clock
    % this is not publishable
    alpha =      0.01804;
    beta  =       11.44;
    t     = beta.*D ./ (alpha .* (beta - D));
    
    
    
    
    clear vector data array matrix
    
    biomass=sum(PHY,1);
    
    scatter(x,y,biomass.^(2/3).*100+1e-9,'k','filled')
    
    box on
    axis square
    title(gen_pars.save_years(nyr))

    fname=[output_dir '/Figures/Foodweb_' num2str(nyr,'%04i') '.png'];
    export_fig(fname,'-r200')
end
    

% wd=cd([output_dir '/Figures/']);
% !ffmpeg -framerate 14 -i Figures/Frame_%04d.png  -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" evolution.mp4
% cd(wd)
% disp(' ')
% disp('Example ffmpeg command to be executed in ''Figures'' directory ...')
% disp('ffmpeg -framerate 14 -i Frame_%04d.png  -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" evolution.mp4')
% disp(' ')


return
%%
