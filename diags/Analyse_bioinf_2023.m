clear
clc
disp('Analyse_bioinf_2023')

%% Options

% path to model output files
outdir = '~/GitHub/omg/output/EcoEvo_debugging';

% biomass concentration threshold to include population in analysis
threshold = 1e-4;


% flag to plot plankton C biomass
plot_phy=0;
% flag to plot plankton RGB colour trait
plot_rgb=1;

fig1pos = [   1 88 2560 1249];
fig2pos = [2561 88 2560 1249];

% colours of grid lines
lncolr_phy = ones(1,3).*0.9;
lncolr_rgb = ones(1,3).*0.9;

% movie frame rate
fps = 5;

%% load data

load([outdir '/Parameters.mat']);                         % model parameters
PHY0 = ncread([outdir '/omg_fields_netcdf.nc'],'PHY');    % plankton P biomass (moles per kg)
m_time = ncread([outdir '/omg_fields_netcdf.nc'],'time'); % model time (years)

PHY = PHY0.*1e6; % convert from mol/kg to mmol/m^3

% find number of completed years
i_time = find(PHY>0,1,'last');                % get linear index of last non-zero datum
[~,~,~,~,ntime] = ind2sub(size(PHY),i_time);  % get time index from linear index

%% Load genomic data for last year only

start = [1 1 1 1 1 ntime];
count = [Inf Inf Inf Inf Inf 1];

RGB  = ncread([outdir '/omg_fields_netcdf.nc'],'RGB',start,count);    % plankton RGB colour trait
GEN  = ncread([outdir '/omg_fields_netcdf.nc'],'GENOME',start,count); % plankton Binary Genome (written as integers if float


ngenes=size(GEN,5); % number of 64 bit genes
if ngenes==0
    error('No genome data available');
end

%% get dimensions of netcdf data

% trait dimensions
T1     = parameters.eco_pars.T1;
T2     = parameters.eco_pars.T2;
nT1    = parameters.eco_pars.nT1;
nT2    = parameters.eco_pars.nT2;
Tname1 = parameters.eco_pars.Tname1;
Tname2 = parameters.eco_pars.Tname2;

% spatial dimensions
nlat   = size(parameters.ocn_pars.lat,1);
nlon   = size(parameters.ocn_pars.lon,1);

% total number of local populations in metacommunity
npopns = nT1*nT2*nlat*nlon;

% number of rgb traits
nrgb   = parameters.eco_pars.nrgb;

% get meta data
n_sub_tsteps        = parameters.gen_pars.n_sub_tsteps;
save_intra_annual_n = parameters.gen_pars.save_intra_annual_n;

% define points for manually drawn axis grid and tick marks
xticks = linspace(0,nT1*nlon,nT1+1)+0.5;
yticks = linspace(0,nT2 *nlat,nT2+1)+0.5;
xticksmiddle = xticks(2:end)-diff(xticks)./2;
yticksmiddle = yticks(2:end)-diff(yticks)./2;

% get biomass and genomes for last time step
biomass = PHY(:,:,:,:,ntime);
genome  = GEN(:,:,:,:,:);
RGB_vect  = reshape(RGB,numel(biomass),size(RGB,5));


%%

% get index of biomass > threshold
isextant=find(biomass>threshold);

% Calculate pairwise divergence (in years) from Hamming distances
distances = genome_distances(biomass,genome,parameters,isextant);

%%

% Create hierarchical cluster tree
tree = linkage(distances);

% set divergence threshold (years) that define taxa
Cutoff = 100;

% Cluster populations together if less than 'Cutoff' years divergence
clusters = cluster(tree,'Cutoff',Cutoff,'Criterion','distance');

% Find optimal leaf order
leafOrder = optimalleaforder(tree,distances);

% grid traits and locations so that values can be extracted where required
[X_lon,X_lat,X_t1,X_t2] = ndgrid(parameters.ocn_pars.lon,...
    parameters.ocn_pars.lat,...
    parameters.eco_pars.T1,...
    log10(parameters.eco_pars.T2));

% Find cluster IDs
taxa  = clusters(leafOrder)';

% extract locations and phenotypes in order of leaves
X_lat = X_lat(isextant(leafOrder))';
X_lon = X_lon(isextant(leafOrder))';
X_t1  = X_t1(isextant(leafOrder))';
X_t2  = X_t2(isextant(leafOrder))';

% Get RGB gene colours in order of leaves
X_RGB = RGB_vect(isextant(leafOrder),:);

% normalise RGB between 0 and 1
X_RGB=normalize(X_RGB,'range');

% loop through clusters and find mean of all RGB colours within cluster
% these mean colours will be used to differentiate taxa
TCLR=zeros(size(X_RGB));
for i=unique(taxa)
    ii = find(taxa==i);
    clr = mean(X_RGB(ii,:),1);
    for j=1:3
        TCLR(ii,j) = clr(j);
    end
end

% pad 1st dimension of X_RGB and TCLR so they can be shown as images
X_RGB =  permute(repmat(X_RGB,[1 1 2]),[3 1 2]);
TCLR =  permute(repmat(TCLR,[1 1 2]),[3 1 2]);

%%
RGB_phenotype(:,:,1) = repmat(linspace(0,1,nT1),[nT2 1]);
RGB_phenotype(:,:,2) = repmat(linspace(0,1,nT2)',[1 nT2]);
RGB_phenotype(:,:,3) = repmat(linspace(1,0,nT1),[nT2 1]);

clf
imagesc(RGB_phenotype)

%%
% Draw dendrogram
f3 = figure(3);
clf

logtime=false;

set(f3,'color','w');   
ax0 = axes('Position',[0.15 0.300 0.8 0.700]);
ax1 = axes('Position',[0.15 0.275 0.8 0.025]);
ax2 = axes('Position',[0.15 0.250 0.8 0.025]);
ax3 = axes('Position',[0.15 0.225 0.8 0.025]);
ax4 = axes('Position',[0.15 0.200 0.8 0.025]);
ax5 = axes('Position',[0.15 0.175 0.8 0.025]);
ax6 = axes('Position',[0.15 0.150 0.8 0.025]);

axes(ax0)
[H,T,outperm] = dendrogram(tree,0,...
    'Reorder',leafOrder,...
    'ColorThreshold',Cutoff);
set(ax0,'XTick',[]);
box off
if logtime
    yl=ylim;
    set(gca,'YScale','Log')
    ylim([1 yl(2)])
end

% set dendrogram clusters to mean RGB colour within that group
for i=1:size(H,1) % loop through all leaves
    if logtime; H(i).YData(H(i).YData==0)=1; end
    if ~ismember(H(i).Color,[0 0 0],  'rows')
        H(i).Color = squeeze(TCLR(1,round(H(i).XData(1)),:))';
        H(i).LineWidth=1.5;
    end
end

% Taxa colours
axes(ax1);
imagesc(TCLR)
set(gca,'XTick',[],'YTick',[]);
ylabel('Taxon','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',16)

% RGB colours
axes(ax2);
imagesc(X_RGB)
set(gca,'XTick',[],'YTick',[]);
ylabel('RGB','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',16)

% Lon
axes(ax6);
imagesc(X_lon)
set(gca,'XTick',[],'YTick',[]);
ylabel('Longitude','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',16)
set(gca,'Colormap',hsv(numel(X_lon)));

% Latitude
axes(ax5);
imagesc(X_lat)
set(gca,'XTick',[],'YTick',[]);
ylabel('Latitude','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',16)
set(gca,'Colormap',turbo(numel(X_lat)));

% Trait 1
axes(ax4);
imagesc(X_t1)
set(gca,'XTick',[],'YTick',[]);
ylabel(Tname1,'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',16)
if strmatch(Tname1,'Trophic')
    set(gca,'Colormap',greenmag(numel(X_t1)));
elseif strmatch(Tname1,'Thermal')
    set(gca,'Colormap',redblue(numel(X_t1)));
end
caxis([min(T1) max(T1)])

% Trait 2
axes(ax3);
imagesc(X_t2)
set(gca,'XTick',[],'YTick',[]);
ylabel(Tname2,'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',16)
set(gca,'Colormap',jet(numel(X_t2)));
caxis(log10([min(T2) max(T2)]))

%%
Dsquare=squareform(distances./2);

f4=figure(4);
f4.Position=[154 401 800 800];
set(f4,'color','w');   
clf
ax0 = axes('Position',[0.15 0.2 0.70 0.70]);

imagesc(ax0,Dsquare(outperm,outperm))
set(ax0,'XTick',[],'YTick',[]);

ax1 = axes('Position',[0.150 0.170 0.700 0.025]); % Lat
ax2 = axes('Position',[0.150 0.145 0.700 0.025]); % Lon
ax3 = axes('Position',[0.150 0.905 0.700 0.025]); % TCLR
ax4 = axes('Position',[0.150 0.930 0.700 0.025]); % RGB
ax5 = axes('Position',[0.120 0.200 0.025 0.700]); % T1
ax6 = axes('Position',[0.855 0.200 0.025 0.700]); % T2

imagesc(ax1,X_lat)
set(ax1,'XTick',[],'YTick',[]);
ylabel(ax1,'Latitude','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',16)
ax1.Colormap=turbo(numel(parameters.ocn_pars.lon));

imagesc(ax2,X_lon)
set(ax2,'XTick',[],'YTick',[]);
ylabel(ax2,'Longitude','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',16)
ax2.Colormap=hsv(numel(parameters.ocn_pars.lat));

imagesc(ax3,TCLR)
set(ax3,'XTick',[],'YTick',[]);
ylabel(ax3,'Taxa','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',16)

imagesc(ax4,X_RGB)
set(ax4,'XTick',[],'YTick',[]);
ylabel(ax4,'RGB','Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',16)

imagesc(ax5,X_t1')
set(ax5,'XTick',[],'YTick',[]);
ylabel(ax5,Tname1,'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',16)
if strmatch(Tname1,'Trophic')
    ax5.Colormap=greenmag(numel(parameters.eco_pars.T1));
elseif strmatch(Tname1,'Thermal')
    ax5.Colormap=redblue(numel(parameters.eco_pars.T1));
end
caxis([min(T1) max(T1)])

imagesc(ax6,X_t2')
set(ax6,'XTick',[],'YTick',[]);
ylabel(ax6,Tname2,'Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','FontSize',16)
ax6.YAxisLocation='Right';
ax6.Colormap=turbo(numel(parameters.eco_pars.T1));

caxis(log10([min(T2) max(T2)]))

hc=colorbar(ax0,'SouthOutside','Position',[0.150 0.100 0.700 0.025]);
hc.Label.String='Divergence (years)';
hc.FontSize=16;


%%


for i=1:nrgb
    taxa_grid=zeros(size(biomass));
    taxa_grid(isextant(leafOrder))=squeeze(TCLR(1,:,i));
    taxon_RGB(:,:,:,:,i)=taxa_grid;
end


% rearrange dimensions - order: lat, T2, lon, T1
taxa_map=permute(taxon_RGB,[2 4 1 3 5]);
% combine size+lon (rows) and troph+lat (columns)
taxa_map=reshape(taxa_map,nlat*nT2,nlon*nT1,nrgb);

% rearrange dimensions - order: lat, T2, lon, T1
biomass_map=permute(biomass,[2 4 1 3]);
% combine size+lon (rows) and troph+lat (columns)
biomass_map=reshape(biomass_map,nlat*nT2,nlon*nT1);

biomass_map(biomass_map<threshold)=threshold;
land_mask = (log10(biomass_map)-log10(threshold))./(max(log10(biomass_map(:)))-log10(threshold));
% land_mask(isnan(land_mask))=1;

f5 = figure(5);                                                             % change active figure
f5 = clf(f5);                                                               % clear figure
% f5.Position=fig1pos;                                                      % Position figure
set(f5,'color','w');                                                        % set figure background to white




% Plot individual populations as T1xT2 matrix
ax1a = axes('Parent',f5);                                                   % setup axes
imagesc(taxa_map,'AlphaData',land_mask);                                % plot biomass

set(gca,'color',[0 0 0]);                                              % set axis background colour to grey for land
axis xy                                                                     % flip vertical axis
ax1a.Position=[1/(nT1+1)+0.05 ...
    1/(nT2+1)+0.05 ...
    nT1/(nT1+1)-0.125 ...
     nT2/(nT2+1)-0.1];                                            % expand axes within window
ax1a.XTick = xticksmiddle;                                                  % set x ticks
ax1a.YTick = yticksmiddle;                                                  % set y ticks
ax1a.TickLength=[0 0];                                                      % make ticks invisble
ax1a.XTickLabel=num2str('');                                                % no x tick labels
ax1a.YTickLabel=num2str('');                                                % no y tick labels
hold on                                                                     % hold figure for plotting lines
plot([xticks(1) xticks(end)],[yticks;yticks],'color',lncolr_phy,'LineW',1); % plot horizontal grid lines
plot([xticks;xticks],[yticks(1) yticks(end)],'color',lncolr_phy,'LineW',1); % plot vertical grid lines
% caxis([log10(threshold) log10(ceil(nanmax(PHY_a(:))))])                     % set colour scale

title('Taxa','FontSize',16)                                % add title
hold off

%%

%%
figure(1)
clf
set(gcf,'Color','w')
set(gca,'Position',[0 0 1 1])

[sx,sy,sz]=sphere;
scl=0.01;

axis equal
axis(repmat([-scl 1+scl],[1 3]))
axis off
camlight
set(gca, 'Projection','perspective')
hold on

for i=1:size(X_RGB,2)
    x=sx.*scl+X_RGB(1,i,1);
    y=sy.*scl+X_RGB(1,i,2);
    z=sz.*scl+X_RGB(1,i,3);

    sh=surf(x,y,z);
    sh.FaceColor=squeeze(TCLR(1,i,:))';
    sh.EdgeColor='none';

    if rem(i,10)==0
        drawnow
    end
end


%%


taxon_ID=zeros(size(biomass));
taxon_ID(isextant(leafOrder)) = taxa;

% % rearrange dimensions - order: lat, T2, lon, T1
% taxon_ID_map=permute(taxon_ID,[2 4 1 3]);
% % combine size+lon (rows) and troph+lat (columns)
% taxon_ID_map=reshape(taxon_ID_map,nlat*nT2,nlon*nT1);

for i=1:nlat
    % TAXONOMIC DIVERSITY
    % get all taxon IDs at latitude i (lon, T1 and T2);
    tmp = squeeze(taxon_ID(:,i,:,:));
    % PHENOTYPIC DIVERSITY
    % count number of non-zero phenotypes across longitudes (dimension 1 of tmp)
    pheno_lat(i)=sum(any(tmp,1),'all');
    % T1 DIVERSITY
    % count number of non-zero phenotypes across longitudes and T2 (dimensions 1 and 3 of tmp)
    T1_lat(i)=sum(any(squeeze(any(tmp,1)),2));
    % T2 DIVERSITY
    % count number of non-zero phenotypes across longitudes and T2 (dimensions 1 and 2 of tmp)
    T2_lat(i)=sum(any(squeeze(any(tmp,1)),1));

    % find non-zero and non-NaN taxa within latitude
    tmp2 = tmp(find(tmp~=0 & ~isnan(tmp)));
    % taxonomic diversity is the number of unique taxonomic IDs
    taxon_lat(i)=numel(unique(tmp2));
    


    for j=1:nlon
         % get all taxon IDs at latitude i and longitude j (T1 and T2);
        local_taxa_grid = squeeze(taxon_ID(i,j,:,:));
        % find non-zero and non-NaN taxa within grid cell
        extant_local_taxa = local_taxa_grid(find(local_taxa_grid~=0 & ~isnan(local_taxa_grid)));

        % phenotypic diversity is the number of occupied phenotypes
        pheno_div(i,j)=numel(extant_local_taxa);
        % taxonomic diversity is the number of unique taxonomic IDs
        taxon_div(i,j)=numel(unique(extant_local_taxa));

        % thermal optima diversity
        T2_div(i,j)=sum(any(local_taxa_grid,1));
        % size class diversity
        T1_div(i,j)=sum(any(local_taxa_grid,2));
    end
end

maxdiv1=max([pheno_lat';taxon_lat';pheno_div(:);taxon_div(:)]);

maxdiv2=max([T1_div(:)' T2_div(:)' taxon_div(:)' pheno_div(:)']);
%% Rank Abundance

figure(99)
clf
for i=1:nlat
    for j=1:nlon
         % get all biomasses at latitude i and longitude j (T1 and T2);
        local_biomass_grid = squeeze(biomass(i,j,:,:));
        Vol = unique(parameters.eco_pars.V)';

        local_abundance_grid = local_biomass_grid./Vol;

         % get all taxon IDs at latitude i and longitude j (T1 and T2);
        local_taxa_grid = squeeze(taxon_ID(i,j,:,:));


        unq_taxa = unique(local_taxa_grid);
        unq_taxa = unq_taxa(unq_taxa>0)';

        if numel(unq_taxa)>0
            clear local_RankAbund
            for k=1:numel(unq_taxa)
                local_RankAbund(k) = sum(local_abundance_grid(local_taxa_grid==unq_taxa(k)));

            end
            local_RankAbund = sort(local_RankAbund,'descend');
            RankAbund(i,j,1:k) = local_RankAbund;

            loglog(local_RankAbund,'Color',[0 0 0 0.1])
            hold on
        end
    end
end


%%
figure(2)
clf
set(gcf,'Color','w',...
        'Position',[27 534 1420 803])

nclrs=10*ceil(maxdiv2/10);
nyaxs=10*ceil(maxdiv1/10);
lnclr=[0 0.6 0.6];

coloraxis=[0 nclrs];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Taxonomic Diversity
ax1=axes;
axesm ('MapProjection','putnins5',...
       'Frame', 'off',...
       'Grid', 'off',...
       'MapLonLimit',parameters.ocn_pars.lon_edges([1 end]));
pcolorm(parameters.ocn_pars.lat_edges,parameters.ocn_pars.lon_edges,taxon_div')
ch=colorbar('Location','westoutside');
set(ax1,'position',[0.035 0.55 0.4 0.39])
axis off
title('Taxonomic Diversity')
caxis(coloraxis)
colormap(ax1,turbo(nclrs))
ch.Position(1:2)=ch.Position(1:2)+[0.025 -0.2];

ax2=axes;
stairs(parameters.ocn_pars.lat_edges,taxon_lat([1:end end]),'Color',lnclr,'LineW',1.5)
xlim(parameters.ocn_pars.lat([1 end]))
set(ax2,'position',[0.425    0.5925    0.075    0.306],...
    'linewidth',1.5,...
    'XTick',-90:30:90,...
    'YTick',0:20:120,...
    'YAxisLocation','Right',...
    'XDir','Reverse',...
    'Color','none')
box on
xlim(parameters.ocn_pars.lat_edges([1 end]))
ylim([0 nyaxs])
ylabel('n Taxa')
camroll(-90)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total Phenotypic Diversity
ax3=axes;
axesm ('MapProjection','putnins5',...
       'Frame', 'off',...
       'Grid', 'off',...
       'MapLonLimit',parameters.ocn_pars.lon_edges([1 end]));
pcolorm(parameters.ocn_pars.lat_edges,parameters.ocn_pars.lon_edges,pheno_div')
set(ax3,'position',[0.5 0.55 0.4 0.39])
axis off
title('Phenotypic Diversity')
caxis(coloraxis)
colormap(ax3,turbo(nclrs))
% ch=colorbar('Location','westoutside');
% ch.Position(1)=ch.Position(1)+0.05;

ax4=axes;
stairs(parameters.ocn_pars.lat_edges,pheno_lat([1:end end]),'Color',lnclr,'LineW',1.5)
set(ax4,'position',[0.895    0.5925    0.075    0.306],...
    'linewidth',1.5,...
    'XTick',-90:30:90,...
    'YTick',0:20:120,...
    'YAxisLocation','Right',...
    'XDir','Reverse',...
    'Color','none')
xlim(parameters.ocn_pars.lat_edges([1 end]))
ylim([0 nyaxs])
ylabel('n Phenotypes')
camroll(-90)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T1 Diversity
ax5=axes;
axesm ('MapProjection','putnins5',...
       'Frame', 'off',...
       'Grid', 'off',...
       'MapLonLimit',parameters.ocn_pars.lon_edges([1 end]));
pcolorm(parameters.ocn_pars.lat_edges,parameters.ocn_pars.lon_edges,T1_div')
set(ax5,'position',[0.035 0.15 0.4 0.39])
axis off
title([Tname1 ' Diversity'])
caxis(coloraxis)
colormap(ax5,turbo(parameters.eco_pars.nT1))
% ch=colorbar('Location','westoutside');
% ch.Position(1)=ch.Position(1)+0.05;

ax6=axes;
stairs(parameters.ocn_pars.lat_edges,T1_lat([1:end end]),'Color',lnclr,'LineW',1.5)
set(ax6,'position',[0.425    0.205    0.075    0.28],...
    'linewidth',1.5,...
    'XTick',-90:30:90,...
    'YAxisLocation','Right',...
    'XDir','Reverse',...
    'Color','none')
xlim(parameters.ocn_pars.lat_edges([1 end]))
box on
ylim([0 nyaxs])
ylabel('n Phenotypes')
camroll(-90)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T2 Diversity
ax7=axes;
axesm ('MapProjection','putnins5',...
       'Frame', 'off',...
       'Grid', 'off',...
       'MapLonLimit',parameters.ocn_pars.lon_edges([1 end]));
pcolorm(parameters.ocn_pars.lat_edges,parameters.ocn_pars.lon_edges,T2_div')
set(ax7,'position',[0.5 0.15 0.4 0.39])
axis off
title([Tname2 ' Diversity'])
caxis(coloraxis)
colormap(ax7,turbo(parameters.eco_pars.nT2))
% ch=colorbar('Location','westoutside');
% ch.Position(1)=ch.Position(1)+0.05;

ax8=axes;
stairs(parameters.ocn_pars.lat_edges,T2_lat([1:end end]),'Color',lnclr,'LineW',1.5)
set(ax8,'position',[0.895    0.205    0.075    0.28],...
    'linewidth',1.5,...
    'XTick',-90:30:90,...
    'YAxisLocation','Right',...
    'XDir','Reverse',...
    'Color','none')
xlim(parameters.ocn_pars.lat_edges([1 end]))
box on
ylim([0 nyaxs])
ylabel('n Phenotypes')
camroll(-90)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(findall(gcf,'-property','FontSize'),'FontSize',14)


