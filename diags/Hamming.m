function Hamming(output_dir,threshold,nsp,nclust)


if nclust>nsp & nsp~=0
    error('number of clusters must be <= number of terminal leaf nodes')
end
disp(['Processing output from ' output_dir])

% lget matObj for matfile data
matObj = matfile([output_dir '/OUTPUT.mat'])
eco_pars=matObj.eco_pars;
gen_pars=matObj.gen_pars;

% load saved data
load('/Users/baw103/GitHub/omg/TM_data/rw_lo_annual/A.mat')
load([output_dir '/matrix_vars.mat'])
vi=v_index.i(Ib);
vj=v_index.j(Ib);
vk=v_index.k(Ib);

gen_fcns=general_functions;

% load('/Users/baw103/GitHub/omg/TM_data/rw_lo_annual/OCEAN.mat')

%% Define grid coordinates
% get grid edge coordinates
Xe=OCEAN.lon_edges([unique(vi);max(vi)+1]);
Ye=OCEAN.lat_edges([unique(vj);max(vj)+1]);
% mesh grid of edges
[xe ye]=meshgrid(Xe,Ye);

% get grid midpoint coordinates
Xm=OCEAN.lon(unique(vi));
Ym=OCEAN.lat(unique(vj));
% mesh grid of midpoints
[xm ym]=meshgrid(Xm,Ym);


% calculate water flow vectors
B=A(Ib,Ib);
[xlon ylat]=meshgrid(OCEAN.lon(3:end-2),OCEAN.lat);
xn=B'*xlon(:);
yn=B'*ylat(:);
u=xn-xlon(:);
v=yn-ylat(:);


%% Prepare figure axes
% figure(1)
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
% drawnow
% gap = [0.01 0.005];
% marg_h = 0.030;
% marg_w = 0.02;

% color axis limits
maxPO4=0.1;
minPHY=0.0001;
maxPHY=1;


nphy=eco_pars.nsize*eco_pars.ntroph;


% make colormap for clusters (latin hypercube)
cmap=lhsdesign(nclust,3);
%%

PO4=matObj.PO4;
i_timeslice=find(~cellfun(@isempty,PO4))'; % find years with output data

disp(['Processing output from ' num2str(numel(i_timeslice)) ' time-slices'])
if ~exist([output_dir '/Figures/'],'dir')
    wd=cd(output_dir);
    !mkdir Figures
    cd(wd)
end


%% Calculate global metacommunity distance matrix
% maybe need to find a way to do this locally/regionally as well
nyr=i_timeslice(end);


disp(['Year ' num2str(gen_pars.save_output(nyr)+0.5)])
% for nyr=fliplr(i_timeslice)%(end)

PHY=cell2mat(matObj.PHY(nyr,1));
usedata = find(PHY>threshold);    % Extract indices of non-negligible biomass 

disp(['Analysing ' num2str(numel(usedata)) ' extant local populations.'])

% define EiE matrix grid coordinates
ESD     = (6.*eco_pars.V./pi)'.^(1/3);
ESDgrid = repmat(ESD,[size(PHY,1) 1]);
TROgrid = 1-repmat(eco_pars.trophic',[size(PHY,1) 1]);
LATgrid=repmat(vj,[1 eco_pars.ntroph*eco_pars.nsize]);
LONgrid=repmat(vi,[1 eco_pars.ntroph*eco_pars.nsize]);

% threshold = sum(PHY,2)./100; % threshold defined relative to community biomass at each location
% threshold = 1e-3;


% get trait and geographic coordinates of useable data
H_sz = ESDgrid(usedata);
H_tr = TROgrid(usedata);
H_lt = LATgrid(usedata);
H_ln = LONgrid(usedata);

N_sz=(log10(H_sz)-min(log10(H_sz)))./range(log10(H_sz));
N_tr=H_tr;
N_lt=(H_lt-min(H_lt))./range(H_lt);
N_ln=(H_ln-min(H_ln))./range(H_ln);

[~,~,xind]=unique(H_ln);
[~,~,yind]=unique(H_lt);
[~,~,sind]=unique(H_sz);
[~,~,tind]=unique(H_tr);


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HACKHACKHACKHACKHACKHACKHACKHACKHACKHACKHACKHACKHACKHACKHACKHACKHACKHACK!
% hackish compensation for slowdown of binary clock
% this is probably not publishable
alpha =      0.01804;
beta  =       11.44;
t     = beta.*D ./ (alpha .* (beta - D));
% HACKHACKHACKHACKHACKHACKHACKHACKHACKHACKHACKHACKHACKHACKHACKHACKHACKHACK!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create tree
% tree = linkage(D);
tree = linkage(t,'average');
%%
% Probably too many large eigenvalues to use Multidimensional Scaling
% [Y e] = cmdscale(D); % use multi-dimensional scaling to distribute nodes in 3D (colour) space according to genetic distance
% Y = (Y-min(Y))./range(Y); % normalise
% cmap=Y(:,1:3);
% xyz =Y(:,4:6);

% number of clusters to differentiate

% define cutoff as value between nclust nodes 
cutoff = median([tree(end-(nclust-1),3) tree(end-(nclust-2),3)]);

% generate tree
figure(1)
clf
subplot(1,3,1)
[~,rankorder] = sortrows([H_tr H_sz],'descend'); % rank nodes by trophic strategy and size
[leafOrder]   = sort_branches(tree,rankorder);
[H, nodes, outperm] = dendrogram(tree,nsp,'ColorThreshold',cutoff,'Orientation','left','ReOrder',leafOrder);
% [H, nodes, outperm] = dendrogram(tree,nsp,'Orientation','left'); % 
set(gcf,'Position',[53         523        1338         822]);
ax0=gca;
set(ax0,'YTick',[],...
        'YDir','reverse');

np=max(outperm);

% get nclust clusters
T = cluster(tree,'maxclust',nclust);
% get unique xcoord-node-cluster triplets
ync=unique([nodes T],'rows');
ync=[(1:nsp)' ync(outperm,:)];

hold on
xl=xlim;

H_orig=H; % save original dendrogram image data (H=H_orig;)
Dcrop =0;

% change colours to match maps
for i=1:numel(H)
    if ~ismember(H(i).Color,[0 0 0],'rows') % if node is not black
        % find an integer node within range of branch i
        ii=round(median(H(i).YData));
        % find cluster/s corresponding to branch
        ci=T(outperm(ii));
        H(i).Color=cmap(unique(ci),:); % change colour to row of cmap
    end
    H(i).LineWidth=0.75;
    %     squash branches after x=Dcrop onto x=Dcrop
    %     H(i).XData = max(H(i).XData,[Dcrop Dcrop Dcrop Dcrop]);
end
xlim([Dcrop xl(2)])

spstr=cellstr(strcat(num2str(sind,'%02i'),num2str(tind,'%02i')));
[H,spind] = prune_tree(tree,spstr,H);

set(gca,'FontSize',14)
xlabel('Estimated divergence (years)','FontSize',14)
title({'(a) Phylogenetic tree derived from neutral genome',' '})
%% Grow some fruit on the tree!
m=size(tree,1);
dots=[];
tic
for i=1:m % loop through branches as appearing in tree and H
    iparent=m+1+i; % find parent branch
    [ip,~,~]=find(tree(:,1:2)==iparent); % find parent branch
    % if branch is collapsed (i.e. [0.9 0.9 0.9]) and parent is not
    if ismember(H(i).Color,[0.9 0.9 0.9],'rows') && ~ismember(H(ip).Color,[0.9 0.9 0.9],'rows') % if collapsed and parent is not...
%         return
        ii=round(median(H(i).YData)); % find an integer node within range of child branch
        ci=T(outperm(ii)); % get color from lhc colormap
        cmp=cmap(ci,:);
        
        dots=[dots;[max(H(i).XData) median(H(i).YData) cmp]];
    % or if branch is not collapsed but one or both nodes are terminal (x==0)
    elseif ~ismember(H(i).Color,[0.9 0.9 0.9],'rows') && min(H(i).XData)==0

        ii=round(median(H(i).YData));
        ci=T(outperm(ii));    
        cmp=cmap(ci,:);
        
        tln = find(H(i).XData==0); % identify terminal leaf nodes
        
        cmp=repmat(cmp,numel(tln),1);
        dots=[dots;[H(i).XData(tln)' H(i).YData(tln)' cmp]];
    end
end

% scatter(dots(:,1),dots(:,2),7.5,dots(:,3:5),'filled');

[Bdt,Idt]=sort(dots(:,2));
toc

%%

Clstrs=zeros(size(PHY));
Clstrs(usedata)=T; % map clusters


clear nodemap phytmap 
jp=0;
for i=1:eco_pars.nsize
    for j=1:eco_pars.ntroph
        jp=jp+1;
        clr0 = gen_fcns.v2f(Clstrs(:,jp),vi,vj,vk);
        nodemap{i,j}(:,:) = squeeze(clr0(end,:,:));
        
        phy0 = gen_fcns.v2f(PHY(:,jp),vi,vj,vk);
        phytmap{i,j}(:,:) = squeeze(phy0(end,:,:));
    end
end
nodemap=cell2mat(nodemap);
phytmap=cell2mat(phytmap);

xtck=linspace(1,size(nodemap,2)+1,eco_pars.nsize+1);
ytck=linspace(1,size(nodemap,1)+1,eco_pars.ntroph+1);
othr=[zeros(1,eco_pars.nsize+1);ones(1,eco_pars.nsize+1)];

figure(2)
isz=find(nodemap<=0);
isn=find(isnan(nodemap));
nodemap(isz)=NaN;
nodemap(isn)=-inf;
pcolor(nodemap)
shading flat
axis xy
caxis([0 nclust])
colormap([0.75 0.75 0.75; cmap])
hold on
plot([xtck;xtck],othr.*size(nodemap,1),'k-')
plot(othr.*size(nodemap,2),[ytck;ytck],'k-')
set(gca,'XTick',[],'YTick',[])

fname=[output_dir '/Figures/OMG_Clusters_' num2str(nyr,'%04i') '.png'];
export_fig(fname,'-r300','-transparent')


figure(3)
isz=find(phytmap<=threshold);
isn=find(isnan(phytmap));
logphyt=log10(phytmap);
logphyt(isz)=NaN;
logphyt(isn)=-inf;
pcolor(logphyt)
hold on
plot([xtck;xtck],othr.*size(nodemap,1),'k-')
plot(othr.*size(nodemap,2),[ytck;ytck],'k-')
set(gca,'XTick',[],'YTick',[])
shading flat
axis xy
colormap([0.75 0.75 0.75; parula])
colorbar

fname=[output_dir '/Figures/OMG_Biomass_' num2str(nyr,'%04i') '.png'];
export_fig(fname,'-r300','-transparent')
%%

figure(1)

subplot(1,3,1)
% Make colormap
x=[0 0 0.5 1 1]'; % x cordinates of corners and centre
y=[0 1 0.5 0 1]'; % y cordinates of corners and centre
r=[0 1 1 0 1]';   % corresponding r values
g=[1 1 1 0 0]';   % corresponding g values
b=[0 0 1 1 1]';   % corresponding b values
xl=linspace(0,1,31); % unique values of trait x
yl=linspace(0,1,31); % unique values of trait y
[xx yy]=meshgrid(xl,yl); % generate trait grid
Fr=scatteredInterpolant(x,y,r); % r values across grid
Fg=scatteredInterpolant(x,y,g); % g values across grid
Fb=scatteredInterpolant(x,y,b); % b values across grid
rgb=[Fr(xx(:),yy(:)) Fg(xx(:),yy(:)) Fb(xx(:),yy(:))]; % rgb array

trind=N_tr.*30+1;
szind=N_sz.*30+1;
data=(trind-1).*31+szind;
ax1=axes;
imagesc(data(outperm))
ax1.Position=[0.3585 0.11 0.015 0.8150];
colormap(ax1,rgb)
box on
ax1.XTick=[];
ax1.YTick=[];

ax2=axes;
imagesc(Bdt)
colormap(ax2,dots(Idt,3:5))
ax2.Position=[0.3435 0.11 0.015 0.8150];
box on
ax2.XTick=[];
ax2.YTick=[];

%%
ax3=axes;
image(permute(reshape(rgb,[31,31,3]),[2 1 3]));
axis square

ax3.Position=[0.125 0.49    0.04    0.8150];
set(gca,'XTick',[1 16 31],'XTickLabel',{'A','M','H'})
set(gca,'YTick',[1 31],'YTickLabel',{'Small','Large'})
tp=[1 16 31];
set(ax3,'XAxisLocation','bottom','YAxisLocation','left','FontSize',12)
axis xy
xlabel({'Trophic','strategy'},'FontSize',12)
ylabel('Size (\mum)','FontSize',12)




%%

% Geographically-orientated food-web plot

% place trait coordinates within spatial grid
xfw=xind+N_tr.*0.9+0.05;
yfw=yind+N_sz.*0.9+0.05;

%
[xx yy]=meshgrid(Xe,Ye);
lon=interp1(1:numel(Xe),Xe,xfw);
lat=interp1(1:numel(Ye),Ye,yfw);

biomass=full(PHY(usedata));

data=matObj.PHY;
data=full(data{nyr});
data=sum(data,2);
data(data==0)=NaN;
data=gen_fcns.v2f(data,vi,vj,vk);
data=squeeze(data(end,:,:));
Cbio=data(:,3:end); % remove land cells


cls=flipud(ync(:,2));
figure(1)
ax0=subplot(10,6,[3 42])
ax0.Position=[0.3991-0.03    0.3615+0.075    0.4730    0.5];
ax=axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLonLim',minmax(Xe'));
pcolorm(Ye,Xe,log10(Cbio))
plotm(yy,xx,'-k')
plotm(yy',xx','-k')




hc=colorbar('SouthOutside');
set(hc,...
    'Ticks',log10([0.001 0.002 0.005 0.01 0.02 0.03 0.05 0.1 0.2 0.3]),...
    'TickLabels',{'0.001','0.002','0.005','0.01','0.02','0.03','0.05','0.1','0.2','0.3'},...
    'FontSize',12,...
    'Position',[0.3945 0.4106 0.4215 0.0195])
caxis(log10([0.001 0.3]))
cmp=flipud(gray(128));
colormap(ax,cmp(1:64,:))

scatterm(lat,lon,sqrt(biomass).*200,cmap(T,:),'filled')
% scatterm(lat,lon,sqrt(biomass)'.*300,'k','LineWidth',0.1)

% grid on
axis normal
axis off
title('(b) Genetic clustering within local food-webs','FontSize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ndmap =nodemap;
phymap=phytmap;
rg=1:18;
cg=3:16;

ax4=subplot(7,5,[28 33]);
ax4.Position=[0.4556-0.030-0.03    0.1100    0.1237    0.2086];
ax=axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLonLim',minmax(Xe'));
pcolorm(Ye,Xe,ndmap(7*18+rg,cg))
hold on
biom=phymap(7*18+rg,cg);
% contourm(Ym,Xm,log10(biom),-6:0,'w')
uB=u.*sqrt(biom(:));
vB=v.*sqrt(biom(:));
qh=quiverm(ylat(:)-vB,xlon(:)-uB,vB,uB,'w');
qh(1).LineWidth=1;qh(2).LineWidth=1;
axis xy
plotm(yy,xx,'-k')
plotm(yy',xx','-k')
caxis([0 nclust])
colormap(ax4,[0.2 0.2 0.2; cmap])
grid on
axis off
set(gca,'XTickLabel','','YTickLabel','',...
        'XTick',0.5:14.5,'YTick',0.5:18.5)
title({['(c) 0.88 ' char(956) 'm Phytoplankton'],' '},'FontSize',14)

ax5=subplot(7,5,[29 34]);
ax5.Position=[0.618-0.045-0.03    0.1100    0.1237    0.2086];
ax=axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLonLim',minmax(Xe'));
pcolorm(Ye,Xe,ndmap(19*18+rg,21*16+cg))
hold on
biom=phymap(19*18+rg,21*16+cg);
% contourm(Ym,Xm,log10(biom),-6:0,'w')
uB=u.*sqrt(biom(:));
vB=v.*sqrt(biom(:));
qh=quiverm(ylat(:)-vB,xlon(:)-uB,vB,uB,'w');
qh(1).LineWidth=1;qh(2).LineWidth=1;
axis xy
plotm(yy,xx,'-k')
plotm(yy',xx','-k')
caxis([0 nclust])
colormap(ax5,[0.2 0.2 0.2; cmap])
grid on
axis off
set(gca,'XTickLabel','','YTickLabel','',...
        'XTick',0.5:14.5,'YTick',0.5:18.5)
title({['(d) 88 ' char(956) 'm Mixotroph'],' '},'FontSize',14)

ax6=subplot(7,5,[30 35]);
ax6.Position=[0.7813-0.06-0.03    0.1100    0.1237    0.2086];
ax=axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLonLim',minmax(Xe'));
pcolorm(Ye,Xe,ndmap(28*18+rg,30*16+cg))
hold on
biom=phymap(28*18+rg,30*16+cg);
% contourm(Ym,Xm,log10(biom),-6:0,'w')
uB=u.*sqrt(biom(:));
vB=v.*sqrt(biom(:));
qh=quiverm(ylat(:)-vB,xlon(:)-uB,vB,uB,'w');
qh(1).LineWidth=1;qh(2).LineWidth=1;
axis xy
plotm(yy,xx,'-k')
plotm(yy',xx','-k')
caxis([0 nclust])
colormap(ax6,[0.2 0.2 0.2; cmap])
grid on
axis off
set(gca,'XTickLabel','','YTickLabel','',...
        'XTick',0.5:14.5,'YTick',0.5:18.5)
title({['(e) 2,700 ' char(956) 'm Zooplankton'],' '},'FontSize',14)


% fname=[output_dir '/Figures/OMG_Phylogeny_' num2str(nyr,'%04i') '.png'];
set(gcf,'color','w');
fname=['~/GoogleDrive/MatrixMetacommunityModel/Figures/Genetic_distance_' num2str(nyr,'%04i') '.png'];
export_fig(fname,'-r300')

%% 

H_lt = LATgrid(usedata);
H_ln = LONgrid(usedata);

BIO=PHY(usedata);

% table of location indices for each extant population
extloc=[H_ln,H_lt];
[~,IA,siteind] = unique(extloc,'rows'); % index of unique locations

[~,latsrt]=sortrows([vj vi]);

% BRAY-CURTIS
BC=zeros(size(B));
for ii=1:size(B,1)
    for jj=1:size(B,1)
        ii_ind=find(siteind==ii); % index for populations at site 1
        jj_ind=find(siteind==jj); % index for populations at site 2
        
        BIO1=full(BIO(ii_ind));   % biomass of populations at site 1
        BIO2=full(BIO(jj_ind));   % biomass of populations at site 2
        Btot=sum([BIO1;BIO2]);    % total biomass at two sites
        
        T1=T(ii_ind); % species (clusters) at site 1
        T2=T(jj_ind); % species (clusters) at site 2
        
        [T1,~,ic1]=unique(T1); % get index of unique species clusters at site 1
        [T2,~,ic2]=unique(T2); % get index of unique species clusters at site 2
        
        BIO1 = accumarray(ic1,BIO1); % accumulate biomass of non-unique species
        BIO2 = accumarray(ic2,BIO2); % accumulate biomass of non-unique species
        
        [~,I1,I2] = intersect(T1,T2); % site specific indices of intersection
        
        % Bray curtis dissimilarity index
        BC(ii,jj)=1-sum(2.*min([BIO1(I1) BIO2(I2)],[],2))./Btot;
    end
end
BC(BC<0)=0;


% BRAY-CURTIS (transport and functional group)
[~,pp]=meshgrid(1:numel(Ib),1:961);
pp=pp';
functyp=pp(usedata);
BCf=zeros(size(B));
for ii=1:size(B,1)
    for jj=1:size(B,1)
        ii_ind=find(siteind==ii); % index for populations at site 1
        jj_ind=find(siteind==jj); % index for populations at site 2
        
        BIO1=full(BIO(ii_ind));   % biomass of populations at site 1
        BIO2=full(BIO(jj_ind));   % biomass of populations at site 2
        Btot=sum([BIO1;BIO2]);    % total biomass at two sites
        
        P1=functyp(ii_ind); % species (phenotypes) at site 1
        P2=functyp(jj_ind); % species (phenotypes) at site 2
        
        [P1,~,ic1]=unique(P1); % get index of unique species clusters at site 1
        [P2,~,ic2]=unique(P2); % get index of unique species clusters at site 2
        
        BIO1 = accumarray(ic1,BIO1); % accumulate biomass of non-unique species
        BIO2 = accumarray(ic2,BIO2); % accumulate biomass of non-unique species
        
        [~,I1,I2] = intersect(P1,P2); % site specific indices of intersection 
        
        % Bray curtis dissimilarity index
        BCf(ii,jj)=1-sum(2.*min([BIO1(I1) BIO2(I2)],[],2))./Btot;
    end
end
BCf(BCf<0)=0;

%%

G=digraph(1./B); % shortestpath transport time
t_dist = distances(G); % transport time matrix

G=digraph(1./(B*(eye(numel(Ib)).*sum(PHY,2)))); % biomass weighted shortest path time
t_distB = distances(G); % transport time matrix


% eucd=sqrt((vi-vi').^2+(vj-vj').^2); % euclidian distance based on rectangular grid (i.e. wrong)
% latd=abs(ym(:)-ym(:)'); % absolute latitudinal distance

% Calculate shortest connected path between gridpoints
latrad=deg2rad(ym(:)); % latitude in radians
lonrad=deg2rad(xm(:)); % longitude in radians
% calculate central angle between all gridpoints
deltasig = acos(min(sin(latrad).*sin(latrad') + cos(latrad).*cos(latrad').*cos(lonrad-lonrad'),1));
r=6.371e3; % radius Earth (km)
grtcrc = r.*deltasig; % great-circle distance between all gridpoints
% index of all directly connected points
adj=abs(vi-vi')==1 & abs(vj-vj')==0 ... % 1 different in vi
  | abs(vi-vi')==0 & abs(vj-vj')==1 ... % or in vj
  | abs(vi-vi')==1 & abs(vj-vj')==1;    % or in both
grtcrc = grtcrc .* adj; % restrict to distances between contiguous points
G=digraph(grtcrc); % graph of connected points
shrtst = distances(G); % find shortest path through connected points

figure(33)
clf

subplot(221)
axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLatLim',minmax(Ye'),'MapLonLim',minmax(Xe'));
nsrc=1;
pcolorm(Ye,Xe,reshape(median(log10(t_dist),1),18,14)),colorbar
colormap(jet)
axis off
plotm(yy,xx,'-k')
plotm(yy',xx','-k')
quiverm(ylat(:)-v,xlon(:)-u,v,u,'w');
title('Median shortest transport path dispersal time (log_{10}[n timestep])')

subplot(222)
axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLatLim',minmax(Ye'),'MapLonLim',minmax(Xe'));
nsrc=1;
pcolorm(Ye,Xe,reshape(median(log10(t_dist),2),18,14)),colorbar
colormap(jet)
axis off
plotm(yy,xx,'-k')
plotm(yy',xx','-k')
quiverm(ylat(:)-v,xlon(:)-u,v,u,'w');
title('Median shortest transport path arrival time (log_{10}[n timestep])')

subplot(223)
axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLatLim',minmax(Ye'),'MapLonLim',minmax(Xe'));
nsrc=1;
tmp=reshape((median(t_dist,1)+median(t_dist',1))./2,18,14);
pcolorm(Ye,Xe,log10(tmp)),colorbar
colormap(jet)
axis off
plotm(yy,xx,'-k')
plotm(yy',xx','-k')
quiverm(ylat(:)-v,xlon(:)-u,v,u,'w');
title('Median shortest transport path arrival/dispersal time (log_{10}[n timestep])')

subplot(224)
axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLatLim',minmax(Ye'),'MapLonLim',minmax(Xe'));
nsrc=1;
pcolorm(Ye,Xe,reshape(median(shrtst,1),18,14)),colorbar
colormap(jet)
% caxis([0.2 0.35])
axis off
plotm(yy,xx,'-k')
plotm(yy',xx','-k')
quiverm(ylat(:)-v,xlon(:)-u,v,u,'w');
title('Median shortest grid-point connectivity (km)')

%%
cf=figure(44)
clf
set(cf,'Position',[58 501 1420 844]);


BCv=BC(:);
BCv(BCv==0)=NaN;

[Nbin,edges,ibin] = histcounts(shrtst(:)./1000,0:25);
for i=unique(ibin)'
    cellBC{i}=BCv(ibin==i);
end
subplot(122)
boxplot(BCv,ibin)
xlabel('Distance (10^3 km)')
ylabel('Bray-Curtis dissimilarity')

% subplot(222)
% [Nbin,edges,ibin] = histcounts(log10(t_dist(:)),25);
% boxplot(BCv,ibin)
% set(gca,'XTickLabel',num2str(round(10.^edges'),'%d'))
% xlabel('Shortest path time (log10(days))')
% ylabel('Bray-Curtis dissimilarity')

% Calculate diversity for each box
% species richness
Tgrid=zeros(numel(Ib),961);
Tgrid(usedata)=T;
for i=1:numel(Ib) % find number of unique clusters in each location
    ncloc(i)=numel(find(unique(Tgrid(i,:))));
    nPloc(i)=numel(find(Tgrid(i,:)));
end


% Calculate diversity for each latitude band
Pgrid=repmat(1:961,numel(Ib),1).*(Tgrid~=0);
zonal=[];
zc=[];
BCinlat=[];
BCfinlat=[];
for j=unique(vj)' % FOR EACH LATITUDINAL BAND
    jj=find(vj==j);
    nclat(j)=numel(find(unique(Tgrid(jj,:))));
    nPlat(j)=numel(find(unique(Pgrid(j,:))));
    
    BClat(j,:)=mean(BC(jj,jj),1);
    
    BCinlat=[BCinlat;mean(BC(vj==j,vj==j),1)];
    BCfinlat=[BCfinlat;mean(BCf(vj==j,vj==j),1)];
    
    zonal=[zonal;Tgrid(jj,any(Tgrid(jj,:)))';NaN.*ones(10,14)];
    zc=[zc size(zonal,1)];
end

% clf
% zonal(zonal==0)=NaN;
% pcolor(1:14,1:size(zonal,1),zonal)
% shading flat
% hold on
% plot([ones(size(zc));ones(size(zc)).*14],[zc;zc],...
%       'k','LineWidth',1)
% set(gca,'YTick',[1 zc],'YTickLabel',Ye)
% colormap(lhsmap)
% imagesc(BClat)







abslat=abs(ym(:));

hs=subplot(221);
hspos=hs.Position;
% [Nbin,edges,ibin] = histcounts(ncloc,0:25);
% boxplot(abs(ym(:)),ncloc)
% hold on

shvr1=randn(size(abslat))./3;
shvr2=randn(size(abslat))./3;
sz=sqrt(biom(:)).*500;

scatter(abslat+shvr1,ncloc'+shvr2,sz,BCinlat(:),'filled',...
        'MarkerEdgeColor','k','LineWidth',0.01)
hold on
plot(abs(Ym(unique(vj)')),nclat./3,'k','LineWidth',1)
colormap(parula),box on
cb=colorbar('SouthOutside');
cb.Label.String='Dissimilarity with same latitide';
cb.Label.FontSize=11;
set(gca,'Position',hspos)
cb.Position(2) = cb.Position(2) -0.03;
caxis([0 1])
colormap(jet)
ylabel('n genetic clusters')
xlabel('absolute latitude')
title('dot area \propto local biomass')

% for i=1:numel(Ib)
%     cl=unique(Tgrid(i,:));
%     cl=cl(find(cl));
%     x=abs(yn(i))+(1:numel(cl))./10;
%     y=repmat(ncloc(i),1,numel(cl));
%     scatter(x,y,10,lhsmap(cl,:),'filled')
% end


subplot(223)
scatter(abslat+shvr1,nPloc'+shvr2,sz,BCfinlat(:),'filled',...
        'MarkerEdgeColor','k','LineWidth',0.01)
hold on
plot(abs(Ym(unique(vj)')),nPlat,'k','LineWidth',1)
colormap(parula),box on
caxis([0 1])
colormap(jet)
ylabel('n phenotypes')
xlabel('absolute latitude')




%%
figure(22)
clf

% n genetic clusters
subplot(421)
ax=axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLonLim',minmax(Xe'));
pcolorm(Ye,Xe,reshape(ncloc,18,14))
contourm(Ym,Xm,reshape(full(sum(PHY,2)),18,14),5,'k')
axis off
caxis(minmax(ncloc))
quiverm(ylat(:)-v,xlon(:)-u,v,u,'w');
title('n genetic clusters')
colormap(jet)
colorbar

% n phenotypic groups
subplot(422)
ax=axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLonLim',minmax(Xe'));
pcolorm(Ye,Xe,reshape(nPloc,18,14))
contourm(Ym,Xm,reshape(full(sum(PHY,2)),18,14),5,'k')
axis off
caxis(minmax(nPloc))
quiverm(ylat(:)-v,xlon(:)-u,v,u,'w');
title('n phenotypic groups')
colorbar

% mean genetic Bray-Curtis
subplot(423)
ax=axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLonLim',minmax(Xe'));
pcolorm(Ye,Xe,reshape(mean(BC,1),18,14))
contourm(Ym,Xm,reshape(full(sum(PHY,2)),18,14),5,'k')
axis off
caxis([0 1])
quiverm(ylat(:)-v,xlon(:)-u,v,u,'w');
title('mean genetic Bray-Curtis dissimilarity')
colorbar

% mean phenotypic Bray-Curtis
subplot(424)
ax=axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLonLim',minmax(Xe'));
pcolorm(Ye,Xe,reshape(mean(BCf,1),18,14))
contourm(Ym,Xm,reshape(full(sum(PHY,2)),18,14),5,'k')
axis off
caxis([0 1])
quiverm(ylat(:)-v,xlon(:)-u,v,u,'w');
title('mean phenotypic Bray-Curtis dissimilarity')
colorbar

% mean shortest transport time
subplot(425)
axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLatLim',minmax(Ye'),'MapLonLim',minmax(Xe'));
nsrc=1;
pcolorm(Ye,Xe,reshape(log10(mean(t_dist,1)),18,14)),colorbar
axis off
quiverm(ylat(:)-v,xlon(:)-u,v,u,'w');
title({'Mean shortest transport path dispersal time ','(log_{10}[n timestep])'})
caxis([3.5 4])

subplot(426)
axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLatLim',minmax(Ye'),'MapLonLim',minmax(Xe'));
nsrc=1;
pcolorm(Ye,Xe,reshape(log10(mean(t_dist,2)),18,14)),colorbar
axis off
quiverm(ylat(:)-v,xlon(:)-u,v,u,'w');
title({'Mean shortest transport path arrival time','(log_{10}[n timestep])'})
caxis([3.5 4])

subplot(427)
axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLatLim',minmax(Ye'),'MapLonLim',minmax(Xe'));
nsrc=1;
tmp=reshape((mean(t_dist,1)+mean(t_dist',1))./2,18,14);
pcolorm(Ye,Xe,log10(tmp)),colorbar
axis off
quiverm(ylat(:)-v,xlon(:)-u,v,u,'w');
title({'Mean shortest transport path arrival/dispersal time','(log_{10}[n timestep])'})
caxis([3.5 4])

% mean shortest path through ocean grid (km)
subplot(428)
ax=axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLonLim',minmax(Xe'));
pcolorm(Ye,Xe,reshape(mean(shrtst,1),18,14))
% contourm(Ym,Xm,reshape(full(sum(PHY,2)),18,14),5,'k')
axis off
% rm(ylat(:)-v,xlon(:)-u,v,u,'w');
title('median shortest path through ocean grid (km)')
colorbar
caxis([8000 14000])



set(gcf,'Color','w')



%%
figure(99)

nclst=40;
lhsmap=lhsdesign(nclst,3);

% shortest grid-point connectivity
subplot(231)
tree=linkage(squareform(shrtst)); % convert distance matrix to match 'pdist' output
clst = cluster(tree,'MaxClust',nclst);
ax=axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLonLim',minmax(Xe'));
pcolorm(Ye,Xe,reshape(clst,18,14))
hold on
axis off
quiverm(ylat(:)-v,xlon(:)-u,v,u,'w');
title('shortest path through ocean grid')
colormap(lhsmap)

% transport shortest-path distance
subplot(232)
G=digraph(1./B); % graph inverse of transport
t_dist = distances(G); % transport shortest-path distance
D=mean([squareform(t_dist);squareform(t_dist')]);
tree=linkage(D);
clst = cluster(tree,'MaxClust',nclst);
[JJ,II]=sort(clst);
ax=axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLonLim',minmax(Xe'));
pcolorm(Ye,Xe,reshape(clst,18,14))
hold on
axis off
quiverm(ylat(:)-v,xlon(:)-u,v,u,'w');
title('transport shortest path')

% genetic Bray-Curtis community dissimilarity
subplot(233)
BC(find(eye(numel(Ib))))=0;
D=mean([squareform(BC);squareform(BC')]);
tree=linkage(D);
clst = cluster(tree,'MaxClust',nclst);
ax=axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLonLim',minmax(Xe'));
pcolorm(Ye,Xe,reshape(clst,18,14))
hold on
axis off
quiverm(ylat(:)-v,xlon(:)-u,v,u,'w');
title('genetic community dissimilarity')

% phenotypic Bray-Curtis community dissimilarity
subplot(234)
BCf(find(eye(numel(Ib))))=0;
D=mean([squareform(BCf);squareform(BCf')]);
tree=linkage(D);
clst = cluster(tree,'MaxClust',nclst);
ax=axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLonLim',minmax(Xe'));
pcolorm(Ye,Xe,reshape(clst,18,14))
hold on
axis off
quiverm(ylat(:)-v,xlon(:)-u,v,u,'w');
title('phenotypic community dissimilarity')

% distance between phenotypic pairs
subplot(235)
PHYext=PHY;
D = pdist(PHY,'euclidean'); % transport distance
tree=linkage(D);
clst = cluster(tree,'MaxClust',nclst);
ax=axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLonLim',minmax(Xe'));
pcolorm(Ye,Xe,reshape(clst,18,14))
hold on
axis off
quiverm(ylat(:)-v,xlon(:)-u,v,u,'w');
title('euclidean (biomass) distance between extant pairs')

    


%%
% function [H,spstr] = prune_tree(tree,spstr,H)
%     m = size(spstr,1);
%     for i=1:m-1
%         n1=tree(i,1);
%         n2=tree(i,2);
%         
%         if strcmp(spstr{n1},spstr{n2})
%             spstr{m+i}=spstr{n1};
%             H(i).Color=[1 1 1 0];
%         else
%             spstr{m+i}=0;
%         end
%     end
% end







