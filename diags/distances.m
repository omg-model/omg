function distances(output_dir)

nyr=1200;%i_timeslice(end);

disp(['Processing output from ' output_dir])

% lget matObj for matfile data
matObj = matfile([output_dir '/OUTPUT.mat'])
eco_pars=matObj.eco_pars;
gen_pars=matObj.gen_pars;

% load saved data
load('/Users/baw103/GitHub/omg/TM_data/rwlla_annual/surfaceTM.mat')
load([output_dir '/matrix_vars.mat'])
vi=v_index.i(Ib);
vj=v_index.j(Ib);
vk=v_index.k(Ib);

gen_fcns=general_functions;

load('/Users/baw103/GitHub/omg/TM_data/rwlla_annual/OCEAN.mat')

%% Define grid coordinates
% get grid edge coordinates
Xe=OCEAN.lon_edges(3:end-2);
Ye=OCEAN.lat_edges;
% mesh grid of edges
[xe ye]=meshgrid(Xe,Ye);

% get grid midpoint coordinates
Xm=OCEAN.lon(3:end-2);
Ym=OCEAN.lat;
% mesh grid of midpoints
[xm ym]=meshgrid(Xm,Ym);

%% Calculate and plot shortest path distances
% transport coordinate midpoints to determine effect of 1 transport step
x_1=B'*xm(:);
y_1=B'*ym(:);

% calculate transport vectors
u=x_1-xm(:);
v=y_1-ym(:);

% optinoal squishing
f_squish=1; % set to 1 for no effect, 0 gives direction only
u=sign(u).*nthroot(abs(u),f_squish);
v=sign(v).*nthroot(abs(v),f_squish);

% bi-directional graph of non-zero elements in transport matrix 
% Edges are effectively distance between nodes,
% so the 1./B is used, such that stronger fluxes have shorter distances (i.e. higher weights)
G=digraph(1./B,'OmitSelfLoops');
t_dist=distances(G); % matrix of all shortestpath distances

%% Plot transport fluxes
figure(1)
clf

% use a square number for nexamp
nexamp=16; % nuber of example cases
ip=round(linspace(1,252,nexamp));
for i=1:nexamp
    subplot(sqrt(nexamp),sqrt(nexamp),i)
    % get transport distances for point ippoint
	up_dist   = t_dist(ip(i),:)'; % vector of upstream distances
	dn_dist   = t_dist(:,ip(i)); % vector of arrival distance
	min_dist  =  min([up_dist dn_dist],[],2);
	mean_dist = mean([up_dist dn_dist],2);
    
    % reshape to map array
    map=reshape(min_dist,18,14);
    
    % plot
    ax=axesm ('flatplrs', 'Frame', 'on', 'Grid', 'off','MapLonLimit',[-140 140]);
    pcolorm(Ye,Xe,log10(map))
    axis tight
    axis off
    caxis([2 5])
    colorbar
    
    % plot flow vectors
    hold on
    quiverm(y_1,x_1,v,u,'k');
    
    % plotlocation of point ip
    plotm(y_1(ip(i)),x_1(ip(i)),'ro')
    drawnow
end
%% Plot median upstream and downstream times from and to all locations

figure(2)
clf

% plot median global upstream arrival time
subplot(211)
ax=axesm ('flatplrs', 'Frame', 'on', 'Grid', 'off','MapLonLimit',[-140 140]);
vect=median(t_dist,2);
map=reshape(vect,18,14);
pcolorm(Ye,Xe,log10(map))
hold on
% plot flow vectors
hq=quiverm(y_1,x_1,v,u,'k');
colorbar
caxis([3 5])
axis tight
axis off
title('Median global upstream arrival time')

% plot median global downstream arrival time
subplot(212)
ax=axesm ('flatplrs', 'Frame', 'on', 'Grid', 'off','MapLonLimit',[-140 140]);
vect=median(t_dist,1);
map=reshape(vect,18,14);
pcolorm(Ye,Xe,log10(map))
hold on
% plot flow vectors
hq=quiverm(y_1,x_1,v,u,'k');
colorbar
caxis([3 5])
axis tight
axis off
title('Median global downstream arrival time')



%% Plot global shortest-path distance matrix
figure(3)
clf

[~,Ilat]=sort(ym(:)); % sort spatial grid by latitude (default sort is by longitude)

imagesc(0:14:252,0:14:252,log10(t_dist(Ilat,Ilat)))
grid on
set(gca,'XTick',0:14:252,'XTickLabel',num2str(round(Ye)),...
        'YTick',0:14:252,'YTickLabel',num2str(round(Ye)));
axis xy
axis square
colorbar
title('global shortest-path distance matrix')
xlabel('Upstream latitude')
ylabel('Downstream latitude')

%% plot transport clustering, biomass and transport vectors
% load all plankton biomass for year = nyr
PHY0=cell2mat(matObj.PHY(nyr,1));

threshold = 1e-6; 
usedata = find(PHY0>threshold);    % Extract indices of non-negligible biomass 

% get trait info
ESD     = (6.*eco_pars.V./pi)'.^(1/3); % Diameter for all populations
ESDgrid = repmat(ESD,[size(PHY0,1) 1]); % Diameter grided by space
TROgrid = repmat(eco_pars.trophic',[size(PHY0,1) 1]); % Trophic class grided by space


% get trait and geographic coordinates of useable data
H_sz = ESDgrid(usedata);
H_tr = 1-TROgrid(usedata);
% lookup table of traits 
traits=[1-eco_pars.trophic ESD'];

[C,IA,IC] = unique([H_tr H_sz],'rows'); %

ncl=20;

figure(5)
clf

for iextant=1:numel(IA) % for all extant phenotypes
    % get phenotype index for extant population 'iextant'
    iphen=find(ismember(traits,C(iextant,:),'rows'));
    % extract biomass of phenotype iextant
    PHY=PHY0(:,iphen);
    % index biomass > threshold
    isgd=find(PHY>1e-6);
    
    % Calculate distances weighted by 1/(biomass transported)
    DX=1./(B.*PHY);    
    G=digraph(DX,'OmitSelfLoops');
    b_dist=distances(G);
    
    tree=linkage(b_dist(isgd,isgd));
    clst = cluster(tree,'MaxClust',ncl);
%     clst = cluster(tree,'Cutoff',1.15);ncl=numel(unique(clst))
    
    cmap=lhsdesign(ncl,3);
    subplot(5,6,iextant)
    ax=axesm ('gortho', 'Frame', 'on', 'Grid', 'off','MapLonLimit',[-140 140]);
    pcolorm(Ye,Xe,reshape(log10(PHY),18,14))
    caxis([log10(threshold) -1])
    hold on
    scatterm(ym(isgd),xm(isgd),10,cmap(clst,:),'filled')
    quiverm(y_1,x_1,v,u,'k');
    colormap(flipud(gray))
	axis off
    title({['#' num2str(iextant) ': '],[char(964) '_{auto} = ' num2str(C(iextant,1))],['ESD = ' num2str(C(iextant,2)) ' \mum']})
end

%% calculate genetic distances
GENOME=full(cell2mat(matObj.GENOME(nyr,1)));
nbit=53; % log2(FlIntMax) = 53
str='';
for i=1:eco_pars.ngenes
    bin_gene = dec2bin(GENOME(usedata,i),nbit); % convert double value to binary string
    str=[str bin_gene];
end
% convert to numeric values (1/0)
binstr=str-'0';
ww=nthroot(2,10).^(0:eco_pars.ngenes-1);
w=reshape(repmat(ww,nbit,1),1,[]); % weight genes to have different mutation rates
w=w./mean(w);
binstr=binstr./w; % multiply bits by inverse mutation weights (more likely worth less)

%% Plot genetic distance as a function of different metrics of transport distance
% e.g.
% raw shortest-path
% biomass weighted shortest path
% absolute latitudinal distance
% spatial euclidean distance

figure(6)
clf

for iextant=1:numel(IA) % for all extant phenotypes
    % get phenotype index for extant population 'iextant'
    iphen=find(ismember(traits,C(iextant,:),'rows'));
    % extract biomass of phenotype iextant
    PHY=PHY0(:,iphen);
    % index biomass > threshold
    isgd=find(PHY>1e-6);
    
    % Calculate transport distances
%     DX=1./B;        % unweighted
    DX=1./(B.*PHY); % weighted by 1/(biomass transported)
    G=digraph(DX,'OmitSelfLoops');
    t_dist=distances(G);
    t_dist=t_dist(isgd,isgd); % extract for non-zero locations

    % Calculate genetic distances
    spcstr=binstr(IC==iextant,:); % find genes of current phenotype
    g_dist = squareform(pdist(spcstr,'CityBlock')./size(spcstr,2));

%     % define full sized matrices filled with NaNs
%     t_full=zeros(size(B)).*NaN;
%     g_full=t_full;
%     % place t and g matrices in full ones
%     t_full(isgd,isgd)=t_dist;
%     g_full(isgd,isgd)=g_dist;
%     imagesc(Ilat,Ilat,log10([t_full(Ilat,Ilat) g_full(Ilat,Ilat)]))
%     grid on
%     set(gca,'XTick',0:14:252,'XTickLabel',num2str(round(Ye)),...
%         'YTick',0:14:252,'YTickLabel',num2str(round(Ye)));
%     axis xy
%     axis square
%     colorbar
%     title('global shortest-path distance matrix')
%     xlabel('Upstream latitude')
%     ylabel('Downstream latitude')

    
    subplot(5,6,iextant)
    biom=repmat(full(PHY(isgd)),numel(isgd),1);
    scatter(t_dist(:),g_dist(:),5,log10(biom(:)),'filled')
    set(gca,'XScale','Log','YScale','Log')
	box on
    ylabel('genetic')
    xlabel('shortest-path')
    axis square
%     axis([1e2 1e11 1e-4 1e2])
    caxis([log10(threshold) -1])
    colorbar
    
    logX{iextant} =log10(t_dist(:));
    logY{iextant} =log10(g_dist(:));
    
    isr=min([~isinf(logX{iextant}) ~isnan(logX{iextant}) ~isinf(logY{iextant}) ~isnan(logY{iextant})],[],2);
    
    % linear fir weighted by local population size
    [fitobject,gof] = fit(logX{iextant}(isr),logY{iextant}(isr),'poly1','Weights',biom(isr));
    hold on
    plot(10.^logX{iextant}(isr),10.^fitobject(logX{iextant}(isr)),'r','LineWidth',2)
    
    title(['r = '  num2str(gof.rsquare)])
    
end

%%
figure(7)
clf
figure(8)
clf
rgb=csquare(31);

nbin=10;
binedg=linspace(3,11,nbin+1);

cmax=max(cellfun('size',logX,1));
cmin=min(cellfun('size',logX,1));
clr=round(63.*(cellfun('size',logX,1)-cmin)./(cmax-cmin))+1;
cmap=jet(64);

for iextant=1:numel(IA) % foar all extant phenotypes
    clear Cl A
    
    for i=1:nbin
        inband=find(logX{iextant}>binedg(i) & logX{iextant}<=binedg(i+1));
        Cl{i}=logY{iextant}(inband)';
    end
    
    maxLengthCell=max(cellfun('size',Cl,2));    
    for i=1:length(Cl)
        for j=cellfun('size',Cl(i),2)+1:maxLengthCell
            Cl{i}(j)=NaN;   %zeropad the elements in each cell array with a length shorter than the maxlength
        end
    end
    A=cell2mat(Cl')'; %A is your matrix
    
    figure(7)
    subplot(5,6,iextant)
    
    hsdata=find(max(~isnan(A)));
    violin(A(:,hsdata));
    xval=binedg(hsdata);
    axis([0.5-hsdata(1) nbin+0.5-hsdata(1) -5 3])
    set(gca,'XTick',(1:nbin)-hsdata(1),'XTickLabel',cellstr(num2str((binedg)','%4.2f\n'))','XTickLabelRotation',90)    
    legend off    
    title(['#' num2str(iextant) ': ' char(964) '_{auto} = ' num2str(C(iextant,1)) ', ESD = ' num2str(C(iextant,2)) ' \mum'])
    
    figure(8)    
    iphen=find(ismember(traits,C(iextant,:),'rows'));
    
    plot(hsdata,nanmedian(A(:,hsdata)),'Color',cmap(clr(iextant),:),'LineWidth',2)
    hold on
    scatter(hsdata,nanmedian(A(:,hsdata)),2.*sqrt(C(iextant,2)),cmap(clr(iextant),:),'filled')
    text(hsdata(1),nanmedian(A(:,hsdata(1))),num2str(iextant),'Color',cmap(clr(iextant),:))
end




%%

% make grid of all lons and lats in model (size of EiE matrix)
% vi and vj are spatial indices TM 
LONgrid=repmat(OCEAN.lon(vi),[1 eco_pars.ntroph*eco_pars.nsize]); % 
LATgrid=repmat(OCEAN.lat(vj),[1 eco_pars.ntroph*eco_pars.nsize]); %
% load lat/lon data for extant populations
% lookup table of unique locations with biomass > threshold
H_lt = LATgrid(usedata);
H_ln = LONgrid(usedata);

% table of location indices for each extant population
extloc=[H_ln,H_lt];
[~,~,SC] = unique(extloc,'rows'); % index of unique locations
% SC maps unique locations to usedata index (expand B to match usedata)

% get binary genome
GENOME=full(cell2mat(matObj.GENOME(nyr,1)));
nbit=53; % log2(FlIntMax) = 53
str='';
for i=1:eco_pars.ngenes
    bin_gene = dec2bin(GENOME(usedata,i),nbit); % convert double value to binary string
    str=[str bin_gene];
end
binstr=str-'0';
ww=nthroot(2,10).^(0:eco_pars.ngenes-1);
w=reshape(repmat(ww,nbit,1),1,[]); % weight genes to have different mutation rates
w=w./mean(w);
binstr=binstr./w; % multiply bits by inverse mutation weights (more likely worth less)
D = squareform(pdist(binstr,'CityBlock')./size(binstr,2)); % genetic distance
% convert to years (empirical)
alpha =      0.01804;
beta  =       11.44;
g_dist     = beta.*D ./ (alpha .* (beta - D));



[xy,latsrt]=sortrows([vj(SC) vi(SC)]);

G=digraph(B); % shortestpath distance
t_dist = distances(G); % transport distance
t_dist = t_dist(SC,SC); % expand transport distance for all phenotypes

xg=[1 find(diff(xy(:,2)))' numel(usedata)];
xgl=reshape([xg;xg;NaN(1,length(xg))],1,length(xg)*3);
yg=[1 numel(usedata)];
ygl = repmat([yg NaN],1,length(xg));

figure(9)

subplot(211)
imagesc(log10(t_dist(latsrt,latsrt)));
set(gca,'XTick',[1 find(diff(xy(:,1)))' numel(usedata)],'XTicklabel',cellstr(num2str(Ye,'%2.0f')),...
        'YTick',[1 find(diff(xy(:,1)))' numel(usedata)],'YTicklabel',cellstr(num2str(Ye,'%2.0f')))
hold on
plot(xgl,ygl,'-','Color',[0.5 0.5 0.5])
plot(ygl,xgl,'-','Color',[0.5 0.5 0.5])
set(gca,'GridAlpha',1,'GridColor','k')
grid on
axis square

subplot(212)
imagesc(log10(g_dist(latsrt,latsrt)))
set(gca,'XTick',[1 find(diff(xy(:,1)))' numel(usedata)],'XTicklabel',cellstr(num2str(Ye,'%2.0f')),...
        'YTick',[1 find(diff(xy(:,1)))' numel(usedata)],'YTicklabel',cellstr(num2str(Ye,'%2.0f')))
hold on
plot(xgl,ygl,'-','Color',[0.5 0.5 0.5])
plot(ygl,xgl,'-','Color',[0.5 0.5 0.5])
set(gca,'GridAlpha',1,'GridColor','k')
grid on
axis square
caxis([2 4])


% CC is list of spatially sorted unique coordinates
% SA is index of first occurence of each location in [H_lt,H_ln]
% H_lt(SA) is spatially sorted latitudes (with unique associated longitudes)
% H_ln(SA) is spatially sorted longitudes (with unique associated latitudes)
% B(SC,SC) maps unique locations to order of usedata

















% N.B. str and binstr are already in usedata order

%% %%%%%%%%%%%%%%%%%%%%%%%%%%

nbin=15;
% define distance bins
binedg=linspace(-5,-1,nbin+1);
for i=1:nbin
    inband=find(log10(t_dist)>binedg(i) & log10(t_dist)<=binedg(i+1));
    Cl{i}=log10(g_dist(inband)');
end

maxLengthCell=max(cellfun('size',Cl,2));    
for i=1:length(Cl)
    for j=cellfun('size',Cl(i),2)+1:maxLengthCell
        Cl{i}(j)=NaN;   %zeropad the elements in each cell array with a length shorter than the maxlength
    end
end
A=cell2mat(Cl')'; %A is your matrix

figure(10)
hsdata=find(max(~isnan(A)));
violin(A(:,hsdata));
xval=binedg(hsdata);
% axis([0.5-hsdata(1) nbin+0.5-hsdata(1) -5 11])
set(gca,'XTick',(1:nbin)-hsdata(1),'XTickLabel',cellstr(num2str((binedg)','%4.2f\n'))','XTickLabelRotation',90)
legend off



    
    
%% 
% 
% 
% figure(7)
% clf
% 
% nbin=20;
% xh=linspace(2,11,nbin);
% yh=linspace(-4,2,nbin);
% 
% for iextant=1:28
%     
%     iphen=find(ismember(traits,C(iextant,:),'rows'));
%     
%     BIO=PHY0(:,iphen);
%     isgd=find(BIO>threshold);
%     BIO(BIO<threshold)=threshold;
%     
%    
%     DX=1./(B.*BIO);
%     G=digraph(DX,'OmitSelfLoops');
%     sp=distances(G);
%     sp=sp(isgd,isgd);
%     
%     spcstr=binstr(IC==iextant,:); % find genes of current phenotype
%     D = squareform(pdist(spcstr,'CityBlock')./size(spcstr,2));
%     
%     figure(2)
%     subplot(4,7,iextant)
%     N = hist3([log10(D(:)),log10(sp(:))], {yh xh});
%     imagesc(xh,yh,log10(N))
%     hold on
%     logD =log10(D(:));
%     logSP=log10(sp(:));
%     scatter(logSP,logD,2,'k.','MarkerEdgeAlpha',0.01,'MarkerFaceAlpha',0.01)
%     axis xy
%     ylabel('genetic distance')
%     xlabel('shortest path distance')
%     axis square
%     
%     isr=min([~isinf(logSP) ~isnan(logSP) ~isinf(logD) ~isnan(logD)],[],2);
%     [r,~] = corrcoef(logD(isr),logSP(isr));
%     [p,S] = polyfit(logD(isr),logSP(isr),1);
%     [y_fit,delta] = polyval(p,logD(isr),S); 
%     plot(y_fit,logD(isr),'r-','LineWidth',2)
%     plot(y_fit+2*delta,logD(isr),'m--',y_fit-2*delta,logD(isr),'m--')
%     title(r(2))
% end
% 
% 
% % isgd=find(PHY>1e-6);























