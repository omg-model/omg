function Cluster_data_3(output_dir,threshold)
nsp=0;

disp(['Processing output from ' output_dir])

% get matObj for matfile data
matObj = matfile([output_dir '/OUTPUT.mat']);
eco_pars=matObj.eco_pars;
gen_pars=matObj.gen_pars;

load([output_dir '/matrix_vars.mat'])
vi=v_index.i(Ib);
vj=v_index.j(Ib);
vk=v_index.k(Ib);

gen_fcns=general_functions;

%%


%%

PO4=matObj.PO4;
i_timeslice=find(~cellfun(@isempty,PO4))'; % find years with output data

iyear=round(gen_pars.save_output(i_timeslice(end)));
disp(['Processing output from year ' num2str(iyear)])
if ~exist([output_dir '/Figures/'],'dir')
    wd=cd(output_dir);
    !mkdir Figures
    cd(wd)
end
  
%%
% for nyr=i_timeslice
nyr=i_timeslice(end);

PHY=cell2mat(matObj.PHY(nyr,1));
RGB=cell2mat(matObj.RGB(nyr,1));

% define EiE matrix grid coordinates
ESD     = (6.*eco_pars.V./pi)'.^(1/3);
ESDgrid = repmat(ESD,[size(PHY,1) 1]);
TROgrid = 1-repmat(eco_pars.trophic',[size(PHY,1) 1]);
LATgrid=repmat(vj,[1 eco_pars.ntroph*eco_pars.nsize]);
LONgrid=repmat(vi,[1 eco_pars.ntroph*eco_pars.nsize]);

% threshold = sum(PHY,2)./100; % threshold defined relative to community biomass at each location
% threshold = 1e-3;

[i j k] = find(PHY>threshold);    % Extract biomass distributions
usedata = sub2ind(size(PHY),i,j); % find populations above biomass threshold

disp([num2str(numel(usedata)) ' local populations above threshold biomass'])

% get trait and geographic coordinates of useable data
H_sz = ESDgrid(usedata);
H_tr = TROgrid(usedata);
H_lt = LATgrid(usedata);
H_ln = LONgrid(usedata);

% extract RGB data for populations above threshold
clear X
X=RGB(usedata,:);

% calculate mean square distances
f_msd = @(XI,XJ)( mean((XI-XJ).^2,2) );
msD = pdist(X, @(Xi,Xj) f_msd(Xi,Xj)); % [length = nchoosek(numel(usedata),2) ]

timeD = msD./2./96;
tic;tree = linkage(timeD,'average');toc
% tic;leafOrder = optimalleaforder(tree,timeD);toc
leafOrder = [];
%% Plot dendrogram
% Polar Dendrogram
CThresh=850;
% [H,nodes,outperm] = polardendrogram(tree,nsp,...
%                      'ColorThreshold',CThresh);
[H,nodes,outperm] = polardendrogram(tree,nsp,...
                     'ColorThreshold',CThresh,...
                     'Reorder',leafOrder);

set(gcf,'Position',[53    50   924   932])

% reorder to match dendrogram   
sz=H_sz(outperm);
tr=H_tr(outperm);
lt=H_lt(outperm);
ln=H_ln(outperm);

% combine traits and locations in matrix
Z = [tr log10(sz) lt ln];
Z = padarray(Z,[1 1],'post'); 

r = linspace(1.1,1.4,size(Z,2));
theta = linspace(0.01*2*pi,0.99*2*pi,size(Z,1));

[THETA R] = ndgrid(theta,r);
[X Y]     = pol2cart(THETA,R);

rl=max(r);
cmap={'greenmag','parula','redblue','hsv'};
clabels={'trophic strategy','size','latitude','longitude'};

for i=1:size(Z,2)-1
    if i>1
        ax=axes;
    else
        ax=gca;
    end
    
    surf(X(:,i:i+1),Y(:,i:i+1),-ones(size(Z,1),2),Z(:,i:i+1))
    view(2)
    shading flat
    colormap(ax,cmap{i})
    caxis(minmax(Z(:,i)'))
    
    axis([-rl rl -rl rl]);

    axis off
    axis square
    
    axcl=axes;
    h=colorbar('South');
    h.TickLabels=[];
    h.Position=[0.09 0.12-0.025*i   0.25    0.02];
    h.Label.String=clabels{i};
    h.Label.HorizontalAlignment='left';    
    h.Label.VerticalAlignment='middle';
    h.Label.Position=[1.025 0.5 0];
    h.Label.FontSize=11;
    
    colormap(axcl,cmap{i});
    axis off
    
end

% plot rgb gene
[THETA R] = ndgrid(theta,[1 1+diff(r(1:2))]);
[X Y]     = pol2cart(THETA,R);
% rgb=full(RGB(usedata,1:3));
coeff = pca(full(RGB(usedata,:))');
rgb=coeff(:,1:3); % set rgb to 3 principle components of full rgb gene
for i=1:3
    rgb(:,i)=(rgb(:,i)-min(rgb(:,i)))./range(rgb(:,i)); % normalise
end
rgb=rgb(outperm,:); % sort to match terminal leaf node order
ax=axes;
Z=repmat((1:size(X,1))',[1 2]);
surf(X,Y,-ones(size(Z)),Z)
view(2)
shading flat
colormap(ax,rgb(:,1:3))
axis([-rl rl -rl rl]);

axis off
axis square


for i=[1 1+diff(r(1:2)) r]
    axes(ax)
    hold on
    hp=polar(linspace(0.01*2*pi,0.99*2*pi,100),ones(1,100).*i,'k-');
    hp.LineWidth=1;
end
hp=polar([0.01*2*pi 0.01*2*pi],[r(1) r(end)],'k-');hp.LineWidth=1;
hp=polar([0.99*2*pi 0.99*2*pi],[r(1) r(end)],'k-');hp.LineWidth=1;
hp=polar([0.01*2*pi 0.01*2*pi],[1 1+diff(r(1:2))],'k-');hp.LineWidth=1;
hp=polar([0.99*2*pi 0.99*2*pi],[1 1+diff(r(1:2))],'k-');hp.LineWidth=1;


% calculate branching frequency
% Transform cartesian tree to polar coordinates
XData = cell2mat(get(H,'XData'));
YData = cell2mat(get(H,'YData'));
[TData RData] = cart2pol(XData,YData);
[Ndiv,~] = histcounts(RData(:,2),100); % 'age' of branching points
% Ndiv     = Ndiv./cumsum(Ndiv); % normalise branching frequency by cumulative number of branchings


ax=axes;
uistack(ax,'bottom') % put axis for central histogram at the bottom
Zhist = repmat(Ndiv,[size(Z,1) 1]);
t = linspace(-theta(1),theta(1),size(Zhist,1));
r = linspace(0,1,size(Zhist,2));
[THETA R] = ndgrid(t,r);
[X Y]     = pol2cart(THETA,R);

surf(X,Y,-ones(size(Zhist)),log10(Zhist))
h=colorbar('South');
h.TickLabels=[];
h.Position=[0.09 0.12   0.25    0.02];
h.Label.String='branching frequency';
h.Label.HorizontalAlignment='left';
h.Label.VerticalAlignment='middle';
h.Label.Position=[4.05 0.5 0];
h.Label.FontSize=11;

axis([-rl rl -rl rl]);
hold on

colormap(ax,flipud(hot))
shading flat
axis off
axis square





% cluster colours
nclstr=35;
Tclstrs = cluster(tree,'maxclust',nclstr);
iclstrs = Tclstrs(outperm);

ax=axes;
uistack(ax,'bottom') % put axis for central histogram at the bottom
nClust = repmat(iclstrs,[1 500]);
r = linspace(0,1,size(nClust,2));
theta = linspace(0.01*2*pi,0.99*2*pi,size(nClust,1));
[THETA R] = ndgrid(theta,r);
[X Y]     = pol2cart(THETA,R);

RMin=min(RData,[],2);
TD=TData;
TD(TD<0)=TD(TD<0)+2.*pi;

for i=1:nclstr
    ii=find(iclstrs==i);   % get clockwise index of cluster
        Arc=minmax(theta(ii)); % find min and max theta
        
        jj=find(min(TD,[],2)>Arc(1) & max(TD,[],2)<Arc(2));
        
        kk=find(r<min(RMin(jj)));
        nClust(ii,kk)=NaN;
end
    
hs=surf(X,Y,-ones(size(nClust)).*0.5,nClust);
hs.FaceAlpha=0.2;
h=colorbar('South');
h.TickLabels=[];
h.Position=[0.09 0.12   0.25    0.02];
h.Label.String='branching frequency';
h.Label.HorizontalAlignment='left';
h.Label.VerticalAlignment='middle';
h.Label.Position=[4.05 0.5 0];
h.Label.FontSize=11;

axis([-rl rl -rl rl]);
hold on

colormap(ax,lines)
shading flat
axis off
axis square






fname=[output_dir '/Figures/OMG_Polar_' num2str(nyr,'%04i') '.png'];
export_fig(fname,'-r300','-transparent')

%%

eie=zeros(size(PHY));%.*NaN;
eie(usedata)=Tclstrs;

clear clustrs biomass

ii=0;
for i=1:eco_pars.nsize
    for j=1:eco_pars.ntroph
        ii=ii+1;                            % next row in PHY matrix
        
        tmp = eie(:,ii);
        pbio=PHY(:,ii);
        
        data=gen_fcns.v2f(tmp,vi,vj,vk); % convert to 36x36 GEnIE matrix
        clustrs{i,j}=squeeze(data(end,:,:));    % place in ntraitsxntraits cell array
        
        data=gen_fcns.v2f(pbio,vi,vj,vk); % convert to 36x36 GEnIE matrix
        biomass{i,j}=squeeze(data(end,:,:));    % place in ntraitsxntraits cell array
        
%         vector=zeros(934,1);
%         vector(x(y==ii))=nodes(y==ii);
%         data=gen_fcns.v2f(vector.*extant,vi,vj,vk); % convert to 36x36 GEnIE matrix
%         array{i,j}=squeeze(data(end,:,:));    % place in ntraitsxntraits cell array

    end
end
clustrs=cell2mat(clustrs);
clustrs(isnan(clustrs))=NaN;      % set land to inf
clustrs(clustrs<=0)=inf;          % set zeros to NaN
clustrs=clustrs+1;

biomass=cell2mat(biomass);
biomass(isnan(biomass))=inf;      % set land to inf
biomass(biomass<0)=NaN;          % set zeros to NaN

lins=lines;
lins(1,:)=[0.25 0.25 0.25];
lins(end,:)=[1 1 1];
jt=jet(nclstr+2);
jt(1,:)=[0.25 0.25 0.25];
jt(end,:)=[1 1 1];

figure(99)
imagesc(clustrs)
axis xy
colormap(lins)

figure(100)
imagesc(biomass)
axis xy
colormap(jt)

%%

clear clustrs biomass

ii=0;
for i=1:eco_pars.nsize
    for j=1:eco_pars.ntroph
        ii=ii+1;                            % next row in PHY matrix
        
        tmp = eie(:,ii);
        pbio=PHY(:,ii);
        
        data=gen_fcns.v2f(tmp,vi,vj,vk); % convert to 36x36 GEnIE matrix
        clustrs{i,j}=squeeze(data(end,:,:));    % place in ntraitsxntraits cell array
        
        data=gen_fcns.v2f(pbio,vi,vj,vk); % convert to 36x36 GEnIE matrix
        biomass{i,j}=squeeze(data(end,:,:));    % place in ntraitsxntraits cell array
        
%         vector=zeros(934,1);
%         vector(x(y==ii))=nodes(y==ii);
%         data=gen_fcns.v2f(vector.*extant,vi,vj,vk); % convert to 36x36 GEnIE matrix
%         array{i,j}=squeeze(data(end,:,:));    % place in ntraitsxntraits cell array

    end
end
clustrs=cell2mat(clustrs);
clustrs(isnan(clustrs))=NaN;      % set land to inf
clustrs(clustrs<=0)=inf;          % set zeros to NaN
clustrs=clustrs+1;

biomass=cell2mat(biomass);
biomass(isnan(biomass))=inf;      % set land to inf
biomass(biomass<0)=NaN;          % set zeros to NaN

lins=lines;
lins(1,:)=[0.25 0.25 0.25];
lins(end,:)=[1 1 1];
jt=jet(nclstr+2);
jt(1,:)=[0.25 0.25 0.25];
jt(end,:)=[1 1 1];

figure(99)
imagesc(clustrs)
axis xy
colormap(lins)

figure(100)
imagesc(biomass)
axis xy
colormap(jt)

%%
figure(3)

for nyr=i_timeslice(end)
    
    disp(['Processing ' num2str(nyr) ' of '  num2str(numel(i_timeslice))])
    
    
    clf
    [I,J,V] = find(PHY);
    
    
    field=zeros(1,934);
    for i=unique(I)'
        array{i}=[J(I==i) V(I==i)];
        field(i)=sum(V(I==i));
    end
    
    data=gen_fcns.v2f(field,vi,vj,vk);          % convert to 36x36 GEnIE matrix
    data=squeeze(data(end,:,:));
    data(isnan(data))=Inf;
    data(data<=0)=NaN;
    b=imagesc(0.5:35.5,0.5:35.5,log10(data));
    set(b,'AlphaData',~isnan(data))
    axis xy
    colormap(cmap);
    caxis([-4 1])
    
    hold on
    Val=[];
    x=[];
    y=[];
    for i=unique(I)'
        j=array{i}(:,1);
        Val=[Val;array{i}(:,2)];
        troph=1-eco_pars.trophic(j);
        psize=(log10(eco_pars.V(j))-min(log10(eco_pars.V)))/range(log10(eco_pars.V));
        x=[x;vi(i)+0.05+troph.*0.9];
        y=[y;vj(i)+0.05+psize.*0.9];
    end
    Val(Val<0)=NaN;
    scatter(x-1,y-1,nthroot(Val,3)*250,'k','filled')
    
    % Plot lines dividing global domains
    hold on
    plot([0:36;0:36],[zeros(1,37);ones(1,37).*36],'linestyle','-','color','k')
    plot([zeros(1,37);ones(1,37).*36],[0:36;0:36],'linestyle','-','color','k')
        
    
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    title(num2str(nyr))
    drawnow
    
    fname=[output_dir '/Figures/Foodweb/Frame_' num2str(nyr,'%04i') '.png'];
    export_fig(fname)
end



