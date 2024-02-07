function Food_web_plot(output_dir,threshold)
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

% for nyr=i_timeslice
figure(1)
clf
colormap(jet)

data=zeros(961,i_timeslice(end));


data=full(cell2mat(matObj.PHY(:,1)'));
data=sum(data,1)';
data=reshape(data,[961,i_timeslice(end)]);
data(data==0)=NaN;
data=nthroot(data,2).*100;


i_year=ones([size(data,1) 1]).*i_timeslice;%round(gen_pars.save_output(i_timeslice(:)));
x=      repmat(1-eco_pars.trophic,1267,1);%.*ones(1,i_timeslice(end));
y=log10(repmat(eco_pars.V,1267,1));%.*ones(1,i_timeslice(end));

scatter3(i_year(:),x,y,data(:),repmat(rgb,1267,1),'filled')
hold on
scatter3(i_year(:),x,y,data(:).*0.8,'k','filled')


view(60,70)
ax = gca;
ax.BoxStyle = 'full';
xlim([0 i_timeslice(end)])

%%

data=full(cell2mat(matObj.PHY(:,1)'));
data=sum(data,1)';
data=reshape(data,[31,31,i_timeslice(end)]);

x=round(gen_pars.save_output(i_timeslice));
y=1-unique(eco_pars.trophic);
z=log10(unique(eco_pars.V));

[xx yy zz]=meshgrid(x,y,z);

v=permute(data,[1 3 2]);
%%
clf
for i=1200;%i_timeslice(end)
        
%     x1=x(1:i)+max(x)-x(i); % fix to t=max(x) (rhs)
    x1=x(1:i);             % fix to t=0 (lhs)
    v1=v(:,1:i,:);
    
    clf
    ax1 = axes;
    
    xmn=-3.9465;
    xmx=11.0535;
    patch([0 0 0 0]      ,[1 0 0 1],[xmn xmn xmx xmx],ones(1,3).*0.8)
    patch([0 0 x(i) x(i)],[1 0 0 1],ones(1,4).*xmn,ones(1,3).*0.8)
    patch([0 0 x(i) x(i)],[0 0 0 0],[xmx xmn xmn xmx],ones(1,3).*0.8)
    
    p = patch(isosurface(x1,y,z,v1,1e-2));
    set(gca,'YDir','Reverse')
    
    c1=ceil(mean(1-p.YData(1,:),1)'.*31);
    c1(c1==0)=1;
    c2=mean(p.ZData(1,:),1)';
    c2=round((c2-min(z))./(max(z)-min(z)).*30+1);
    ci=(c2-1).*31+c1;
    
    cscat=1:961;
    
    
    p.FaceVertexCData=ci;
    colormap(rgb)
    p.FaceColor = 'flat';
    p.EdgeColor = 'none';
    axis tight
    camlight
    lighting gouraud
    material METAL
    caxis([1 961])
    
    xlim([1 max(x)])
    ylim([0 1])
    zlim([xmn xmx])
%     set(gca,'XScale','Log')
    view(40,20);
    
    camproj('orthographic')
    camlight('headlight')
    
    hold on
    
    x1=squeeze(xx(:,1,:));
    y1=squeeze(yy(:,1,:));
    z1=squeeze(zz(:,1,:));
    xn=squeeze(xx(:,i,:));
    yn=squeeze(yy(:,i,:));
    zn=squeeze(zz(:,i,:));
    vn=100.*nthroot(squeeze(v(:,i,:)),2);
    vn(vn<=0)=NaN;
    
%     mesh(x1,y1,z1,'FaceAlpha',0,'EdgeColor','k')
    scatter3(x1(:),y1(:),z1(:),4,rgb,'filled')
    h(1)=scatter3(xn(:),yn(:),zn(:),vn(:)+sqrt(vn(:).*10),'k','filled');
    h(2)=scatter3(xn(:),yn(:),zn(:),vn(:),rgb,'filled');
    
    box on
    xlim([0 max(x)])
    
    xlabel('Time (years)')
    ylb=ylabel('Autotrophic coefficient','HorizontalAlignment','Left','Rotation',15.5,'Position',[3707.1    0.5   -5.1])
    zlabel('Log_{10}(Size)','Rotation',0,'HorizontalAlignment','Right')
    
    % move scatter plots to identical hodden axis on top of the first
    axHidden = axes('Visible','off','hittest','off'); % Invisible axes
    linkprop([ax1 axHidden],{'CameraPosition' 'XLim' 'YLim' 'ZLim' 'Position' 'YDir'}); % The axes should stay aligned
    set(h,'Parent',axHidden); % Put the text in the invisible Axes
    
    set(gcf,'color','w');
    drawnow
    fname=[output_dir '/Figures/Tree_' num2str(i,'%04i') '.png'];
%     fname=['~/GoogleDrive/MatrixMetacommunityModel/Figures/Tree_' num2str(i,'%04i') '.png'];
    export_fig(fname,'-nocrop','-r300')
end
whos
return
%%

subplot(221)

axesm ('mollweid', 'Frame', 'off', 'Grid', 'off');
data=PO4{1200};


%%

% define EiE matrix grid coordinates
ESD     = (6.*eco_pars.V./pi)'.^(1/3);
ESDgrid = repmat(ESD,[size(PHY,1) 1]);
TROgrid = 1-repmat(eco_pars.trophic',[size(PHY,1) 1]);
LATgrid=repmat(vj,[1 eco_pars.ntroph*eco_pars.nsize]);
LONgrid=repmat(vi,[1 eco_pars.ntroph*eco_pars.nsize]);

% threshold = sum(PHY,2)./100; % threshold defined relative to community biomass at each location
% threshold = 1e-3;

[i j] = find(PHY>threshold);      % Extract biomass distributions
usedata = sub2ind(size(PHY),i,j); % find populations above biomass threshold

disp([num2str(numel(usedata)) ' local populations above threshold biomass'])

% resample data based on weights
weits = full(PHY(usedata));
% plot rarefaction curve
% % % % % % % % % % % % % % figure(666)
% % % % % % % % % % % % % % set(gca,'XScale','Log','YScale','Log')
% % % % % % % % % % % % % % hold on
% % % % % % % % % % % % % % for c=round(10.^(1:0.1:8))
% % % % % % % % % % % % % %     [~,idx] = datasample(weits,c,'Weights',weits);
% % % % % % % % % % % % % %     n=numel(unique(idx));
% % % % % % % % % % % % % %     plot(c,n,'k.');
% % % % % % % % % % % % % %     drawnow
% % % % % % % % % % % % % % end
% % % % % % % % % % % % % % plot(xlim,[numel(usedata) numel(usedata)])


nsample = 1e4;
[~,idx] = datasample(weits,nsample,'Weights',weits);

% get traits and geographic coordinates of useable data
H_sz = ESDgrid(usedata(idx));
H_tr = TROgrid(usedata(idx));
H_lt = LATgrid(usedata(idx));
H_ln = LONgrid(usedata(idx));

% calculate mean square distances
f_msd = @(XI,XJ)( mean((XI-XJ).^2,2) ); 
msD = pdist(RGBX, @(Xi,Xj) f_msd(Xi,Xj)); % [length = nchoosek(numel(usedata),2) ]

timeD = msD./2./96;
tic;tree = linkage(timeD,'average');toc
% tic;leafOrder = optimalleaforder(tree,timeD);toc
leafOrder = [];
%% Plot dendrogram
% Polar Dendrogram
CThresh=850;
[H,nodes,outperm] = polardendrogram(tree,nsp,...
                     'ColorThreshold',CThresh);

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
coeff = pca(full(RGBX)');
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
nclstr=40;
weits = full(PHY(usedata));
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
% Geographically orientated food-web plot

figure(3)
cmap=colormap(gray(84));
cmap=cmap(21:end,:);
cmap(end,:)=[1 1 1].*0.75;
lns=jet(nclstr);

for nyr=i_timeslice(end)
    
    disp(['Processing ' num2str(nyr) ' of '  num2str(numel(i_timeslice))])
    
    
    clf
%     threshold=1e-5;
    [I,J,V] = find(PHY.*(PHY>threshold)); % find non-zero biomass
    IND = sub2ind(size(PHY),I,J);
    
    coeff = pca(full(RGB(IND,:))');
    CLR=coeff(:,1:3); % set rgb to 3 principle components of full rgb gene
    CLR = (CLR-min(CLR(:)))./range(CLR(:));
    
    field=zeros(1,934);
    for i=unique(I)' % for all locations...
        array{i}=[J(I==i) V(I==i) CLR(I==i,:) iclstrs(I==i)]; % get trait-index, biomass value and rgb triple
        field(i)=sum(V(I==i));      % and sum all biomass for location (for color plot
    end
    
    data=gen_fcns.v2f(field,vi,vj,vk); % convert biomass to 36x36 GEnIE matrix
    data=squeeze(data(8,:,:));
    data(isnan(data))=Inf;
    data(data<=0)=NaN;
    b=imagesc(0.5:35.5,0.5:35.5,log10(data));
    set(b,'AlphaData',~isnan(data)) % set NaNs to zero alpha
    axis xy
    colormap(cmap);
    caxis([-4 1])
    
    hold on
    Val=[];
    rgb=[];
    cls=[];
    x=[];
    y=[];
    for i=unique(I)'
        j=array{i}(:,1);
        Val=[Val;array{i}(:,2)];
        rgb=[rgb;array{i}(:,3:5)];
        cls=[cls;array{i}(:,6)];
        troph=1-eco_pars.trophic(j);
        psize=(log10(eco_pars.V(j))-min(log10(eco_pars.V)))/range(log10(eco_pars.V));
        x=[x;vi(i)+0.05+troph.*0.9];
        y=[y;vj(i)+0.05+psize.*0.9];
    end
    Val(Val<0)=NaN;
    lns
    scatter(x-1,y-1,nthroot(Val,3)*250,lns(cls,:),'filled')
    
    % Plot lines dividing global domains
    hold on
    plot([0:36;0:36],[zeros(1,37);ones(1,37).*36],'linestyle','-','color','k')
    plot([zeros(1,37);ones(1,37).*36],[0:36;0:36],'linestyle','-','color','k')
        
    
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    title(num2str(nyr))
    drawnow
    
%     fname=[output_dir '/Figures/Foodweb/Frame_' num2str(nyr,'%04i') '.png'];
%     export_fig(fname)
end

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
        clustrs{i,j}=squeeze(data(8,:,:));    % place in ntraitsxntraits cell array
        
        data=gen_fcns.v2f(pbio,vi,vj,vk); % convert to 36x36 GEnIE matrix
        biomass{i,j}=squeeze(data(8,:,:));    % place in ntraitsxntraits cell array
        
%         vector=zeros(934,1);
%         vector(x(y==ii))=nodes(y==ii);
%         data=gen_fcns.v2f(vector.*extant,vi,vj,vk); % convert to 36x36 GEnIE matrix
%         array{i,j}=squeeze(data(8,:,:));    % place in ntraitsxntraits cell array

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



