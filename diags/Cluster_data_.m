function Cluster_data(output_dir,nsp)


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

ESD     = (6.*eco_pars.V./pi)'.^(1/3);
ESDgrid = repmat(ESD,[size(PHY,1) 1]);

% threshold = sum(PHY,2)./100; % threshold defined relative to community biomass at each location
threshold = 1e-3;


[x y z] = find(PHY>threshold);    % Extract biomass distributions
usedata = sub2ind(size(PHY),x,y); % find populations above biomass threshold

disp([num2str(numel(usedata)) ' local populations above threshold biomass'])

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

figure(99)
clf
set(gcf,'units','normalized','outerposition',[0.0219 0 0.4176 0.9840])
% drawnow

% color axis limits
nphy=eco_pars.nsize*eco_pars.ntroph;

subplot(1,36,[1:11.9])
% plot dendrogram
CThresh=1000;
[H,nodes,outperm] = dendrogram(tree,nsp,...
                               'Orientation','left',...
                               'ColorThreshold',CThresh,...
                               'reorder',leafOrder);
% [~,nodes,outperm] = dendrogram(tree,nsp,'Orientation','left','reorder',leafOrder);
box on
xlabel('Estimated time before end of simulation (years)')
nd = size(RGB,2); % number of rgb tracers
title({'Phylogenetic tree based on mean squared distance';['after ' num2str(gen_pars.save_output(nyr)+0.5) ' years in ' num2str(nd) ' dimensions']})
xl=xlim; % get x axis limits

% Generate colormap for dendrogram leaf lbels
cnsp=nsp+1;
cmap=colormap(jet(cnsp));
cmap(end,:)=[1 1 1].*0.25;
% get colours from first three rgb tracers
rgb=full(X(:,1:3)); 
rgbl = max(abs(rgb));
rgbl(rgbl==0)=1;
rgb=(rgb+rgbl)./(2.*rgbl);
% scatter plot colours at right-hand edge of dendrogram (i.e. x=xl(1))
hold on
for i=unique(nodes)'
    ii=find(nodes==i);
    nclr(i,:)=full(median(rgb(ii,:),1));    
    scatter(min(xl)-range(xl)./100,find(outperm==i),25,nclr(i,:),'filled');
%     scatter(min(xl),i,100,cmap(outperm(i),:),'filled');
end
set(gca,'YTick',[]);
set(gca,'YTickLabel','');
xl=xlim;
xlim([min(xl)-range(xl)./100 xl(2)]);

XData = cell2mat(get(H,'XData'));
AllX=XData;
YData = cell2mat(get(H,'YData'));
% only below C threshold
ii=find(XData(:,2)<=CThresh);
XData = XData(ii,:);
YData = YData(ii,:);
% sort by Y instead of X
[~,isrt] = sort(min(YData,[],2));
YData = YData(isrt,:);
XData = XData(isrt,:);

lineColours = cell2mat(get(H(isrt),'Color'));

[unqclr,ia,ic]=unique(lineColours,'rows');
[~,~,cind]=unique(ic,'stable');

% cmap=[rem(1:numel(ia),2)' zeros(numel(ia),1) ~rem(1:numel(ia),2)'];
% cmap=prism(numel(ia));
cmap=lines(numel(ia));

for i=ii'    
    H(isrt(i)).Color = cmap(cind(i),:);
end
set(H,'LineWidth',0.75)
%%
[Ndiv,EDGES] = histcounts(AllX(:,2),100);

%%
  
% plot cluster diagram in 6 dimensions
figure(999)
scatter3(X(:,4),X(:,5),X(:,6),25,rgb,'filled')

%%
% scatter3(X(:,1),X(:,2),X(:,3),25,c,'filled')
figure(99)
disp(['Processing ' num2str(nyr) ' of '  num2str(numel(i_timeslice)) ', with ' num2str(nnz(PHY)) ' non-zero local populations'])

clear rgb vector data array matrix biomass

for i=1:nd
    newrgb{i}=reshape(RGB(:,i),numel(vi),eco_pars.jpmax);
end

ii=0;
for i=1:eco_pars.nsize
    for j=1:eco_pars.ntroph
        ii=ii+1;                            % next row in PHY matrix
        
        % get plankton biomass distribution
        pbio = PHY(:,ii);
        extant = pbio>threshold;
        data=gen_fcns.v2f(pbio,vi,vj,vk);  % convert to 36x36 GEnIE matrix
        biomass{i,j}=squeeze(data(8,:,:)); % place in ntraitsxntraits cell array
        
        % get terminal leaf node index
        vector=zeros(934,1);
        vector(x(y==ii))=nodes(y==ii);
        data=gen_fcns.v2f(vector.*extant,vi,vj,vk); % convert to 36x36 GEnIE matrix
        array{i,j}=squeeze(data(8,:,:));    % place in ntraitsxntraits cell array
        
        for irgb=1:3
            vector=newrgb{irgb}(:,ii);             % extract vector for PHY ii
            data=gen_fcns.v2f(vector.*extant,vi,vj,vk); % convert to 36x36 GEnIE matrix
            rgb{i,j,irgb}=squeeze(data(8,:,:));    % place in ntraitsxntraits cell array
        end
    end
e
% array of terminal leaf node indices
matrix=cell2mat(array);         % catenate cell-array to (ntraitsx36)x(ntraitsx36) matrix
matrix(isnan(matrix))=inf;      % set land to inf
matrix(matrix<=0)=NaN;          % set zeros to NaN
% array of biomass
biomass=cell2mat(biomass);
biomass(isnan(biomass))=inf;      % set land to inf
biomass(biomass<=0)=NaN;          % set zeros to NaN
%%

clear nclr clr H_lat H_lon H_sz H_tr
for i=1:size(H,1)
    for j=H(i).YData(H(i).XData==0)
        clr(j,:)=H(i).Color;
    end
end

sz=repmat((1:eco_pars.nsize)',[1,eco_pars.nsize] );
sz=kron(sz,ones(36,36));
tr=repmat( 1:eco_pars.ntroph ,[eco_pars.ntroph,1]);
tr=kron(tr,ones(36,36));

lon=repmat(1:36,[36,1] );
lon=repmat(lon,[eco_pars.ntroph,eco_pars.nsize]);
lat=repmat((1:36)',[1,36] );
lat=repmat(lat,[eco_pars.ntroph,eco_pars.nsize]);

if nsp==0;
    nsp=size(AllX,1);
end

for i=1:nsp
    ii = find(matrix==i);
    jj = find(outperm==i);
    H_sz(jj,:)  = accumarray(sz(ii), biomass(ii),[eco_pars.nsize ,1])';
    H_tr(jj,:)  = accumarray(tr(ii), biomass(ii),[eco_pars.ntroph,1])';
    H_lon(jj,:) = accumarray(lon(ii), biomass(ii),[36,1])';
    H_lat(jj,:) = accumarray(lat(ii), biomass(ii),[36 ,1])';
end

% H_lon=H_lon./sum(H_lon,1);
% H_lat=H_lat./sum(H_lat,1);


subplot(1,36,[13:18])
cla
imagesc(log10(H_sz))
axis xy
caxis([-4 1])
ax=gca;
ax.XTick=1:eco_pars.nsize;
labels=num2str(eco_pars.unq_ESD,1);
labels(2:2:end,:) = NaN;
set(gca,'XTickLabel',labels,'XTickLabelRotation',90);
set(gca,'YTick',[]);
% set(gca,'YTickLabel','');
hx=xlabel('Plankton ESD (\mum)');
hx.Position(2)=-50;
title('Size classes')
grid on
hold on;
for i=1:nsp;  plot(xlim,[i i],'Color',[clr(i,:) 0.05],'LineWidth',2); end

subplot(1,36,[19:24])
cla
imagesc(log10(H_tr))
axis xy
caxis([-4 1])
set(gca,'XTick',1:eco_pars.ntroph);
labels=num2str(unique(eco_pars.trophic),'%4.2f');
labels(2:2:end,:) = NaN;
set(gca,'XTickLabel',labels,'XTickLabelRotation',90);
set(gca,'YTick',[]);
% set(gca,'YTickLabel','');
hx=xlabel('Plankton trophic strategy (A to H)');
hx.Position(2)=-50;
title({'Trophic strategy'})
grid on
hold on;
for i=1:nsp;  plot(xlim,[i i],'Color',[clr(i,:) 0.05],'LineWidth',2); end

subplot(1,36,[25:30])
cla
imagesc(log10(H_lat))
axis xy
caxis([-4 1])
set(gca,'XTick',1:36);
labels=num2str((1:36)',2);
labels(2:2:end,:) = NaN;
set(gca,'XTickLabel',labels,'XTickLabelRotation',90);
set(gca,'YTick',[]);
% set(gca,'YTickLabel','');
hx=xlabel('Latitude band ');
hx.Position(2)=-50;
title('Latitude')
grid on
hold on;
for i=1:nsp;  plot(xlim,[i i],'Color',[clr(i,:) 0.05],'LineWidth',2); end

subplot(1,36,[31:36])
cla
imagesc(log10(H_lon))
axis xy
caxis([-4 1])
set(gca,'XTick',1:36);
labels=num2str((1:36)',2);
labels(2:2:end,:) = NaN;
set(gca,'XTickLabel',labels,'XTickLabelRotation',90);
set(gca,'YTick',[]);
% set(gca,'YTickLabel','');
hx=xlabel('Longitude band');
hx.Position(2)=-50;
ylabel('''Species'' (clusters in dendrogram)')
set(gca,'YAxisLocation','right');
title('Longitude')
grid on
hold on;
for i=1:nsp;  plot(xlim,[i i],'Color',[clr(i,:) 0.05],'LineWidth',2); end

colormap(flipud(hot(64)))

set(gcf,'color','w');
fname=[output_dir '/Figures/Bioinformatics_' num2str(nyr,'%04i') '.png'];
export_fig(fname,'-m2')


%%
% Polar Dendrogram
CThresh=850;
[H,nodes,outperm] = polardendrogram(tree,nsp,...
                     'ColorThreshold',CThresh,...
                     'reorder',leafOrder);


for i=1:nsp
    ii = find(matrix==i);
    jj = find(outperm==i);
    H_sz(jj,:)  = accumarray(sz(ii), biomass(ii),[eco_pars.nsize ,1])';
    H_tr(jj,:)  = accumarray(tr(ii), biomass(ii),[eco_pars.ntroph,1])';
    H_lon(jj,:) = accumarray(lon(ii), biomass(ii),[36,1])';
    H_lat(jj,:) = accumarray(lat(ii), biomass(ii),[36 ,1])';
end

sz=log10(eco_pars.unq_ESD);
tr=unique(eco_pars.trophic);
lt=linspace(0,1,36)';
ln=linspace(0,1,36)';
wsz=H_sz./sum(H_sz,2); % get size weights
wtr=H_tr./sum(H_tr,2);  % get trophic weights
wlt=H_lat./sum(H_lat,2);  % get latitude weights
wln=H_lon./sum(H_lon,2);  % get longitude weights
mnsz=wsz*sz; % weighted mean size
mnsz=(mnsz-min(mnsz))./range(mnsz);
mntr=wtr*tr; % weighted mean trophic
mnlt=wlt*lt; % weighted mean latitude
mnln=wln*ln; % weighted mean longitude


Z = [mntr mnsz mnlt mnln];
Z = padarray(Z,[1 1],'post');

r = linspace(1.02,1.3,size(Z,2));
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
    
    axis([-rl rl -rl rl]);

    axis off
    axis square
    
    axcl=axes;
    h=colorbar('South');
    h.TickLabels=[];
    h.Position=[0.105 0.125-0.025*i   0.25    0.02];
    h.Label.String=clabels{i};
    h.Label.HorizontalAlignment='left';    
    h.Label.VerticalAlignment='middle';
    h.Label.Position=[1.025 0.5 0];
    h.Label.FontSize=11;
    
    colormap(axcl,cmap{i});
    axis off
    
end

for i=[1 r]
    axes(ax)
    hold on
    hp=polar(linspace(0.01*2*pi,0.99*2*pi,100),ones(1,100).*i,'k-');
    hp.LineWidth=1;
end
hp=polar([0.01*2*pi 0.01*2*pi],[r(1) r(end)],'k-');
hp.LineWidth=1;
hp=polar([0.99*2*pi 0.99*2*pi],[r(1) r(end)],'k-');
hp.LineWidth=1;


ax=axes;
uistack(ax,'bottom')
Zhist = repmat(fliplr(Ndiv),[size(Z,1) 1]);
r = linspace(0,1,size(Zhist,2));
[THETA R] = ndgrid(theta,r);
[X Y]     = pol2cart(THETA,R);

surf(X,Y,-ones(size(Zhist)),log10(Zhist))
h=colorbar('South');
h.TickLabels=[];
h.Position=[0.105 0.125-0.025*5   0.25    0.02];
h.Label.String='branching density';
h.Label.HorizontalAlignment='left';
h.Label.VerticalAlignment='middle';
h.Label.Position=[4.055 0.5 0];
h.Label.FontSize=11;

axis([-rl rl -rl rl]);
hold on

colormap(ax,flipud(hot))
shading flat
axis off
axis square



set(gcf,'Position',[53    50   924   932])
fname=[output_dir '/Figures/OMG_Polar_' num2str(nyr,'%04i') '.png'];
export_fig(fname,'-r300','-transparent')





% wd=cd([output_dir '/Figures/']);
% !ffmpeg -framerate 14 -i Figures/Frame_%04d.png  -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" evolution.mp4
% cd(wd)
drawnow


%%













