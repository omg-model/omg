function Cluster_data(output_dir,nsp)


disp(['Processing output from ' output_dir])

% load saved data
load([output_dir '/OUTPUT.mat'])

load([output_dir '/matrix_vars.mat'])
vi=v_index.i(Ib);
vj=v_index.j(Ib);
vk=v_index.k(Ib);

gen_fcns=general_functions;

%%


%% Prepare figure axes
figure(99)
set(gcf,'units','normalized','outerposition',[0 0.1 0.9 0.9])
% drawnow


    
i_timeslice=find(~cellfun(@isempty,PO4))'; % find years with output data

% color axis limits
nphy=eco_pars.nsize*eco_pars.ntroph;
%%

disp(['Processing output from ' num2str(numel(i_timeslice)) ' time-slices'])
if ~exist([output_dir '/Figures/'],'dir')
    wd=cd(output_dir);
    !mkdir Figures
    cd(wd)
end
  

%%
% for nyr=i_timeslice
nyr=i_timeslice(end);

locsum=sum(PHY{nyr},2); % community biomass at each location
threshold = locsum./100; % threshold defined relative to community biomass at each location
% threshold=1e-2;

[x y z] = find(PHY{nyr}>threshold);    % Extract biomass distributions
usedata = sub2ind(size(PHY{nyr}),x,y); % find populations above biomass threshold

% extract RGB data for populations above threshold
clear X
X=RGB{nyr}(usedata,:);

% calculate mean square distances
f_msd = @(XI,XJ)( mean((XI-XJ).^2,2) );
msD = pdist(X, @(Xi,Xj) f_msd(Xi,Xj)); % [length = nchoosek(numel(usedata),2) ]

timeD = msD./2./96;
<<<<<<< HEAD
tic;tree = linkage(timeD,'average');toc
% tic;leafOrder = optimalleaforder(tree,timeD);toc
leafOrder = [];
%% Plot dendrogram
=======
tree = linkage(timeD,'average');
% leafOrder = optimalleaforder(tree,eucD)
>>>>>>> parent of 8f01a61... new diagnostics

%
subplot(131)
% plot dendrogram
[H,nodes,outperm] = dendrogram(tree,nsp,...
                               'Orientation','left',...
                               'ColorThreshold',0);
% [~,nodes,outperm] = dendrogram(tree,nsp,'Orientation','left','reorder',leafOrder);
box on
xlabel('Estimated time before end of simulation (years)')
nd = size(RGB{1},2); % number of rgb tracers
title(['Phylogenetic tree based on mean squared distance after ' num2str(gen_pars.save_output(nyr)+0.5) ' years in ' num2str(nd) ' dimensions'])
xl=xlim; % get x axis limits

% Generate colormap for dendrogram leaf lbels
cmap=colormap(jet(nsp+1));
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
     scatter(min(xl),find(outperm==i),100,nclr(i,:),'filled');
%     scatter(min(xl),i,100,cmap(outperm(i),:),'filled');
end
set(gca,'YTickLabel','');
% get colors of dendrogram lines
for i=1:size(H,1)
    for j=H(i).YData(H(i).XData==0)
        clr(j,:)=H(i).Color;
    end
end
<<<<<<< HEAD
set(H,'LineWidth',0.75)
%%
[Ndiv,EDGES] = histcounts(AllX(:,2),100);
=======
>>>>>>> parent of 8f01a61... new diagnostics

% plot cluster diagram in 6 dimensions
figure(999),scatter3(X(:,4),X(:,5),X(:,6),25,rgb,'filled')
clear rgb
%%
% scatter3(X(:,1),X(:,2),X(:,3),25,c,'filled')
figure(99)
disp(['Processing ' num2str(nyr) ' of '  num2str(numel(i_timeslice)) ', with ' num2str(nnz(PHY{nyr})) ' non-zero local populations'])

clear vector data array matrix biomass

for i=1:nd
    newrgb{i}=reshape(RGB{nyr}(:,i),numel(vi),eco_pars.jpmax);
end

ii=0;
for i=1:eco_pars.nsize
    for j=1:eco_pars.ntroph
        ii=ii+1;                            % next row in PHY matrix
        
        pbio = PHY{nyr}(:,ii);
        extant = pbio>threshold;
        
        data=gen_fcns.v2f(pbio,vi,vj,vk); % convert to 36x36 GEnIE matrix
        biomass{i,j}=squeeze(data(8,:,:));    % place in ntraitsxntraits cell array
        
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
end
matrix=cell2mat(array);         % catenate cell-array to (ntraitsx36)x(ntraitsx36) matrix
matrix(isnan(matrix))=inf;      % set land to inf
matrix(matrix<=0)=NaN;          % set zeros to NaN
biomass=cell2mat(biomass);
biomass(isnan(biomass))=inf;      % set land to inf
biomass(biomass<=0)=NaN;          % set zeros to NaN
%%
sz=repmat((1:eco_pars.nsize)',[1,eco_pars.nsize] );
sz=kron(sz,ones(36,36));
tr=repmat( 1:eco_pars.ntroph ,[eco_pars.ntroph,1]);
tr=kron(tr,ones(36,36));

lon=repmat(1:36,[36,1] );
lon=repmat(lon,[eco_pars.ntroph,eco_pars.nsize]);
lat=repmat((1:36)',[1,36] );
lat=repmat(lat,[eco_pars.ntroph,eco_pars.nsize]);

<<<<<<< HEAD
if nsp==0;
    nsp=size(AllX,1);
end

=======
>>>>>>> parent of 8f01a61... new diagnostics
for i=1:nsp
    ii = find(matrix==i);
    jj = find(outperm==i);
    H_sz(jj,:) = accumarray(sz(ii), biomass(ii),[eco_pars.nsize ,1])';
    H_tr(jj,:) = accumarray(tr(ii), biomass(ii),[eco_pars.ntroph,1])';
    H_lon(jj,:) = accumarray(lon(ii), biomass(ii),[36,1])';
    H_lat(jj,:) = accumarray(lat(ii), biomass(ii),[36 ,1])';
end

% H_lon=H_lon./sum(H_lon,1);
% H_lat=H_lat./sum(H_lat,1);


subplot(1,36,[13:18])
cla
imagesc(log10(H_tr))
axis xy
caxis([-3 1])
set(gca,'XTick',1:eco_pars.ntroph);
set(gca,'XTickLabel',num2str(unique(eco_pars.trophic),2),'XTickLabelRotation',90);
set(gca,'YTick',[]);
% set(gca,'YTickLabel','');
xlabel('Plankton trophic strategy (A to H)')
ylabel('''Species'' (clusters in dendrogram)')
title('Log_{10} biomass by trophic strategy and ''species''')
grid on
hold on;
for i=1:nsp;  plot(xlim,[i i],'Color',clr(i,:)); end

subplot(1,36,[19:24])
cla
imagesc(log10(H_sz))
axis xy
caxis([-3 1])
set(gca,'XTick',1:eco_pars.nsize);
set(gca,'XTickLabel',num2str(eco_pars.unq_ESD,2),'XTickLabelRotation',90);
set(gca,'YTick',[]);
% set(gca,'YTickLabel','');
xlabel('Plankton ESD (\mum)')
title('Log_{10} biomass by size class and ''species''')
grid on
hold on;
for i=1:nsp;  plot(xlim,[i i],'Color',clr(i,:)); end

subplot(1,36,[25:30])
cla
imagesc(log10(H_lon))
axis xy
caxis([-3 1])
set(gca,'XTick',1:36);
% set(gca,'XTickLabel',num2str(unique(eco_pars.trophic),2));
set(gca,'XTickLabelRotation',90);
set(gca,'YTick',[]);
% set(gca,'YTickLabel','');
xlabel('Longitude band')
title('Log_{10} biomass by longitude and ''species''')
grid on
hold on;
for i=1:nsp;  plot(xlim,[i i],'Color',clr(i,:)); end

subplot(1,36,[31:36])
cla
imagesc(log10(H_lat))
axis xy
caxis([-3 1])
set(gca,'XTick',1:36);
% set(gca,'XTickLabel',num2str(unique(eco_pars.trophic),2));
set(gca,'XTickLabelRotation',90);
set(gca,'YTick',[]);
% set(gca,'YTickLabel','');
xlabel('Latitude band')
title('Log_{10} biomass by latitude and ''species''')
grid on
hold on;
for i=1:nsp;  plot(xlim,[i i],'Color',clr(i,:)); end

colormap(flipud(hot))

%%
<<<<<<< HEAD
% Polar Dendrogram
CThresh=850;
[H,nodes,outperm] = polardendrogram(tree,nsp,...
                     'ColorThreshold',CThresh,...
                     'reorder',leafOrder);

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




=======

drawnow
    
%     fname=[output_dir '/Figures/Species_' num2str(nyr,'%04i') '.png'];
%     export_fig(fname)

    
>>>>>>> parent of 8f01a61... new diagnostics

% wd=cd([output_dir '/Figures/']);
% !ffmpeg -framerate 14 -i Figures/Frame_%04d.png  -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" evolution.mp4
% cd(wd)
disp(' ')
disp('Example ffmpeg command to be executed in ''Figures'' directory ...')
disp('ffmpeg -framerate 14 -i Frame_%04d.png  -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" evolution.mp4')
disp(' ')


%%













