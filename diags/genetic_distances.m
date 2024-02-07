function genetic_distances(output_dir,threshold,nsp,nclust)


if nclust>nsp & nsp~=0
    error('number of clusters must be <= number of terminal leaf nodes')
end
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


load('../TM_data/rwlla_annual/OCEAN.mat')


%%
X=OCEAN.lon_edges(3:end-2);
Y=OCEAN.lat_edges;
[xx yy]=meshgrid(X,Y);

x=OCEAN.lon(3:end-2);
y=OCEAN.lat;
[xlon ylat]=meshgrid(x,y);
xn=B'*xlon(:);
yn=B'*ylat(:);
u=xn-xlon(:);
v=yn-ylat(:);

% u=sign(u).*sqrt(abs(u));
% v=sign(v).*sqrt(abs(v));

% bi-directional graph of non-zero elements in transport matrix 
% Edges are weighted by the inverse values, such that stronger fluxes have higher weights
G=digraph(1./B,'OmitSelfLoops');
sp=distances(G);

figure(1)
clf

ip=round(linspace(1,252,16));
for i=1:16
    subplot(4,4,i)
    
%     vect=sp(:,ip(i));
%     vect=sp(ip(i),:)';
    vect=min([sp(:,ip(i)) sp(ip(i),:)'],[],2);
%     vect=mean([sp(:,ip(i)) sp(ip(i),:)'],2);
    
    map=reshape(vect,18,14);
    
    ax=axesm ('flatplrs', 'Frame', 'on', 'Grid', 'off','MapLonLimit',[-140 140]);
    pcolorm(Y,X,log10(map))
    axis tight
    axis off
%     caxis([2 4])
    
    
    % plot flow vectors
    hold on
    hq=quiverm(yn,xn,v,u,'k');
    plotm(yn(ip(i)),xn(ip(i)),'ro')
%     hq(1).LineWidth=1;
%     hq(2).LineWidth=1;
end

figure(2)
clf

subplot(211)
ax=axesm ('flatplrs', 'Frame', 'on', 'Grid', 'off','MapLonLimit',[-140 140]);
vect=median(sp,1);
map=reshape(vect,18,14);
pcolorm(Y,X,log10(map))
hold on
hq=quiverm(yn,xn,v,u,'k');
colorbar
% caxis([log10(2500) log10(1e5)])
axis tight
axis off
title('Median global dispersal time')

subplot(212)
ax=axesm ('flatplrs', 'Frame', 'on', 'Grid', 'off','MapLonLimit',[-140 140]);
vect=median(sp,2);
map=reshape(vect,18,14);
pcolorm(Y,X,log10(map))
hold on
hq=quiverm(yn,xn,v,u,'k');
colorbar
% caxis([log10(2500) log10(7e4)])
axis tight
axis off
title('Median global arrival time')



%%


clf
imagesc(log10(sp))
grid on
set(gca,'XTick',0:18:252,'XTickLabel',num2str((0:14)'),...
        'YTick',0:18:252,'YTickLabel',num2str((0:14)'));
view(-90,90)
%%
clear sp G DX isgd PHY ncl
PHY=cell2mat(matObj.PHY(nyr,1));
nplank=899; % 218,899
PHY=PHY(:,nplank);
isgd=find(PHY>1e-6);

DX=1./(B.*PHY);

G=digraph(DX,'OmitSelfLoops');
sp=distances(G);

x=1:14;
y=1:18;
[xx yy]=meshgrid(x,y);

ncl=10;

tree=linkage(sp(isgd,isgd));

figure(3)
clf
subplot(211)
[H, nodes, outperm] = dendrogram(tree,ncl);

subplot(212)
map=outperm(nodes);

scatter(xx(isgd),yy(isgd),50,map,'filled')
hold on
contour(reshape(log10(PHY),18,14),-3:-1,'k')
colormap(lhsdesign(ncl,3))
caxis([1 ncl])

%%
PHY=cell2mat(matObj.PHY(nyr,1));

traits=[1-eco_pars.trophic ESD'];
[C,IA,IC] = unique([H_tr,H_sz],'rows');

x=1:14;
y=1:18;
[xx yy]=meshgrid(x,y);

binstr=str-'0';
ww=nthroot(2,10).^(0:eco_pars.ngenes-1);
w=reshape(repmat(ww,nbit,1),1,[]); % weight genes to have different mutation rates
w=w./mean(w);
binstr=binstr./w; % multiply bits by inverse mutation weights (less likely worth more)

nbin=20;
xh=linspace(2,11,nbin);
yh=linspace(-4,2,nbin);

for iextant=1:28
    
    iphen=find(ismember(traits,C(iextant,:),'rows'));
    
    BIO=PHY(:,iphen);
    isgd=find(BIO>threshold);
    BIO(BIO<threshold)=threshold;
    
    figure(1)
    subplot(4,7,iextant)
    pcolor(log10(reshape(BIO,18,14)))
    caxis([log10(threshold) -1])
%     title(strcat(char(vpa(C(iextant,1),5)),'___',char(vpa(C(iextant,2),5))))
    
    
    DX=1./(B.*BIO);
    G=digraph(DX,'OmitSelfLoops');
    sp=distances(G);
    sp=sp(isgd,isgd);
    
    spcstr=binstr(IC==iextant,:); % find genes of current phenotype
    D = squareform(pdist(spcstr,'CityBlock')./size(spcstr,2));
    
    figure(2)
    subplot(4,7,iextant)
    N = hist3([log10(D(:)),log10(sp(:))], {yh xh});
    imagesc(xh,yh,log10(N))
    hold on
    logD =log10(D(:));
    logSP=log10(sp(:));
    scatter(logSP,logD,2,'k.','MarkerEdgeAlpha',0.01,'MarkerFaceAlpha',0.01)
    axis xy
    ylabel('genetic distance')
    xlabel('shortest path distance')
    axis square
    
    isr=min([~isinf(logSP) ~isnan(logSP) ~isinf(logD) ~isnan(logD)],[],2);
    [r,~] = corrcoef(logD(isr),logSP(isr));
    [p,S] = polyfit(logD(isr),logSP(isr),1);
    [y_fit,delta] = polyval(p,logD(isr),S); 
    plot(y_fit,logD(isr),'r-','LineWidth',2)
    plot(y_fit+2*delta,logD(isr),'m--',y_fit-2*delta,logD(isr),'m--')
    title(r(2))
end


% isgd=find(PHY>1e-6);























