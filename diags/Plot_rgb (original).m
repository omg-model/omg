function Plot_saved(output_dir)

disp(['Processing output from ' output_dir])

% load saved data
load([output_dir '/OUTPUT.mat'])

load([output_dir '/matrix_vars.mat'])
vi=v_index.i(Ib);
vj=v_index.j(Ib);
vk=v_index.k(Ib);

gen_fcns=general_functions;

%% Prepare figure axes
figure(1)
set(gcf,'units','normalized','outerposition',[0 0 1 1])
drawnow
gap = [0.01 0.005];
marg_h = 0.030;
marg_w = 0.02;

% color axis limits
maxPO4=0.1;
minPHY=0.0001;
maxPHY=1;


nphy=eco_pars.nsize*eco_pars.ntroph;
%%
    
i_timeslice=find(~cellfun(@isempty,PO4))'; % find years with output data

disp(['Processing output from ' num2str(numel(i_timeslice)) ' time-slices'])
if ~exist([output_dir '/Figures/'],'dir')
    wd=cd(output_dir);
    !mkdir Figures
    cd(wd)
end
%%   
% for nyr=i_timeslice
for nyr=fliplr(i_timeslice)%(end)

    PHY=cell2mat(matObj.PHY(nyr,1));
    GENOME=full(cell2mat(matObj.GENOME(nyr,1)));
    spy(GENO)
    
    GENOME=reshape(GENOME,[numel(Ib),nphy,size(GENOME,2)]);
%     for i=1:eco_pars.nsize
%         for j=1:eco_pars.ntroph
%             str='';
%             for i=1:eco_pars.ngenes
%                 
%                 str=[str dec2bin(GENOME(i,j,1),64)];

%     locsum=sum(PHY,2); % community biomass at each location
%     threshold = locsum./1000; % threshold defined relative to community biomass at each location
%     threshold=gen_fcns.v2f(threshold,vi,vj,vk); % convert to 36x36 GEnIE matrix
%     threshold=squeeze(threshold(8,:,:));
    threshold=1e-9;

    
    disp(['Processing ' num2str(nyr) ' of '  num2str(numel(i_timeslice)) ', with ' num2str(nnz(PHY{nyr})) ' non-zero local populations'])
    
    clear vector data array matrix
    
    rgb=full(RGB{nyr});
    rgb=reshape(rgb,[numel(Ib),nphy,size(rgb,2)]);
    
    ii=0;
    for i=1:eco_pars.nsize
        for j=1:eco_pars.ntroph
            ii=ii+1;                            % next row in PHY matrix
            
            vector=PHY{nyr}(:,ii);             % extract vector for PHY ii
            data=gen_fcns.v2f(vector,vi,vj,vk); % convert to 36x36 GEnIE matrix
            Pdata=squeeze(data(8,:,:));    % 
            arrextant{i,j}=Pdata>1e-12;
            
            for k=1:3
                vector=rgb(:,ii,k);              % extract vector for PHY ii
                data=gen_fcns.v2f(vector,vi,vj,vk); % convert to 36x36 GEnIE matrix
                array{i,j,k}=squeeze(data(8,:,:));    % place in ntraitsxntraits cell array
            end
        end
    end
    extant=cell2mat(arrextant);
    matrix=cell2mat(array);         % catenate cell-array to (ntraitsx36)x(ntraitsx36) matrix
    
    rgblim=max(abs(reshape(matrix,[],3)));
    for j=1:3
        matrix(:,:,j)=(matrix(:,:,j) .* extant+rgblim(j))./(2*rgblim(j));
    end
    
    
    image(matrix)
    axis xy
    
    % Plot lines dividing global domains
    hold on
    plot([0:36:eco_pars.ntroph*36;   0:36:eco_pars.ntroph*36                  ]+0.5,...
         [zeros(1,eco_pars.ntroph+1);ones(1,eco_pars.ntroph+1).*size(matrix,1)],'linestyle','-','color','k')
    plot([zeros(1,eco_pars.nsize+1);ones(1,eco_pars.nsize+1).*size(matrix,2)],...
         [0:36:eco_pars.nsize*36;   0:36:eco_pars.nsize*36                  ]+0.5,'linestyle','-','color','k')
    
    title(gen_pars.save_output(nyr))
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    
    
    
    drawnow
    
    fname=[output_dir '/Figures/RGB_' num2str(nyr,'%04i') '.png'];
    export_fig(fname,'-m2')
end
    

% wd=cd([output_dir '/Figures/']);
% !ffmpeg -framerate 14 -i Figures/RGB_%04d.png  -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" RGB_gene.mp4
% cd(wd)
disp(' ')
disp('Example ffmpeg command to be executed in ''Figures'' directory ...')
disp('ffmpeg -framerate 14 -i Frame_%04d.png  -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" evolution.mp4')
disp(' ')


return
%%
cmap=jet(size(i_timeslice,2));
clf
for nyr=i_timeslice%(end)
    
    disp(['Processing ' num2str(nyr) ' of '  num2str(numel(i_timeslice)) ', with ' num2str(nnz(PHY{nyr})) ' non-zero local populations'])
    
    clear vector data array matrix
    
    ii=0;
    for i=1:eco_pars.nsize
        for j=1:eco_pars.ntroph
            ii=ii+1;                            % next row in PHY matrix
            
            vector=PHY{nyr}(:,ii);             % extract vector for PHY ii
            data=gen_fcns.v2f(vector,vi,vj,vk); % convert to 36x36 GEnIE matrix
            Pdata=squeeze(data(8,:,:));    % 
            arrextant{i,j}=Pdata>1e-12;
            
            vector=RGB0001{nyr}(:,ii);             % extract vector for PHY ii
            data=gen_fcns.v2f(vector,vi,vj,vk); % convert to 36x36 GEnIE matrix
            array{i,j,1}=squeeze(data(8,:,:));    % place in ntraitsxntraits cell array
            
            vector=RGB0002{nyr}(:,ii);             % extract vector for PHY ii
            data=gen_fcns.v2f(vector,vi,vj,vk); % convert to 36x36 GEnIE matrix
            array{i,j,2}=squeeze(data(8,:,:));    % place in ntraitsxntraits cell array
            
            vector=RGB0003{nyr}(:,ii);             % extract vector for PHY ii
            data=gen_fcns.v2f(vector,vi,vj,vk); % convert to 36x36 GEnIE matrix
            array{i,j,3}=squeeze(data(8,:,:));    % place in ntraitsxntraits cell array
            
            vector=RGB0004{nyr}(:,ii);             % extract vector for PHY ii
            data=gen_fcns.v2f(vector,vi,vj,vk); % convert to 36x36 GEnIE matrix
            array{i,j,4}=squeeze(data(8,:,:));    % place in ntraitsxntraits cell array
            
            vector=RGB0005{nyr}(:,ii);             % extract vector for PHY ii
            data=gen_fcns.v2f(vector,vi,vj,vk); % convert to 36x36 GEnIE matrix
            array{i,j,5}=squeeze(data(8,:,:));    % place in ntraitsxntraits cell array
            
            vector=RGB0006{nyr}(:,ii);             % extract vector for PHY ii
            data=gen_fcns.v2f(vector,vi,vj,vk); % convert to 36x36 GEnIE matrix
            array{i,j,6}=squeeze(data(8,:,:));    % place in ntraitsxntraits cell array
            
            vector=AGE{nyr}(:,ii);             % extract vector for PHY ii
            data=gen_fcns.v2f(vector,vi,vj,vk); % convert to 36x36 GEnIE matrix
            arrage{i,j}=squeeze(data(8,:,:));    % place in ntraitsxntraits cell array
        end
    end
    extant=cell2mat(arrextant);
    matrix=cell2mat(array);         % catenate cell-array to (ntraitsx36)x(ntraitsx36) matrix
    age=cell2mat(arrage);
    age(age<0)=0;
    
    for j=1:3
        rgblim(j)=max(abs([prctile(matrix(:,:,j),5) prctile(matrix(:,:,j),95)])) ;
        matrix(:,:,j)=(matrix(:,:,j) .* extant+rgblim(j))./(2*rgblim(j));
    end
    
    axl=25.*sqrt(nyr);
    
    r=matrix(:,:,1);
    g=matrix(:,:,2);
    b=matrix(:,:,3);
    x=matrix(:,:,4) .* extant ./ axl;
    y=matrix(:,:,5) .* extant ./ axl;
    z=matrix(:,:,6) .* extant ./ axl;
    
    clf
    scatter3(x(:),y(:),z(:),25,[r(:) g(:) b(:)],'filled')
%     scatter3(x(~isnan(x)),y(~isnan(x)),z(~isnan(x)),10,cmap(round(age(~isnan(x))+1),:),'filled')
    hold on
    axis([-1 1 -1 1 -1 1])
    
    set(gca,'XTick',[])
    set(gca,'YTick',[])
    set(gca,'ZTick',[])
    box on
    
    title(num2str(nyr))
    
    drawnow
    
    fname=[output_dir '/Figures/XYZRGB_' num2str(nyr,'%04i') '.png'];
    export_fig(fname)
end

