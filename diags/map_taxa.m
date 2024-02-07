function map_taxa(distances,Cutoff,biomass,parameters)

    
   
    
    
    

    nclust = numel(unique(c));
    % fit clusters into 'biomass' sized matrix
    Clus = biomass.*0;
    Clus(kk) = c;
    % rearrange dimensions - order: lat, T2, lon, T1
    Clus = permute(Clus(:,:,:,:),[2 4 1 3]);
    % combine size+lon (rows) and troph+lat (columns)
    Clus = reshape(Clus,nlat*nT2,nlon*nT1);  

    % plot distance matrix
    t_square = squareform(t_est);
    ax3=subplot(223);
    imagesc(t_square(outperm,outperm))
    ax4.Colormap=turbo(64);
    axis square
    colorbar


    cmap=rand(nclust+1,3);
    cmap(1,:)=[0 0 0];
    ClusterRGB = zeros([size(Clus),3]).*NaN;
    for i=1:size(Clus,1)
        for j=1:size(Clus,2)
            if ~isnan(Clus(i,j))
                ClusterRGB(i,j,:)=cmap(Clus(i,j)+1,:);
            end
        end
    end
    ClusterRGB=ClusterRGB.*PHY_n;


    % map clusters
    f4 = figure(4);
    ax4a = axes('Parent',f4);                                                   % setup axes
    land_mask=~isnan(Clus);              % blank array for land mask
    imagesc(ClusterRGB,'AlphaData',land_mask); % plot OTUs
    set(gca,'color',0.75.*[1 1 1]);            % set axis background colour to grey for land
    axis xy
    ax4a.XTick = xticksmiddle;                                                  % set x ticks
    ax4a.YTick = yticksmiddle;
    hold on
    plot([xticks(1) xticks(end)],[yticks;yticks],'color',lncolr_phy,'LineW',1); % plot horizontal grid lines
    plot([xticks;xticks],[yticks(1) yticks(end)],'color',lncolr_phy,'LineW',1); % plot vertical grid lines
    ax4a.XTickLabel=num2str(T1,'%6.2f');                                        % set x tick labels
    ax4a.YTickLabel=num2str(T2,'%6.2f');                                        % set y tick labels
    cmap=rand(nclust,3);
%     cmap=turbo(nclust);
    cmap(1,:)=[0 0 0];
    ax4a.Colormap=cmap;
    caxis([0 nclust]+0.5)
    colorbar
    xlabel(['Year = ' num2str(t)])
    title('Map of OTUs')


end