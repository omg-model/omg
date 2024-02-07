clear

%% Options

% path to model output files
outdir = '~/GitHub/omg/output/OMG_output_2023-07-19_Run_0001';

% minimum plotted biomass
threshold = 1e-6;

%% JGOFS locations
JGOFS{1}.latlon = [+023,-158]; JGOFS{1}.name = 'HOT';
JGOFS{2}.latlon = [+032,-064]; JGOFS{2}.name = 'BATS';
JGOFS{3}.latlon = [+000,-140]; JGOFS{3}.name = 'EQPAC';
JGOFS{4}.latlon = [+016,+062]; JGOFS{4}.name = 'ARABIAN';
JGOFS{5}.latlon = [+047,-019]; JGOFS{5}.name = 'NABE';
JGOFS{6}.latlon = [+050,-145]; JGOFS{6}.name = 'STNP';
JGOFS{7}.latlon = [-051,+068]; JGOFS{7}.name = 'KERFIX';
JGOFS{8}.latlon = [-062,-170]; JGOFS{8}.name = 'APFZ';
JGOFS{9}.latlon = [-070,-180]; JGOFS{9}.name = 'ROSS';

%% load data

load([outdir '/Parameters.mat']);                         % model parameters
PHY = ncread([outdir '/omg_fields_netcdf.nc'],'PHY');     % plankton P biomass (moles per kg)
RGB = ncread([outdir '/omg_fields_netcdf.nc'],'RGB');     % palnkton RGB colour trait
m_time = ncread([outdir '/omg_fields_netcdf.nc'],'time'); % model time (timesteps)

PHY = PHY.*1e6; % convert from mol/kg to mmol/m^3

% find number of completed timesteps
i_time = find(PHY>0,1,'last');                % get linear index of last non-zero datum
[~,~,~,~,ntime] = ind2sub(size(PHY),i_time);  % get time index from linear index

%% get dimensions of netcdf data

% trait dimensions
nT1 = parameters.eco_pars.nT1;
nT2 = parameters.eco_pars.nT2;

% spatial dimensions
nlat   = size(parameters.ocn_pars.lat,1);
nlon   = size(parameters.ocn_pars.lon,1);

% total number of local populations in metacommunity
npopns = nT1*nT2*nlat*nlon;

% number of rgb traits
nrgb   = parameters.eco_pars.nrgb;

Tname1=parameters.eco_pars.Tname1;
Tname2=parameters.eco_pars.Tname2;



% get meta data
T1                  = parameters.eco_pars.T1;
T2                  = parameters.eco_pars.T2;
n_sub_tsteps        = parameters.gen_pars.n_sub_tsteps;
save_intra_annual_n = parameters.gen_pars.save_intra_annual_n;

% define points for manually drawn axis grid and tick marks
xticks = linspace(0,nT1*nlon,nT1+1)+0.5;
yticks = linspace(0,nT2 *nlat,nT2+1)+0.5;
xticksmiddle = xticks(2:end)-diff(xticks)./2;
yticksmiddle = yticks(2:end)-diff(yticks)./2;

% find nearest grid cells to JGOFS sites
for i=1:size(JGOFS,2)
    [~,ilat] = min(abs(JGOFS{i}.latlon(1)-parameters.ocn_pars.lat));
    [~,ilon] = min(abs(JGOFS{i}.latlon(2)-parameters.ocn_pars.lon));
    JGOFS{i}.grid_ilat = ilat;
    JGOFS{i}.grid_ilon = ilon;
    JGOFS{i}.grid_lat  = parameters.ocn_pars.lat(ilat);
    JGOFS{i}.grid_lon  = parameters.ocn_pars.lon(ilon);
end

% load timeseries data
load('/Users/baw1d17/Data/In_situ/Time_series/KLEYPAS/Kleypas.mat')


%% plot data from last year

t_in_year=(ntime-save_intra_annual_n):ntime;
    
%


% rearrange dimensions - order: lat, T2, lon, T1
PHY_p = permute(PHY(:,:,:,:,t_in_year),[2 4 1 3 5]);
% sum across phenotypes
PHY_sum = squeeze(sum(sum(PHY_p,2),4));
% combine size+lon (rows) and troph+lat (columns)
PHY_a = reshape(PHY_p,nlat*nT2,nlon*nT1,save_intra_annual_n+1);
% take log10 of biomass
PHY_l=log10(PHY_a);


%%
chl2P = 1.59*16;

% clf
for j=1:9
    subplot(3,3,j)
    ilat = JGOFS{j}.grid_ilat;
    ilon = JGOFS{j}.grid_ilon;
    site_data = squeeze(PHY_p(ilat,:,ilon,:,:));

    plot(0:save_intra_annual_n,squeeze(sum(site_data(:,1,:),1)).*chl2P,'g-','LineWidth',2)
    hold on
    plot(0:save_intra_annual_n,squeeze(sum(site_data(:,2,:),1)),'m-','LineWidth',2)
    title(JGOFS{j}.name)
    xlim([0 save_intra_annual_n])
    ylim([1e-3 1e1])
    set(gca,'YScale','log')
    set(gca,'XTick',linspace(0,save_intra_annual_n,13))
    set(gca,'XTickLabel',{'J','F','M','A','M','J','J','A','S','O','N','D',''})

    obsdata = eval([JGOFS{j}.name '_data']);
    plot(obsdata.chlday./30,obsdata.chldata,'g.')
end


return

%%




while true
    for t=1:save_intra_annual_n

        axesm ('MapProjection','putnins5',...
       'Frame', 'off',...
       'Grid', 'off',...
       'MapLonLimit',parameters.ocn_pars.lon_edges([1 end]));
        pcolorm(parameters.ocn_pars.lat_edges,parameters.ocn_pars.lon_edges,log10(PHY_sum(:,:,t)))
        title(t)
        colorbar
        caxis([-3 1])

        hold on
        for j=1:9
            plotm(JGOFS{j}.grid_latlon(1),JGOFS{j}.grid_latlon(2),'m.')
        end
        
        drawnow
    end
end
%%


    % normalise between 'threshold' and maximum value
    PHY_n = (PHY_l - log10(threshold)) / (nanmax(PHY_l(:)) - log10(threshold));

    PHY_intT1 = sum(PHY_p,4);
    PHY_intT1 = reshape(PHY_intT1,nlat*nT2,nlon);
    PHY_intT2 = sum(PHY_p,2);
    PHY_intT2 = reshape(PHY_intT2,nlat,nlon*nT1);
    PHY_int   = squeeze(sum(PHY_p,[2 4]));

    % find data below threshold
    [i_extinct] = find(PHY_a<threshold);
    % set data below threshold to threshold
    PHY_a(i_extinct)=threshold;

    % get time in years
    t_time = m_time(t);
    year = floor(t./save_intra_annual_n);
    % create title string
    if year==1
        timestr = [num2str(year) ' year'];
    else
        timestr = [num2str(year) ' years'];
    end
    if save_intra_annual_n~=1
        day  = 360 * rem(t_time,1);
        if day==1
            timestr = [timestr ' and ' num2str(day) ' day'];
        else
            timestr = [timestr ' and ' num2str(day) ' days'];
        end
    end

    if plot_phy

        f1 = figure(1);                                                             % change active figure
        f1 = clf(f1);                                                               % clear figure
        f1.Position=fig1pos;                                                        % Position figure
        colormap(turbo);                                                            % Turbo colormap
        set(f1,'color','w');                                                        % set figure background to white

        % Plot individual populations as T1xT2 matrix
        ax1a = axes('Parent',f1);                                                   % setup axes
        land_mask=ones(size(PHY_a));                                                % blank array for land mask
        land_mask(isnan(PHY_a))=0;                                                  % set land values to zero (ocean=1) 
        imagesc(log10(PHY_a),'AlphaData',land_mask);                                % plot biomass  
        set(gca,'color',0.75*[1 1 1]);                                              % set axis background colour to grey for land  
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
        caxis([log10(threshold) log10(ceil(nanmax(PHY_a(:))))])                     % set colour scale

        title('Phytoplankton biomass','FontSize',16)                                % add title
        text(1.03,max(yticks)*1.015,timestr,'FontSize',16)                          % Add simulation year
        ch1=colorbar('EastOutside','Position',...
                     [ax1a.Position(1)+ax1a.Position(3)+0.01 ...
                      ax1a.Position(2)+ax1a.Position(4)/2 ...
                      0.012 ...
                      ax1a.Position(4)/2],...
                      'FontSize',16);                                               % add colorbar
        ch1.TickLabels=cellstr(""+(10.^(ch1.Ticks)));                               % Change tick labels from log10
        ch1.TickLabels{1} = ['≤' ch1.TickLabels{1}];                                % Add ≤ for minimum on colour scale
        ch1.TickLabels{end} = [ch1.TickLabels{end} ' mmol P m^{-3}'];               % Add biomass units
        hold off

        % plot biomasses integrated across T1
        ax1b = axes('Parent',f1);                                                   % add new axis for integrated T1
        land_mask=ones(size(PHY_intT1));                                            % blank array for land mask
        land_mask(isnan(PHY_intT1))=0;                                              % set land values to zero (ocean=1) 
        imagesc(log10(PHY_intT1),'AlphaData',land_mask);                            % plot biomass  
        set(gca,'color',0.75*[1 1 1]);                                              % set axis background colour to grey for land
        axis xy                                                                     % flip vertical axis     
        ax1b.Position = ax1a.Position;                                              % reposition new axis over main
        panel_width   = ax1a.Position(4)./nT1;                                       % get width of panels in main axis
        ax1b.Position([1 3]) = [ax1a.Position(1)-panel_width-0.005 panel_width];    % set to left of main axis
        ax1b.YTick = yticksmiddle;                                                  % set y ticks
        ax1b.TickLength=[0 0];                                                      % make ticks invisble
        ax1b.XTickLabel=num2str('');                                                % no x tick labels
        ax1b.YTickLabel=num2str(T2,'%6.2f');                                        % set y tick labels
        ax1b.FontSize=14;                                                           % set font size
        hold on                                                                     % hold figure for plotting lines
        plot([xticks(1) xticks(end)],[yticks;yticks],'color',lncolr_phy,'LineW',1); % plot horizontal grid lines
        ylabel(Tname2,'FontSize',16);                                               % label sizes
        caxis([log10(threshold) log10(ceil(nanmax(PHY_a(:))))]);                    % set colour scale

        % plot biomasses integrated across T2
        ax1c = axes('Parent',f1);                                                   % add new axis for integrated T1
        land_mask=ones(size(PHY_intT2));                                            % blank array for land mask
        land_mask(isnan(PHY_intT2))=0;                                              % set land values to zero (ocean=1) 
        imagesc(log10(PHY_intT2),'AlphaData',land_mask);                            % plot biomass  
        set(gca,'color',0.75*[1 1 1]);                                              % set axis background colour to grey for land
        axis xy                                                                     % flip vertical axis     
        ax1c.Position = ax1a.Position;                                              % reposition new axis over main
        panel_height = ax1a.Position(3)./nT2;                                       % get height of panels in main axis
        ax1c.Position([2 4]) = [ax1a.Position(2)-panel_height-0.01 panel_height];   % set underneath main axis
        ax1c.XTick = xticksmiddle;                                                  % set x ticks
        ax1c.YTick=[];                                                              % no y ticks  
        ax1c.TickLength=[0 0];                                                      % make ticks invisble
        ax1c.XTickLabel=num2str(T1,'%6.2f');                                        % set x tick labels
        ax1c.FontSize=14;                                                           % set font size
        hold on                                                                     % hold figure for plotting line
        plot([xticks;xticks],[yticks(1) yticks(end)],'color',lncolr_phy,'LineW',1); % plot vertical grid lines
        xlabel(Tname1,'FontSize',16);                                               % label trophic indices
        caxis([log10(threshold) log10(ceil(nanmax(PHY_a(:))))])                     % set colour scale

        % plot total integrated biomss
        ax1d = axes('Parent',f1);                                                   % add new axis for integrated T1
        land_mask=ones(size(PHY_int));                                              % blank array for land mask
        land_mask(isnan(PHY_int))=0;                                                % set land values to zero (ocean=1) 
        imagesc(log10(PHY_int),'AlphaData',land_mask);                              % plot biomass  
        set(gca,'color',0.75*[1 1 1]);                                              % set axis background colour to grey for land
        axis xy                                                                     % flip vertical axis     
        ax1d.Position = ax1a.Position;                                              % reposition new axis over main
        ax1d.Position = [ax1a.Position(1)-panel_width-0.005 ...
                         ax1a.Position(2)-panel_height-0.01 ...
                         panel_width ...
                         panel_height];   % set underneath main axis
        ax1d.XTick = [];                                                            % set x ticks
        ax1d.YTick = [];                                                            % no y ticks 
        xlabel('Integrated biomass','FontSize',16);                                               % label trophic indices
        caxis([log10(threshold) log10(ceil(nanmax(PHY_a(:))))])                     % set colour scale
      
                                                                

        
        drawnow                                                                     % draw figure
        currFrame = getframe(f1);                                                   % get frame for movie
        writeVideo(vidObj1,currFrame);                                              % write to video
    end

    if plot_rgb

        % rearrange dimensions - order: lat, T2, lon, T1
        RGB_p = permute(RGB(:,:,:,:,:,t),[2 4 1 3 5]);
        % combine size+lon (rows) and troph+lat (columns)
        RGB_a = reshape(RGB_p,nlat*nT2,nlon*nT1,nrgb);
        % pad i_extinct index for ngb colour dimensions
        rgb_extinct = reshape(i_extinct + (0:nrgb-1).*npopns,[],1); % (uses implicit array expansion)
        % set data with PHY below threshold to NaN
        RGB_a(rgb_extinct) = NaN; 
        % remove outliers
        RGB_a = filloutliers(RGB_a,'clip','median'); % Replaces all elements more than 3 scaled MAD from the median
                                                     % with value 3 scaled MAD from the median
        % normalise to range across all colour channels
        RGB_n = (RGB_a - min(RGB_a(:))) ./ (max(RGB_a(:)) - min(RGB_a(:)));
        % weight colours by normalised biomass
        RGB_w = RGB_n.*PHY_n;  
    
        f2 = figure(2);                                                             % change active figure
        f2 = clf(f2);                                                               % clear figure
        f2.Position=fig2pos;                                                        % Position figure
        ax2a = axes('Parent',f2);                                                   % setup axes
        land_mask=ones(size(PHY_a));                                                % blank array for land mask
        land_mask(isnan(PHY_a))=0;                                                  % set land values to zero (ocean=1) 
        imagesc(RGB_w,'AlphaData',land_mask);                                       % plot weighted colours
        set(gca,'color',0.75*[1 1 1]);                                              % set axis background colour to grey for land
        axis xy                                                                     % flip vertical axis    
        set(f2,'color','w');                                                        % set figure background to white  
        ax2a.Position=ax1a.Position;                                                % put axis in same place as Fig 1 main axis
        ax2a.XTick = xticksmiddle;                                                  % set x ticks
        ax2a.YTick = yticksmiddle;                                                  % set y ticks
        ax2a.XTickLabel=num2str(T1,'%6.2f');                                        % set x tick labels
        ax2a.YTickLabel=num2str(T2,'%6.2f');                                        % set y tick labels
        ax2a.TickLength=[0 0];                                                      % make ticks invisble
        ax2a.FontSize=14;                                                           % set font size
        hold on                                                                     % hold figure for plotting lines
        plot([xticks;xticks],[yticks(1) yticks(end)],'color',lncolr_rgb,'LineW',1); % plot horizontal grid lines
        plot([xticks(1) xticks(end)],[yticks;yticks],'color',lncolr_rgb,'LineW',1); % plot vertical grid lines
        ylabel(Tname2,'FontSize',16);                            % label trophic indices
        xlabel(Tname1,'FontSize',16);                                               % label sizes
        title('Neutral colour trait')                                               % add title
        text(1.03,max(yticks)*1.015,timestr,'FontSize',16)                          % Add simulation year
        caxis([log10(threshold) log10(ceil(nanmax(PHY_a(:))))])                     % set colour scale
        ch2=colorbar('EastOutside','Position',ch1.Position,'FontSize',16);          % add colorbar
        ch2.TickLabels=cellstr(""+(10.^(ch2.Ticks)));                               % Change tick labels from log10
        ch2.TickLabels{1} = ['≤' ch2.TickLabels{1}];                                % Add ≤ for minimum on colour scale
        ch2.TickLabels{end} = [ch2.TickLabels{end} ' mmol P m^{-3}'];               % Add biomass units

        ax2b = axes('Parent',f2);                                                   % add new axis for multiple colour scale
        nc = 256;                                                                   % number of discrete colours/shades in colour bar
        cmap=hsv2rgb([linspace(0,1,nc)' ones(nc,1).*0.5 ones(nc,1)]);               % generate rgb colourmap (all hues, 50% saturation)
        cmap=repmat(permute(cmap,[3 1 2]),[nc 1 1]).*linspace(1,0,nc)';             % add another dimension for variable values
        ih2=image(ax2a.XLim,clim(ax2a),cmap);                                       % plot multiple colour scale as an image
        ylim(clim(ax2a))                                                            % match y axis to colour scale
        ax2b.Position = ch2.Position;                                               % reposition new image over old colour scale
        ax2b.YTick = ch2.Ticks;                                                     % set position of y tick marks
        ax2b.YTickLabel = '';                                                       % remove y tick labels
        ax2b.XTick = [];                                                            % remove x tick marks
        
        drawnow                                                                     % draw figure
        currFrame = getframe(f2);                                                   % get frame for movie
        writeVideo(vidObj2,currFrame);                                              % write to video
    end
% end

%% close video objects for movies
if plot_phy
    close(vidObj1);
end

if plot_rgb
    close(vidObj2);
end


