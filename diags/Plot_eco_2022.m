clear

%% Options

% path to model output files
outdir = '~/GitHub/omg/output/OMG_output_2023-05-10_Run_0023';

% minimum plotted biomass
threshold = 1e-6;

% flag to plot plankton C biomass
plot_phy=0;
% flag to plot plankton RGB colour trait
plot_rgb=1;

figpos = [1 58 1440 739];

% colours of grid lines
lncolr_phy = ones(1,3).*0.75;
lncolr_rgb = ones(1,3).*0.25;

% movie frame rate
fps = 10;

%% load data

load([outdir '/Parameters.mat']);                           % model parameters
PHY    = ncread([outdir '/omg_fields_netcdf.nc'],'PHY');    % plankton P biomass (moles per kg)
RGB = ncread([outdir '/omg_fields_netcdf.nc'],'RGB');       % palnkton RGB colour trait
m_time = ncread([outdir '/omg_fields_netcdf.nc'],'time');   % model time (years)

PHY = PHY.*1e6; % convert from mol/kg to mmol/m^3

% find number of completed years
i_time = find(PHY>0,1,'last');                  % get linear index of last non-zero datum
[ntime,~,~,~,~,~] = ind2sub(size(PHY),i_time);  % get time index from linear index

% get dimensions of netcdf data
nsize  = parameters.eco_pars.nsize;
ntroph = parameters.eco_pars.ntroph;
nlat   = size(parameters.ocn_pars.lat,1);
nlon   = size(parameters.ocn_pars.lon,1);
npopns = nsize*nlon*ntroph*nlat;
nrgb   = parameters.eco_pars.nrgb;

% get meta data
ESD                 = parameters.eco_pars.unq_ESD;
trophic             = unique(parameters.eco_pars.trophic);
n_sub_tsteps        = parameters.gen_pars.n_sub_tsteps;
save_intra_annual_n = parameters.gen_pars.save_intra_annual_n;

% define points for manually drawn axis grid and tick marks
xticks = linspace(0,ntroph*nlon,ntroph+1)+0.5;
yticks = linspace(0,nsize *nlat,nsize +1)+0.5;
xticksmiddle = xticks(2:end)-diff(xticks)./2;
yticksmiddle = yticks(2:end)-diff(yticks)./2;



%% open video objects for movies

if plot_phy
    vidObj1           = VideoWriter([outdir '/Biomass.mp4'], 'MPEG-4');
    vidObj1.FrameRate = fps;
    vidObj1.Quality   = 100;
    open(vidObj1);
end

if plot_rgb
    vidObj2           = VideoWriter([outdir '/RGB.mp4'], 'MPEG-4');
    vidObj2.FrameRate = fps;
    vidObj2.Quality   = 100;
    open(vidObj2);
end
%% plot frames year by year

for t=1:ntime

    % rearrange dimensions - order: lat, size, lon, troph
    PHY_p = permute(PHY(t,:,:,:,:),[4 2 5 3 1]);
    % combine size+lon (rows) and troph+lat (columns)
    PHY_a = reshape(PHY_p,nsize*nlon,ntroph*nlat);
    % find data below threshold
    [i_extinct] = find(PHY_a<threshold);
    % set data below threshold to threshold
    PHY_a(i_extinct)=threshold;
    % take log10 of biomass
    PHY_l=log10(PHY_a);
    % normalise between 'threshold' and maximum value
    PHY_n = (PHY_l - log10(threshold)) / (nanmax(PHY_l(:)) - log10(threshold));

    % get time in years
    t_time = m_time(t);
    year = floor(t./save_intra_annual_n);
    % create title string
    if save_intra_annual_n==1
        timestr = [num2str(year) ' years'];
    else
        day  = 360 * rem(t_time,1);
        timestr = [num2str(year) ' years and ' num2str(day) ' days'];
    end

    if plot_phy

        f1 = figure(1);                                                             % change active figure
        f1 = clf(f1);                                                               % clear figure
        f1.Position=figpos;                                            % Position figure
        ax1 = axes('Parent',f1);                                                    % setup axes
        imagesc(log10(PHY_a));                                                      % plot biomass        
        axis xy                                                                     % flip vertical axis     
        set(f1,'color','w');                                                        % set figure background to white
        ax1.Position=[0.05 0.05 0.875 0.915];                                       % expand axes within window
        ax1.XTick = xticksmiddle;                                                   % set x ticks
        ax1.YTick = yticksmiddle;                                                   % set y ticks
        ax1.XTickLabel=num2str(trophic,'%6.2f');                                    % set x tick labels
        ax1.YTickLabel=num2str(str2num(num2str(ESD,2)),'%6.2f');                    % set y tick labels
        ax1.TickLength=[0 0];                                                       % make ticks invisble
        ax1.FontSize=14;                                                            % set font size
        hold on                                                                     % hold figure for plotting lines
        plot([xticks(1) xticks(end)],[yticks;yticks]','color',lncolr_phy,'LineW',1);% plot horizontal grid lines
        plot([xticks;xticks]',[yticks(1) yticks(end)],'color',lncolr_phy,'LineW',1);% plot vertical grid lines
        xlabel('Trophic strategy','FontSize',16);                                   % label trophic indices
        ylabel(['Size (' char(181) 'm)'],'FontSize',16);                            % label sizes
        title('Phytoplankton biomass')                                              % add title
        text(1.03,max(yticks)*1.015,timestr,'FontSize',16)                          % Add simulation year
        caxis([log10(threshold) log10(ceil(nanmax(PHY_a(:))))])                     % set colour scale
        ch1=colorbar('EastOutside','Position',[0.935 0.5075 0.012 0.4575]);         % add colorbar
        ch1.TickLabels=cellstr(""+(10.^(ch1.Ticks)));                               % Change tick labels from log10
        ch1.TickLabels{1} = ['≤' ch1.TickLabels{1}];                                % Add ≤ for minimum on colour scale
        ch1.TickLabels{end} = [ch1.TickLabels{end} ' mmol P m^{-3}'];               % Add biomass units
        hold off
        
        drawnow                                                                     % draw figure
        currFrame = getframe(f1);                                                   % get frame for movie
        writeVideo(vidObj1,currFrame);                                              % write to video
    end

    if plot_rgb

        % rearrange dimensions - order: lat, size, lon, troph
        RGB_p = permute(RGB(t,:,:,:,:,:),[4 2 5 3 6 1]);
        % combine size+lon (rows) and troph+lat (columns)
        RGB_a = reshape(RGB_p,nsize*nlon,ntroph*nlat,nrgb);
        % pad i_extinct index for ngb colour dimensions
        rgb_extinct = reshape(i_extinct + (0:nrgb-1).*npopns,[],1); % (uses implicit array expansion)
        % set data with PHY below threshold to 0
        RGB_a(rgb_extinct) = 0; 
        % normalise to range across all colour channels
        RGB_n = (RGB_a - nanmin(RGB_a(:))) ./ (nanmax(RGB_a(:)) - nanmin(RGB_a(:)));
        % weight colours by normalised biomass
        RGB_w = RGB_n.*PHY_n;  
    
        f2 = figure(2);                                                             % change active figure
        f2 = clf(f2);                                                               % clear figure
        f2.Position=figpos;                                            % Position figure
        ax2a = axes('Parent',f2);                                                   % setup axes
        image(RGB_w);                                                               % plot weighted colours
        axis xy                                                                     % flip vertical axis    
        set(f2,'color','w');                                                        % set figure background to white  
        ax2a.Position=[0.05 0.05 0.875 0.915];                                      % expand axes within window
        ax2a.XTick = xticksmiddle;                                                  % set x ticks
        ax2a.YTick = yticksmiddle;                                                  % set y ticks
        ax2a.XTickLabel=num2str(trophic,'%6.2f');                                   % set x tick labels
        ax2a.YTickLabel=num2str(str2num(num2str(ESD,2)),'%6.2f');                   % set y tick labels
        ax2a.TickLength=[0 0];                                                      % make ticks invisble
        ax2a.FontSize=14;                                                           % set font size
        hold on                                                                     % hold figure for plotting lines
        plot([xticks;xticks]',[yticks(1) yticks(end)],'color',lncolr_rgb,'LineW',1);% plot horizontal grid lines
        plot([xticks(1) xticks(end)],[yticks;yticks]','color',lncolr_rgb,'LineW',1);% plot vertical grid lines
        ylabel(['Size (' char(181) 'm)'],'FontSize',16);                            % label trophic indices
        xlabel('Trophic strategy','FontSize',16);                                   % label sizes
        title('Neutral colour trait')                                               % add title
        text(1.03,max(yticks)*1.015,timestr,'FontSize',16)                          % Add simulation year
        caxis([log10(threshold) log10(ceil(nanmax(PHY_a(:))))])                     % set colour scale
        ch2=colorbar('EastOutside','Position',[0.935 0.5075 0.012 0.4575]);         % add colorbar
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
end

%% close video objects for movies
if plot_phy
    close(vidObj1);
end

if plot_rgb
    close(vidObj2);
end


