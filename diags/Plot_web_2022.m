clear

%% Options

% path to model output files
outdir = '~/GitHub/omg/output/OMG_output_2022-11-22_Run_0001';

% minimum plotted biomass
threshold = 1e-6;

% flag to plot plankton C biomass
plot_phy=1;
% flag to plot plankton RGB colour trait
plot_rgb=0;

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
lat    = parameters.ocn_pars.lat_edges;
lon    = parameters.ocn_pars.lon_edges;
nlat   = size(parameters.ocn_pars.lat,1);
nlon   = size(parameters.ocn_pars.lon,1);
npopns = nsize*nlon*ntroph*nlat;
nrgb   = parameters.eco_pars.nrgb;

% get meta data
ESD                 = parameters.eco_pars.unq_ESD;
trophic             = unique(parameters.eco_pars.trophic);
n_sub_tsteps        = parameters.gen_pars.n_sub_tsteps;
save_intra_annual_n = parameters.gen_pars.save_intra_annual_n;

% define tick locations
xticks = (0:nlon)+0.5;
yticks = (0:nlat)+0.5;

%% open video objects for movies

vidObj3           = VideoWriter([outdir '/FoodWeb.mp4'], 'MPEG-4');
vidObj3.FrameRate = fps;
vidObj3.Quality   = 100;
open(vidObj3);

%% plot frames year by year

for t=1:ntime

    % rearrange dimensions - order: lat, size, lon, troph
    PHY_g = squeeze(sum(sum(PHY(t,:,:,:,:),2),3));

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

    f3 = figure(3);                                                             % change active figure
    f3 = clf(f3);                                                               % clear figure
    f3.Position=[figpos];                                            % Position figure
    ax1 = axes('Parent',f3);                                                    % setup axes
    imagesc(log10(PHY_g));                                                      % plot biomass
    axis xy                                                                     % flip vertical axis
    set(f3,'color','w');                                                        % set figure background to white
    ax1.Position=[0.05 0.05 0.875 0.915];                                       % expand axes within window
    ax1.XTick = xticks;                                                         % set x ticks
    ax1.YTick = yticks;                                                         % set y ticks
    grid on
    ax1.XTickLabel=num2str(lon,'%6.0f');                                        % set x tick labels
    ax1.YTickLabel=num2str(lat,'%6.0f');                                        % set y tick labels
    ax1.TickLength=[0 0];                                                       % make ticks invisble
    ax1.FontSize=14;                                                            % set font size
    hold on                                                                     % hold figure for plotting lines
    xlabel('Longitude','FontSize',16);                                          % label longitude
    ylabel('Latitude','FontSize',16);                                           % label latitude
    title('Total Phytoplankton biomass')                                        % add title
    text(1.03,max(yticks)*1.015,timestr,'FontSize',16)                          % Add simulation year
    caxis([-4 log10(ceil(nanmax(PHY_g(:))))])                     % set colour scale
    ch1=colorbar('EastOutside','Position',[0.935 0.5075 0.012 0.4575]);         % add colorbar
    ch1.TickLabels=cellstr(""+(10.^(ch1.Ticks)));                               % Change tick labels from log10
    ch1.TickLabels{1} = ['≤' ch1.TickLabels{1}];                                % Add ≤ for minimum on colour scale
    ch1.TickLabels{end} = [ch1.TickLabels{end} ' mmol P m^{-3}'];               % Add biomass units
 
    % rearrange dimensions - order: lat, size, lon, troph
    PHY_p = permute(PHY(t,:,:,:,:),[2 4 3 5 1]);
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


    xsub = linspace(0.1,0.9,ntroph)';
    ysub = linspace(0.1,0.9,nsize)';
    x = reshape(xsub + (0:nlon-1),[],1)+0.5;
    y = reshape(ysub + (0:nlat-1),[],1)+0.5;
    [xx yy] = meshgrid(x,y);
    i_extant=find(PHY_n);


    % rearrange dimensions - order: lat, size, lon, troph
    RGB_p = permute(RGB(t,:,:,:,:,:),[2 4 3 5 6 1]);
    % combine size+lon (rows) and troph+lat (columns)
    RGB_a = reshape(RGB_p,nsize*nlon,ntroph*nlat,nrgb);
    % normalise to range across all colour channels
    RGB_n = (RGB_a - nanmin(RGB_a(:))) ./ (nanmax(RGB_a(:)) - nanmin(RGB_a(:)));

    % pad i_extinct index for ngb colour dimensions
    rgb_extant = reshape(i_extant + (0:nrgb-1).*npopns,[],1); % (uses implicit array expansion)
    
    rgbclr = reshape(RGB_n,[],nrgb);

    scatter(xx(i_extant),yy(i_extant),PHY_n(i_extant).^(2/3).*50,rgbclr(i_extant,:),'o','filled')
    



    drawnow                                                                     % draw figure
    currFrame = getframe(f3);                                                   % get frame for movie
    writeVideo(vidObj3,currFrame);                                              % write to video

end

%% close video objects for movies

close(vidObj3);

