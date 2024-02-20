%% Ocean Circulation & Climate Parameters

%% load transport matrix
Ainput = load([gen_pars.TM_path '/A.mat']); % ** needs developing with all matrices **
% and associated matrix variables
load([gen_pars.TM_path '/matrix_vars.mat']); % Matrix index parameters

% FIX FOR INCONSISTENT FORMAT OF SEASONAL TRANSPORT MATRICES
if ~iscell(Ainput) % if not already a cell - relax!
    if isstruct(Ainput) % if a structural array
        nms = char(fieldnames(Ainput));           % get TM matrix names
        [~,idx] = sort(str2num(nms(:,2:end)));    % remove prefix letter and sort by numbers
        if isempty(idx); 
            idx = 1;
            if iscell(Ainput.A)
                Ainput.A=cell2mat(Ainput.A);
            end
        end
        A = struct2cell(Ainput);                  % convert to cell array
        A = A(idx);                               % out in correct order  
    else % assuming it is 1 double array
        A = {Ainput}; % convert to cell
    end
end

% JDW - clear Ainput?

% and associated matrix variables - JDW mini bodge
load([gen_pars.TM_path '/matrix_vars.mat']); % Matrix index parameters

%% 3D climate and circulation variables from GEnIE timeslices:

load([gen_pars.TM_path '/forcing_vars.mat'])
OCEAN.T=forcing_vars.T;
OCEAN.S=forcing_vars.S;
OCEAN.solfor=forcing_vars.solfor;
OCEAN.seaice=forcing_vars.seaice;
OCEAN.MLD=forcing_vars.MLD;
OCEAN.wspeed=forcing_vars.wspeed;
OCEAN.atm_temp=forcing_vars.atm_temp;

% calculate density (rho - using equations in gem_util.90 and gem_cmn.f90 in genie)
OCEAN.rho=1000.0 + (0.7968 * OCEAN.S - 0.0559 * OCEAN.T - 0.0063 * OCEAN.T.^2 + 3.7315E-05 * OCEAN.T.^3);

% set units and precalculate
OCEAN.PAR0              =  OCEAN.solfor ...
                           .* (1-OCEAN.seaice./100).* bgc_pars.parfrac;    % photosynthetically available radiation (PAR)
OCEAN.solfor            = OCEAN.solfor./1368;                              % solar forcing (normalised to ~solar constant)
OCEAN.seaice            = 1-(OCEAN.seaice./100);                           % sea ice
OCEAN.atm_temp          = OCEAN.atm_temp+273.15;                           % atmosphere temperature from deg C to deg K


%% Transport Matrix Metadata

% n time-steps
ocn_pars.n_A    =  numel(A);
% Matrix Indexing
ocn_pars.i      =  v_index.i;
ocn_pars.j      =  v_index.j;
ocn_pars.k      =  v_index.k;
ocn_pars.rk     =  v_index.rk; % inverted depth index
ocn_pars.Ib     =  Ib;
ocn_pars.Ii     =  Ii;
ocn_pars.nb     =  nb;
ocn_pars.Idt    =  1;      
ocn_pars.ni     =  size(ocn_pars.Ib,1);

%% 1D grid metadata

load([gen_pars.TM_path '/grid_vars.mat'])

ocn_pars.zt_edges=grid_vars.zt_edges;
ocn_pars.depth=grid_vars.depth;
ocn_pars.lon=grid_vars.lon;
ocn_pars.lat=grid_vars.lat;
ocn_pars.lon_edges=grid_vars.lon_edges;
ocn_pars.lat_edges=grid_vars.lat_edges;
ocn_pars.A=grid_vars.A;
genie_time=grid_vars.genietime;

ocn_pars.A = double(ocn_pars.A);


%i_v = sub2ind(size(ocn_pars.A),ocn_pars.rk,ocn_pars.j,ocn_pars.i);

%% Generate remaining variables

% OMGWTF! hard-coded parameters that change with grid resolution ...
% OMGWTF! ocn_pars.A     = 393444720640.0; % (m2) equal area grid
% OMGWTF! ocn_pars.atm_A = 393444720640.0; 

% Some TMs do not have area, volume, mass data below the surface,
% so reverting to calculation based on zt_edges

% Find GEnIE grid area
%ocn_pars.A = unique(ocn_pars.A(i_v(ocn_pars.rk==1))); % use areas in surface only
ocn_pars.A = unique(ocn_pars.A(ocn_pars.rk==1)); % use areas in surface only
if numel(ocn_pars.A)>1; error('GEnIE grid is not equal area!!!! (It should be!)'); end
ocn_pars.atm_A = ocn_pars.A; % BAW Assuming atmospheric grid is same area as oceanic grid (as previosly assumed)

% Volume = Area * Depth
ocn_pars.V = diff(ocn_pars.zt_edges).*ocn_pars.A;

% broadcast 1D variables to global vector
ocn_pars.V = ocn_pars.V(ocn_pars.rk); 
ocn_pars.depth  = ocn_pars.depth(ocn_pars.rk); 

% calculate mass using volume loaded from fields_biogem_3d.nc
ocn_pars.M  = 1027.649.*ocn_pars.V;
ocn_pars.rM = 1./ocn_pars.M;

% get day of year for each forcing timestep
ocn_pars.timesteps = round([0+(1/ocn_pars.n_A)/2:(1/ocn_pars.n_A):1].*360)';
%ocn_pars.timesteps = repmat(ocn_pars.timesteps,gen_pars.runtime,1).*[1:1:gen_pars.runtime]';
%OCEAN.timesteps = rem(OCEAN.timesteps,1).*360; % get days of year
%OCEAN.timesteps=unique(round(OCEAN.timesteps));% round to day and find unique

% JDW - why load in time points from genie run?!?! OK so assumes netcdf
% timselices have seasonal saving points same as the TM...OK but iffy
% JDW - amended to match number of matrices (assuming they are evenly
% spaced which why not?!)

% get water columns
ocn_pars.wc=zeros(ocn_pars.nb,1);
count=1;
for ii=1:numel(ocn_pars.lon)
    for jj=1:numel(ocn_pars.lat)
        ind=find(ocn_pars.i==ii & ocn_pars.j==jj);
        if ~isempty(ind)
            ocn_pars.wc(ind,1)=count;
            count=count+1;
        end
    end
end

% get benthic grid cells
for n=1:max(ocn_pars.wc)
    ocn_pars.Iben(n,1)=find(ocn_pars.wc==n,1,'last');
end


%% Forcing data
ocn_pars.z0                = ocn_pars.zt_edges(2);                               % bottom of surface layer (m)
ocn_pars.atm_dh            = 7777.0;                                          % atmosphere thickness of single cell (m)
ocn_pars.atm_V             = ocn_pars.atm_A(:,1).*ocn_pars.atm_dh;                  % atmosphere volume (m3)
ocn_pars.atm_Vtot          = sum(ocn_pars.atm_V(:,1));                           % total atmosphere volume (m3)
ocn_pars.atm_mol           = 1.77e20;                                         % ???
ocn_pars.const_R           = 8.3145;                                          % ideal gas constant
ocn_pars.conv_Pa_atm       = 1/1.01325e5;                                     % see gem_cmn.f90
ocn_pars.conv_atm_mol      = repmat(ocn_pars.atm_V,1,ocn_pars.n_A) ...
                          ./(ocn_pars.conv_Pa_atm ...
                           .*ocn_pars.const_R.*OCEAN.atm_temp);               % convert atm to mol, see atchem.f90
ocn_pars.conv_mol_atm      = 1./ocn_pars.conv_atm_mol;                           % convert mol to atm, see atchem.f90
%OCEAN.gastransfer_a     = 0.31;                                            % ocn_pars.gastransfer_a; % BODGE !?!?!?!

%[OCEAN.Sc,OCEAN.Sol]=bgc_pars.calc_gasexchange_constants(OCEAN.T,OCEAN.S,OCEAN.rho,bgc_pars , I )

% JDW -> check gastransfer_a
% JDW -> check piston velocity calculation happens OK
                

%%


%% Transport matrix correction
if strcmp(bgc_pars.uptake_scheme,'eco')
    % Ecosystem model requires TM correction
    gen_pars.TM_correction = true;
end
if gen_pars.TM_correction 
    for n=1:numel(A)
        % The original transport matrix is given as a concentration flux per timestep. This
        % means that each element transforms an upstream concentration into a
        % downstream concentration, implicitly accounting for changes in grid cell
        % size. We first convert the concentration flux matrix to a mass flux
        % matrix, with each element describing the actual water flux from one grid 
        % cell to another (irrespective of grid size).   
        tm          = A{n}.*ocn_pars.M; % get transport matrix for timestep n

        % This matrix *should* (ideally) have the following properties:
        %   positive-definite
        %   columns should sum to the local water mass
        %       (because each column describes the sink (downstream) distribution for each location)
        %   rows should also sum to the local water mass
        %       (because each row describes the source (upstream) distribution for each location)

    % These properties are not exactly met, for various reasons, 
    % so we apply some corrections...

        % First we get an index of all elements on the diagonal    
        %adindx      = find(eye(size(tm)));
        adindx      = [1:size(tm,1)+1:size(tm,1)*size(tm,1)]'; % faster than find
        % We correct them by subtracting the column sum of the off diagonals from the local volume
        tm(adindx)  = ocn_pars.M - (sum(tm,1)' - tm(adindx));
        % This ensures that the columns of the mass transport matrix sum to exactly
        % the local grid cell volume, such that mass is conserved during transport

        % next we get rid of the negative fluxes by changing their sign and moving them to the transpose
        % (so a negative fux from A to B becomes a positive flux from B to A)

        % first we find subscript indices of negative fluxes
        [i j]       = find(tm<0); 
        % and convert this to a linear index
        indneg      = sub2ind(size(tm),i,j);
        % we place the negative values in their own matrix
        Neg         = spalloc(size(tm,1),size(tm,2),nnz(tm)); % initialise negative matrix
        Neg(indneg) = tm(indneg); % put negatives in negative array

        % we then subtract the negative values from their original locations, 
        % and place equivalent positives in the transpose location
        aa          = tm - Neg + (- Neg');                      % change negatives to positives on transpose
        % we then need to account for these additions to the matrix by adding
        % the equivalent values to the diagonal
        aa(adindx)  = aa(adindx) + sum(Neg,1)' + sum(Neg,2); % add sums of negative rows and columns to diagonal

        % With the matrix corrected, we can extract the surface transport matrix
        bb          = aa(ocn_pars.Ib,ocn_pars.Ib);                    % extract mass transport within surface layer
        bdindx      = find(eye(size(bb)));                      % get coordinates of diagonal

        % In doing this, we removed all vertical fluxes.
        % These must be corrected by assuring the columns correctly sum to the
        % local mass (i.e. vertical fluxes now kept in surface layer). 
        % As above, this is done by subtracting the column sum of the 
        % off-diagonals from the local volume
        bb(bdindx)  = ocn_pars.M(ocn_pars.Ib) + bb(bdindx) - sum(bb,1)'; % correct mass conservation in surface matrix 
        % This keeps unresolved vertical fluxes in surface)

        % Finally, we convert both matrices back to concentration fluxes,
        % and subtract the identity matrix so that the matrix calculation gives
        % the amount of change per time step
        aa          = aa./ocn_pars.M        - speye(size(aa)); % subtract identity matrix for dX/dt
        bb          = bb./ocn_pars.M(ocn_pars.Ib) - speye(size(bb));  % subtract identity matrix for dX/dt

        % Convert from inherent GEnIE timestep (1/96 years) to per day
        tstep_factor = 96/360; % N.B. this will not work if TMs calculated with a different time step!
        % the transport matrices can now be passed back to original arrays
        A{n}        = aa.*tstep_factor;                   
        B{n}        = bb.*tstep_factor;

    end
else % JDW- timestep now in days
    for n=1:numel(A)
        A{n}        = A{n}-speye(size(A{n}));
        B{n}        = 1;
        % Convert from inherent GEnIE timestep (1/96 years) to per day
        tstep_factor = 96/360; % N.B. this will not work if TMs calculated with a different time step!
        % the transport matrices can now be passed back to original arrays
        A{n}        = A{n}.*tstep_factor;
        B{n}        = B{n}.*tstep_factor;
    end
end

    
% save (corrected) transport matrices in ocn_pars
% JDW: OCEAN now contains boundary conditions - should rename to BC or
% something. Everything else is put straight into ocn_pars
%ocn_pars=OCEAN;

ocn_pars.TMA=A ;
ocn_pars.TMB=B'; % BAW why transpose (of cell array)

% JDW: better place for these?
ocn_pars.Idt=1;
ocn_pars.Idt_tot=1;
ocn_pars.Iyr=1;




% find interpolation weighting
% JDW - now done inside dOMGdt
% ocn_pars=gen_fcns.calc_interpolation_weights(ocn_pars,gen_pars.dt);

clear Ii Ib nb

% JDW - what happens to A & B, ocn_pars.A ocn_pars.B
% JDW - clear? aa adindx Neg bbindx tstep_factor indneg idx i j tm
% JDW - all OCEAN is in ocn_pars? yes, so delete OCEAN?!













































