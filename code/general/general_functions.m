function [ gen_fcns ] = general_functions ( )

    gen_fcns.f2v                          = @f2v;
    gen_fcns.v2f                          = @v2f;
    gen_fcns.read_genie_netcdf            = @read_genie_netcdf;
    gen_fcns.read_genie_netcdf_slice	  = @read_genie_netcdf_slice;
    gen_fcns.generate_save_indices        = @generate_save_indices;
    gen_fcns.create_netcdf_output         = @create_netcdf_output;
    gen_fcns.write_netcdf_output          = @write_netcdf_output;
    gen_fcns.initialise_timeseries_output = @initialise_timeseries_output;
    gen_fcns.write_timeseries_output      = @write_timeseries_output;
    gen_fcns.reset_output                 = @reset_output;
    gen_fcns.collate_output               = @collate_output;
    gen_fcns.write_output                 = @write_output;
    gen_fcns.solve_ODEs                   = @solve_ODEs;
    gen_fcns.check_conserve               = @check_conserve;
    gen_fcns.quick_plot                   = @quick_plot;
    gen_fcns.report_progress              = @report_progress;
    gen_fcns.interpolate_TM_vars          = @interpolate_TM_vars;
    gen_fcns.setup_array_indices          = @setup_array_indices;
    gen_fcns.GetFullPath                  = @GetFullPath;
    gen_fcns.load_forcing_data            = @load_forcing_data;
    gen_fcns.calc_interpolation_weights   = @calc_interpolation_weights;
    gen_fcns.check_model_setup            = @check_model_setup;
    gen_fcns.check_restart_setup          = @check_restart_setup;
    gen_fcns.initialise_MatFiles          = @initialise_MatFiles;
    gen_fcns.write_MatFiles               = @write_MatFiles;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
function [vector] = f2v(field,i,j,k)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % f2v
    % Desciption: convert tracer field to vector for use in matrix
    % calculations, based on code from Jake Gebbie.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Nfield = length(i);
    vector = zeros(Nfield,1);
    for ni = 1:Nfield
        vector(ni) = field(k(ni),j(ni),i(ni));
    end

end

%%
function [field] = v2f(vector,i,j,k)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % v2f
    % Desciption: convert tracer vector to field for use in matrix
    % calculations, based on code from Jake Gebbie.
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    NZ = max(k);
    NY = max(j);
    NX = max(i);
    NVEC = length(vector);
    field = nan(NZ,NY,NX);
    for nv = 1:NVEC
        field(k(nv),j(nv),i(nv)) = vector(nv);
    end

end

%%
function [ varargout ] = read_genie_netcdf( path , vec_flag , vec_info , varargin )

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read_genie_netcdf
    % Description: reads netcdf files from GENIE and processes for use with
    % matrices.  Data is last time slice of run, NaNs for missing values and
    % dimensions are permuted (k,j,i).
    %
    % Inputs:
    % path                 - path to folder i.e. exp1/biogem/
    %                        (leave empty if netcdf is sitting in current folder)
    % vec_flag             - flag for 3D/vector output, 1: vector, 0: 3D
    % vec_info             - indexing for vector output as structure
    %                        (i.e. v_index)
    % varargin             - list of variable names to extract from netcdf
    %
    % Outputs:
    % varargin             - names of variable outputs
    %
    % Author: J.D.Wilson 26/03/2014
    %
    % Example:
    %[PO4,SAL,LON,LAT]=read_genie_netcdf('',1,v_index,'ocn_PO4','ocn_sal','lon','lat');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    narg=size(varargin,2);    % number of files to get from netcdf
    
    %newpath=cat(2,path,'fields_biogem_3d.nc');
    newpath=path;
    
    % process and output data
    for n=1:narg
        
        temp_data=extract_data(newpath,varargin{n});
        n_depths=numel(extract_data(newpath,'zt'));
        
        nd=size(temp_data);
        if ndims(temp_data)==4 % if 3D plus time
            temp_data=temp_data(:,:,:,end);             % choose last time slice
            temp_data(temp_data>1e30)=NaN;              % missing values -> NaNs
            temp_data=permute(temp_data,[3 2 1]);   % permute dimensions for vectorising
        elseif ndims(temp_data)==3 % if 2D plus time
            temp_data2=temp_data(:,:,end);             % choose last time slice
            temp_data2(temp_data2>1e30)=NaN;              % missing values -> NaNs
            temp_data2=permute(temp_data2,[2 1]);   % permute dimensions for vectorising
            temp_data=zeros(n_depths,nd(1),nd(2));
            temp_data(1,:,:)=temp_data2; % make 3D for vectorising
        end
        
        % output vector data if selected
        if vec_flag==1
            v_index=vec_info;
            varargout{n}=f2v(temp_data,v_index.i,v_index.j,v_index.rk);
        else
            varargout{n}=temp_data;
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%% nested function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % basic function to read in data from netcdf file
    function [ file_out ] = extract_data ( netcdf_path , variable_name )
        
        
        % netcdf file
        fname=netcdf_path;
        
        % open netcdf file
        ncid=netcdf.open(fname,'NC_NOWRITE');
        
        % retrieve variable data
        file_out=netcdf.getVar(ncid,netcdf.inqVarID(ncid,variable_name));
        
        % close netcdf file
        netcdf.close(ncid);
        
        % end of nested function
    end



end

%%
function [ varargout ] = read_genie_netcdf_slice( path , vec_flag , vec_info , varargin )

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read_genie_netcdf
    % Description: reads netcdf files from GENIE and processes for use with
    % matrices.  Data is last time slice of run, NaNs for missing values and
    % dimensions are permuted (k,j,i).
    %
    % Inputs:
    % path                 - path to folder i.e. exp1/biogem/
    %                        (leave empty if netcdf is sitting in current folder)
    % vec_flag             - flag for 3D/vector output, 1: vector, 0: 3D
    % vec_info             - indexing for vector output as structure
    %                        (i.e. v_index)
    % varargin             - list of variable names to extract from netcdf
    %
    % Outputs:
    % varargin             - names of variable outputs
    %
    % Author: J.D.Wilson 26/03/2014
    %
    % Example vector style:
    %[PO4,SAL,DIC,ALK]=read_genie_netcdf(filenme,1,v_index,'ocn_PO4','ocn_sal','ocn_DIC','ocn_ALK');
    %
    % Example 3D style:
    %[PO4,SAL,LON,LAT]=read_genie_netcdf(filename,0,'','ocn_PO4','ocn_sal','lon','lat');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    narg=size(varargin,2);    % number of files to get from netcdf
    
    %newpath=cat(2,path,'fields_biogem_3d.nc');
    newpath=path;
    
    % process and output data
    for n=1:narg
        
        temp_data=extract_data(newpath,varargin{n});
        temp_data(temp_data>1e30)=NaN;              % missing values -> NaNs
        
        temp_data_copy=temp_data;
        
        % find number of depths
        n_depths=numel(extract_data(newpath,'zt'));
        % find number of timesteps in data
        dtstps = extract_data(newpath,'time');
        nd_tstps = numel(dtstps); % number of time steps in data
        tstps    =       find(dtstps>floor(max(dtstps))) ; % ntimesteps in last year
        ntstps   = numel(find(dtstps>floor(max(dtstps)))); % number of timesteps in last year
        
        % get number of data dimensions
        time_ind=ndims(temp_data);
        % get number of spatial dimensions
        nspdim=ndims(temp_data)-real(nd_tstps>1);
        dimblanks = repmat({':'},1,nspdim); % index for spatial dimensions
        % (temporal dim (last) can be explicit or implicit (if OCEAN.n_A==1)
        for slice=1:ntstps
            
            slicedata=temp_data(dimblanks{:},tstps(slice)); % extract spatial dimensions for time slice
            
            nd=size(temp_data); % get size of spatial grid
            if nspdim==3 % if 3D
                slicedata=permute(slicedata,[3 2 1]); % permute dimensions for vectorising
            elseif nspdim==2 % if 2D
                tmp=permute(slicedata,[2 1]); % permute dimensions for vectorising
                slicedata=zeros(n_depths,nd(1),nd(2));
                slicedata(1,:,:)=tmp; % make 3D for vectorising
            end
            
            % output vector data if selected
            if vec_flag==1
                v_index=vec_info;
                temp_out(:,slice)=f2v(slicedata,v_index.i,v_index.j,v_index.rk);
                varargout{n}=temp_out;
            else
                temp_out(:,:,:,slice)=slicedata;
                varargout{n}=temp_out;
            end
            
            temp_data=temp_data_copy;
        end
        
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%% SUB FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % basic function to read in data from netcdf file
    function [ file_out ] = extract_data ( netcdf_path , variable_name )
        
        
        % netcdf file
        fname=netcdf_path;
        
        % open netcdf file
        ncid=netcdf.open(fname,'NC_NOWRITE');
        
        % retrieve variable data
        file_out=netcdf.getVar(ncid,netcdf.inqVarID(ncid,variable_name));
        
        % close netcdf file
        netcdf.close(ncid);
        
        % end of nested function
    end



end

%%
function [ gen_pars ] = generate_save_indices ( gen_pars )

    % generalising for flexible number of timesteps
    n_dt    = gen_pars.n_dt;
    n_intra = gen_pars.save_intra_annual_n;

    % truncate to within runtime
    timeslice=gen_pars.save_timeslice_output;
    timeslice(timeslice>gen_pars.runtime+gen_pars.start_year)=[]; 

    timeseries=gen_pars.save_timeseries_output;
    timeseries(timeseries>gen_pars.runtime+gen_pars.start_year)=[]; 

    % add final year of specified runtime
    if ~ismember(gen_pars.runtime,gen_pars.save_timeslice_output)
        timeslice=[timeslice gen_pars.runtime];
    end    
    if ~ismember(gen_pars.runtime,gen_pars.save_timeseries_output)
        timeseries=[timeseries gen_pars.runtime];
    end

    % aggregate all output years 
    output=unique([timeslice timeseries])';

    % find number of timesteps for each output plus intra-annual steps
    year_chunks = round(linspace(0,n_dt,n_intra+1)); % divide n_dt into n_intra approximately equal chunks
    tmp = zeros(numel(output)*n_intra,2);
    for n=1:numel(output)
        % get time step index within output year
        tmp2  = (output(n)-1)*n_dt + year_chunks;
        % get consecutive row indices
        i_row = (n-1)*n_intra+1:n*n_intra;
        % collate start and end points for each save interval
        tmp(i_row,:) = [tmp2(1:end-1)+1;tmp2(2:end)]';
    end
    
    % timestep points
    gen_pars.save_ind=tmp;
    % aggregated output timesteps (years)
    gen_pars.save_output=(tmp(:,2)-n_dt/n_intra/2)./n_dt;
    
    % record which type of output is expected at each output point
    tmp2=zeros(numel(gen_pars.save_output),2);

    i_tslice  = 1;
    for n=1:numel(timeslice) 
        tmp2(gen_pars.save_output>timeslice(n)-1 & gen_pars.save_output<timeslice(n),i_tslice)=1; 
    end
    
    i_tseries = 2;
    for n=1:numel(timeseries) 
        tmp2(gen_pars.save_output>timeseries(n)-1 & gen_pars.save_output<timeseries(n),i_tseries)=1; 
    end

    % output types
    gen_pars.output_type=logical(tmp2);

    gen_pars.save_timeseries_output = gen_pars.save_output(gen_pars.output_type(:,i_tseries));
    gen_pars.save_timeslice_output  = gen_pars.save_output(gen_pars.output_type(:,i_tslice ));

    gen_pars.save_years = output;

end


%%
function [ ] = create_netcdf_output ( parameters , OCEAN , functions)

    if parameters.gen_pars.save2netcdf

        fill_val=9.969209968386869E36;          % fill value
    
        filename=strcat(parameters.gen_pars.OMG_output_dir,'/omg_fields_netcdf.nc');
        
        ocn_pars = parameters.ocn_pars;
        gen_pars = parameters.gen_pars;
        bgc_pars = parameters.bgc_pars;
        ind_pars = parameters.ind_pars;
        
        i = ocn_pars.i;
        j = ocn_pars.j;
        k = ocn_pars.rk;

        depth = unique(ocn_pars.depth);
        lat   = ocn_pars.lat;
        lon   = ocn_pars.lon;
        time  = gen_pars.save_timeslice_output; 

        n_intra = gen_pars.save_intra_annual_n;
        time1 = linspace(1/n_intra,1,n_intra)-0.5/n_intra; % times of timeslices in any output year

        n_depth = max(k);
        n_lat   = max(j);
        n_lon   = max(i);  
        n_saves = numel(time); 

        if strcmp(bgc_pars.uptake_scheme,'eco')

            eco_pars = parameters.eco_pars;

            nT1      = eco_pars.nT1;
            nT2      = eco_pars.nT2;
            T1       = eco_pars.T1;
            T2       = eco_pars.T2;
            Tname1   = eco_pars.Tname1;
            Tname2   = eco_pars.Tname2;
            nrgb     = eco_pars.nrgb;
            ngenes   = eco_pars.ngenes;            
        end
        
        % delete file if already exists in output folder
        if exist(filename)==2
            delete(filename)
        end

% create georef variables and attributes    

        nccreate(filename,'time','Dimensions',{'time',n_saves},'Datatype','double','FillValue',fill_val);
        ncwriteatt(filename,'time','standard_name','Time')
        ncwriteatt(filename,'time','units','years')
        ncwriteatt(filename,'time','axis','T')
        ncwrite(filename,'time',time)  

        nccreate(filename,'time1','Dimensions',{'time1',n_intra},'Datatype','double','FillValue',fill_val);
        ncwriteatt(filename,'time1','standard_name','time1')
        ncwriteatt(filename,'time1','long_name','Time in single year')
        ncwriteatt(filename,'time1','units','years')
        ncwriteatt(filename,'time1','axis','T1')
        ncwrite(filename,'time1',time1)
        
        nccreate(filename,'zt','Dimensions',{'zt',n_depth},'Datatype','double','FillValue',fill_val);
        ncwriteatt(filename,'zt','standard_name','Depth')
        ncwriteatt(filename,'zt','units','m')
        ncwriteatt(filename,'zt','axis','Z')
        ncwrite(filename,'zt',depth)

        nccreate(filename,'zt_edges','Dimensions',{'zt_edges',n_depth+1},'Datatype','double','FillValue',fill_val);
        ncwriteatt(filename,'zt_edges','standard_name','Depth')
        ncwriteatt(filename,'zt_edges','units','m')
        ncwriteatt(filename,'zt_edges','axis','Z')
        ncwrite(filename,'zt_edges',ocn_pars.zt_edges)
        
        nccreate(filename,'lon','Dimensions',{'lon',n_lon},'Datatype','double','FillValue',fill_val);
        ncwriteatt(filename,'lon','standard_name','Longitude')
        ncwriteatt(filename,'lon','units','degrees east')
        ncwriteatt(filename,'lon','axis','X')
        ncwrite(filename,'lon',lon)
        
        nccreate(filename,'lat','Dimensions',{'lat',n_lat},'Datatype','double','FillValue',fill_val);
        ncwriteatt(filename,'lat','standard_name','Latitude')
        ncwriteatt(filename,'lat','units','degrees north')
        ncwriteatt(filename,'lat','axis','Y')
        ncwrite(filename,'lat',lat)

        if strcmp(bgc_pars.uptake_scheme,'eco')    

            nccreate(  filename,[Tname1 ' index'],'Dimensions',  {[Tname1 ' index'],nT1},'Datatype','double','FillValue',fill_val);
            ncwriteatt(filename,[Tname1 ' index'],'standard_name',[Tname1 ' index'])
            ncwriteatt(filename,[Tname1 ' index'],'axis',         ['J_' Tname1])
            ncwrite(   filename,[Tname1 ' index'],T1)    

            nccreate(  filename,[Tname2 ' index'],'Dimensions',  {[Tname2 ' index'],nT2},'Datatype','double','FillValue',fill_val);
            ncwriteatt(filename,[Tname2 ' index'],'standard_name',[Tname2 ' index'])
            ncwriteatt(filename,[Tname2 ' index'],'axis',         ['J_' Tname2])
            ncwrite(   filename,[Tname2 ' index'],T2)     
            
%             nccreate(filename,'rgb','Dimensions',{'rgb',nrgb},'Datatype','double','FillValue',fill_val);
%             ncwriteatt(filename,'rgb','standard_name','rgb')
%             ncwriteatt(filename,'rgb','axis','jp')
%             ncwrite(filename,'rgb',1:3)
        end

% create variables for model forcing data

        % Forcing temperature data
        nccreate(filename,'Temperature','Dimensions',{'lon',n_lon,'lat',n_lat,'zt',n_depth,'time1',n_intra},'Datatype','double','FillValue',fill_val);
        ncwriteatt(filename,'Temperature','units','deg C')
        ncwriteatt(filename,'Temperature','long_name','Environmental Forcing: Temperature')

        % Forcing salinity data
        nccreate(filename,'Salinity','Dimensions',{'lon',n_lon,'lat',n_lat,'zt',n_depth,'time1',n_intra},'Datatype','double','FillValue',fill_val);
        ncwriteatt(filename,'Salinity','units','PSU')
        ncwriteatt(filename,'Salinity','long_name','Environmental Forcing: Salinity')
        
        % Forcing seaice data
        nccreate(filename,'Seaice','Dimensions',{'lon',n_lon,'lat',n_lat,'zt',n_depth,'time1',n_intra},'Datatype','double','FillValue',fill_val);
        ncwriteatt(filename,'Seaice','units','%')
        ncwriteatt(filename,'Seaice','long_name','Environmental Forcing: Seaice Percentage Cover')
        
        % Forcing wind speed data
        nccreate(filename,'Windspeed','Dimensions',{'lon',n_lon,'lat',n_lat,'zt',n_depth,'time1',n_intra},'Datatype','double','FillValue',fill_val);
        ncwriteatt(filename,'Windspeed','units','m s-1')
        ncwriteatt(filename,'Windspeed','long_name','Environmental Forcing: Gas exchange windspeed')
        
% write model forcing data

        data = OCEAN.T;
        n_forcing = size(data,2);        
        
        gasXwspeed = sqrt(OCEAN.wspeed./parameters.bgc_pars.gastransfer_a./0.01./8766./OCEAN.seaice);
        
        for t=1:n_intra
            
            % get Forcing data at TimeSlice times
            FORCING.T=functions.OMG_fcns.interpolate_var_at_t(time1(t)*360,ocn_pars,OCEAN.T);
            FORCING.S=functions.OMG_fcns.interpolate_var_at_t(time1(t)*360,ocn_pars,OCEAN.S);
            FORCING.seaice=functions.OMG_fcns.interpolate_var_at_t(time1(t)*360,ocn_pars,OCEAN.seaice);
            FORCING.wspeed=functions.OMG_fcns.interpolate_var_at_t(time1(t)*360,ocn_pars,gasXwspeed);

            
            % Global Temperature at output time t
            data=v2f(FORCING.T,i,j,k);
            data = permute(data,[4 1 2 3]);
            ncwrite(filename,'Temperature',permute(data,[4 3 2 1]),[1 1 1 t]); % add singular first dimension for time = t
        
            % Global Salinity at output time t
            data=v2f(FORCING.S,i,j,k);
            data = permute(data,[4 1 2 3]);
            ncwrite(filename,'Salinity',permute(data,[4 3 2 1]),[1 1 1 t]); % add singular first dimension for time = t
        
            % Global Sea Ice at output time t
            data=v2f(FORCING.seaice,i,j,k);
            data = permute(data,[4 1 2 3]);
            ncwrite(filename,'Seaice',permute(data,[4 3 2 1]),[1 1 1 t]); % add singular first dimension for time = t
        
            % Global Wind Speed at output time t
            data=v2f(FORCING.wspeed,i,j,k);
            data = permute(data,[4 1 2 3]); % add singular first dimension for time = t
            ncwrite(filename,'Windspeed',permute(data,[4 3 2 1]),[1 1 1 t]); % add singular first dimension for time = t
        end


% create variables for model outputs

        % initialise ocean tracers
        for n=1:gen_pars.n_bgc_tracers
            
            var=ind_pars.OCN_names{n};
            name=ind_pars.OCN_long_names{n};
            unit=ind_pars.OCN_units{n};
            
            nccreate(filename,var,'Dimensions',{'lon',n_lon,'lat',n_lat,'zt',n_depth,'time',n_saves},'Datatype','double','FillValue',fill_val);
            ncwriteatt(filename,var,'units',unit)
            ncwriteatt(filename,var,'long_name',name)
               
        end
        
        % initialise sediment/particle tracers
        for n=1:gen_pars.n_particles
            
            var=ind_pars.SED_names{n};
            name=ind_pars.SED_long_names{n};
            unit=ind_pars.SED_units{n};
            
            nccreate(filename,var,'Dimensions',{'lon',n_lon,'lat',n_lat,'zt',n_depth,'time',n_saves},'Datatype','double','FillValue',fill_val);
            ncwriteatt(filename,var,'units',unit)
            ncwriteatt(filename,var,'long_name',name)
               
        end
        
        % initialise carbonate chemistry 
        if gen_pars.save_carbchem
            for n=1:gen_pars.n_carbchem
            
            var=ind_pars.ECC_names{n};
            name=ind_pars.ECC_long_names{n};
            unit=ind_pars.ECC_units{n};
            
            nccreate(filename,var,'Dimensions',{'lon',n_lon,'lat',n_lat,'zt',n_depth,'time',n_saves},'Datatype','double','FillValue',fill_val);
            ncwriteatt(filename,var,'units',unit)
            ncwriteatt(filename,var,'long_name',name)
            end
        end
    
        % initialise ecosystem  
         if strcmp(bgc_pars.uptake_scheme,'eco')
            var ='PHY';
            name='Plankton Community';
            unit='mol kg-1';
            
            nccreate(filename,var,'Dimensions',{'lon',n_lon,'lat',n_lat,[Tname1 ' index'],nT1,[Tname2 ' index'],nT2,'time',n_saves},'Datatype','double','FillValue',fill_val);
            ncwriteatt(filename,var,'units',unit)
            ncwriteatt(filename,var,'long_name',name)
            
            if nrgb>0
                var ='RGB';
                name='Colour Trait';
                unit='-';
                nccreate(filename,var,'Dimensions',{'lon',n_lon,'lat',n_lat,[Tname1 ' index'],nT1,[Tname2 ' index'],nT2,'rgb',nrgb,'time',n_saves},'Datatype','double','FillValue',fill_val);
                ncwriteatt(filename,var,'units',unit)
                ncwriteatt(filename,var,'long_name',name)
            end
            
            if ngenes>0
                var ='GENOME';
                name='Binary genome';
                unit='-';
                nccreate(filename,var,'Dimensions',{'lon',n_lon,'lat',n_lat,[Tname1 ' index'],nT1,[Tname2 ' index'],nT2,'ngenes',ngenes,'time',n_saves},'Datatype','double','FillValue',fill_val);
                ncwriteatt(filename,var,'units',unit)
                ncwriteatt(filename,var,'long_name',name)
            end
        end

    end
    

end

%%
function [ ] = write_netcdf_output (yr,OUTPUT,parameters)

    if parameters.gen_pars.save2netcdf
        
        gen_pars = parameters.gen_pars;
        ocn_pars = parameters.ocn_pars;
        bgc_pars = parameters.bgc_pars;
        I        = parameters.ind_pars;
        % eco_pars = parameters.eco_pars;

        save_timeslice_output   = gen_pars.save_timeslice_output;
        save_timeslice_years    = unique(ceil(gen_pars.save_timeslice_output));
        save_intra_annual_n     = gen_pars.save_intra_annual_n;

        if ismember(yr,save_timeslice_years)

            filename=strcat(parameters.gen_pars.OMG_output_dir,'/omg_fields_netcdf.nc');

            fill_val=9.969209968386869E36;          % fill value

            i=parameters.ocn_pars.i;
            j=parameters.ocn_pars.j;
            k=parameters.ocn_pars.rk;

            save_ind = find(save_timeslice_output<yr & save_timeslice_output>yr-1);

            zero_full = zeros(save_intra_annual_n,max(k),max(j),max(i));

            % write ocean tracers
            for n=1:parameters.gen_pars.n_bgc_tracers

                data_name=parameters.ind_pars.OCN_names{n};
                field = OUTPUT.TimeSlice.(data_name);

                data = zero_full;
                for t=1:save_intra_annual_n
                    data(t,:,:,:)= v2f(field(:,1,t),i,j,k);
                end

                %ncwrite(filename,data_name,data,[save_ind(1) 1 1 1]);
                ncwrite(filename,data_name,permute(data,[4 3 2 1]),[1 1 1 save_ind(1)]);

            end

            % write sediment/particle tracers
            for n=1:parameters.gen_pars.n_particles

                data_name=parameters.ind_pars.SED_names{n};
                field = OUTPUT.TimeSlice.(data_name);

                data = zero_full;
                for t=1:save_intra_annual_n
                    data(t,:,:,:)=v2f(field(:,1,t),i,j,k);
                end

                %ncwrite(filename,data_name,data,[save_ind(1) 1 1 1]);
                ncwrite(filename,data_name,permute(data,[4 3 2 1]),[1 1 1 save_ind(1)]);

            end

            % write carbonate chemistry
            if parameters.gen_pars.save_carbchem
                for n=1:parameters.gen_pars.n_carbchem

                    data_name=parameters.ind_pars.ECC_names{n};
                    field = OUTPUT.TimeSlice.(data_name);

                    data = zero_full;
                    for t=1:save_intra_annual_n
                        data(t,:,:,:)=v2f(field(:,1,t),i,j,k);
                    end
    
                    %ncwrite(filename,data_name,data,[save_ind(1) 1 1 1]);
                    ncwrite(filename,data_name,permute(data,[4 3 2 1]),[1 1 1 save_ind(1)]);

                end
            end

            % write plankton
            if strcmp(parameters.bgc_pars.uptake_scheme,'eco')

                ib      = parameters.ocn_pars.Ib;
                i_sz    = max(i(ib));
                j_sz    = max(j(ib));
                nT1     = parameters.eco_pars.nT1;
                nT2     = parameters.eco_pars.nT2;
                nrgb    = parameters.eco_pars.nrgb;
                ngenes  = parameters.eco_pars.ngenes;

                data_name = 'PHY';
                field = OUTPUT.TimeSlice.(data_name);

                data = zeros(save_intra_annual_n,nT2,nT1,max(j),max(i));
                for t=1:save_intra_annual_n
                    for iT2=1:nT2 % outside loop through trait T2
                        for iT1=1:nT1 % inside loop through trait T1 (nested order of loops not critical)
                            jp = (iT2-1)*nT1+iT1; % increment through T1 then T2 (this is critical)
                            data(t,iT2,iT1,:,:)=v2f(field(:,jp,t),i(ib),j(ib),k(ib));
                        end
                    end
                end
                data(data<0)=0;

                % ncwrite(filename,data_name,data,[save_ind(1) 1 1 1 1]);
                ncwrite(filename,data_name,permute(data,[5 4 3 2 1]),[1 1 1 1 save_ind(1)]);

                if nrgb>0
                    data_name = 'RGB';
                    field = OUTPUT.TimeSlice.(data_name);

                    data = zeros(save_intra_annual_n,nrgb,nT2,nT1,max(j),max(i));
                    for t=1:save_intra_annual_n
                        for iT2=1:nT2 % outside loop through trait T2
                            for iT1=1:nT1 % inside loop through trait T1 (nested order of loops not critical)
                                jp = (iT2-1)*nT1+iT1; % increment through T1 then T2 (this is critical)
                                for i_rgb =1:nrgb
                                    data(t,i_rgb,iT2,iT1,:,:) = v2f(field(:,jp,i_rgb,t),i(ib),j(ib),k(ib));
                                end
                            end
                        end
                    end
                    % ncwrite(filename,data_name,data,[save_ind(1) 1 1 1 1 1]);
                    ncwrite(filename,data_name,permute(data,[6 5 4 3 2 1]),[1 1 1 1 1 save_ind(1)]);
                end

                if ngenes>0
                    data_name = 'GENOME';
                    field = OUTPUT.TimeSlice.(data_name);

                    data = zeros(save_intra_annual_n,ngenes,nT2,nT1,max(j),max(i));
                    for t=1:save_intra_annual_n
                        for iT2=1:nT2 % outside loop through trait T2
                            for iT1=1:nT1 % inside loop through trait T1 (nested order of loops not critical)
                                jp = (iT2-1)*nT1+iT1; % increment through T1 then T2 (this is critical)
                                for i_genes=1:ngenes
                                    data(t,i_genes,iT2,iT1,:,:) = v2f(field(:,jp,i_genes,t),i(ib),j(ib),k(ib));
                                end
                            end
                        end
                    end

                    % ncwrite(filename,data_name,data,[save_ind(1) 1 1 1 1 1]);
                    ncwrite(filename,data_name,permute(data,[6 5 4 3 2 1]),[1 1 1 1 1 save_ind(1)]);
                end

            end
        end  
    end
end


%%
 function [ ] = initialise_timeseries_output ( parameters )

     if parameters.gen_pars.save2text
     
        dir_loc=parameters.gen_pars.OMG_output_dir;
     
        for n=1:parameters.gen_pars.n_bgc_tracers
            
            var=parameters.ind_pars.OCN_names{n};
            name=parameters.ind_pars.OCN_long_names{n};
            unit=parameters.ind_pars.OCN_units{n};
            
            fid=fopen(strcat(dir_loc,'/omg_series_ocn_',var,'.res'),'w');
            fprintf(fid,strcat('%% Time (year) / Inventory (mol) / Global Mean (', unit,') \n'));
            fclose(fid);
            
        end
     
        for n=1:parameters.gen_pars.n_particles
            
            var=parameters.ind_pars.SED_names{n};
            name=parameters.ind_pars.SED_long_names{n};
            unit=parameters.ind_pars.SED_units{n};
            
            fid=fopen(strcat(dir_loc,'/omg_series_ocn_',var,'.res'),'w');
            fprintf(fid,strcat('%% Time (year) / Inventory (mol) / Global Mean (', unit,') \n'));
            fclose(fid);
            
        end
     
        for n=1:parameters.gen_pars.n_carbchem
            
            var=parameters.ind_pars.ECC_names{n};
            name=parameters.ind_pars.ECC_long_names{n};
            unit=parameters.ind_pars.ECC_units{n};
            
            fid=fopen(strcat(dir_loc,'/omg_series_ocn_',var,'.res'),'w');
            fprintf(fid,strcat('%% Time (year) / Inventory (mol) / Global Mean (', unit,') \n'));
            fclose(fid);
            
        end
        
        
        for n=1:parameters.gen_pars.n_atm
            
            var=parameters.ind_pars.ATM_names{n};
            name=parameters.ind_pars.ATM_long_names{n};
            unit=parameters.ind_pars.ATM_units{n};
            
            fid=fopen(strcat(dir_loc,'/omg_series_atm_',var,'.res'),'w');
            fprintf(fid,strcat('%% Time (year)  / Concentration (', unit,') \n'));
            fclose(fid);
            
        end
     end
 end
 
 %%
 
function [ ] = write_timeseries_output (yr,OUTPUT,parameters)

    if parameters.gen_pars.save2text

        gen_pars = parameters.gen_pars;
        ocn_pars = parameters.ocn_pars;
        bgc_pars = parameters.bgc_pars;
        I        = parameters.ind_pars;
        % eco_pars = parameters.eco_pars;

        save_timeseries_output   = gen_pars.save_timeseries_output;
        save_timeseries_years    = unique(ceil(gen_pars.save_timeseries_output));
        save_intra_annual_n     = gen_pars.save_intra_annual_n;

        save_ind = find(save_timeseries_output<yr & save_timeseries_output>yr-1);

        time = save_timeseries_output(save_ind);

        if ismember(yr,save_timeseries_years)

            dir_loc=parameters.gen_pars.OMG_output_dir;

            Ib=parameters.ocn_pars.Ib;
            M=parameters.ocn_pars.M;
            M_total=sum(M);
            A=parameters.ocn_pars.A;

            for n=1:parameters.gen_pars.n_bgc_tracers

                Data_name=parameters.ind_pars.OCN_names{n};
                Write_array = [time OUTPUT.TimeSeries.(Data_name)' OUTPUT.TimeSeries.(Data_name)'./M_total]';

                fid=fopen(strcat(dir_loc,'/omg_series_ocn_',Data_name,'.res'),'a+'); % open with write/append permissions
                fprintf(fid,'%13.5f %17.5e %24.5e \n',Write_array);
                fclose(fid);

            end

            for n=1:parameters.gen_pars.n_particles

                Data_name=parameters.ind_pars.SED_names{n};
                Write_array = [time OUTPUT.TimeSeries.(Data_name)' OUTPUT.TimeSeries.(Data_name)'./M_total]';

                fid=fopen(strcat(dir_loc,'/omg_series_ocn_',Data_name,'.res'),'a+'); % open with write/append permissions
                fprintf(fid,'%13.5f %17.5e %24.5e \n',Write_array);
                fclose(fid);

            end

            for n=1:parameters.gen_pars.n_carbchem

                Data_name=parameters.ind_pars.ECC_names{n};
                Write_array = [time OUTPUT.TimeSeries.(Data_name)' OUTPUT.TimeSeries.(Data_name)'./M_total]';

                fid=fopen(strcat(dir_loc,'/omg_series_ocn_',Data_name,'.res'),'a+'); % open with write/append permissions
                fprintf(fid,'%13.5f %17.5e %24.5e \n',Write_array);
                fclose(fid);

            end

            for n=1:parameters.gen_pars.n_atm

                Data_name=parameters.ind_pars.ATM_names{n};
                Write_array = [time OUTPUT.TimeSeries.(Data_name)']';

                fid=fopen(strcat(dir_loc,'/omg_series_atm_',Data_name,'.res'),'a+'); % open with write/append permissions
                fprintf(fid,'%13.5f %17.5e \n',Write_array);
                fclose(fid);

            end

        end
    end
end


%%
function [ OUTPUT ] = reset_output(parameters)

    I = parameters.ind_pars;
    nb                  = parameters.ocn_pars.nb;
    save_intra_annual_n = parameters.gen_pars.save_intra_annual_n;
    n_dt                = parameters.gen_pars.n_dt;
    n_tracers           = parameters.gen_pars.n_tracers;
    n_bgc_tracers       = parameters.gen_pars.n_bgc_tracers;
    n_particles         = parameters.gen_pars.n_particles;
    n_atm               = parameters.gen_pars.n_atm;
    n_carbchem          = parameters.gen_pars.n_carbchem;

% System state at end of year
    OUTPUT.TRACERS   = zeros(nb,n_tracers);
    if parameters.gen_pars.n_particles>0
        OUTPUT.PARTICLES = zeros(nb,n_particles);
    end
    if parameters.gen_pars.n_atm>0
        OUTPUT.ATM       = zeros(n_atm);
    end
    if parameters.gen_pars.n_carbchem>0
        OUTPUT.ECC=zeros(nb,n_carbchem);
    end
    
% Setup Blank TimeSlices and TimeSeries
    blank_TimeSlice  = zeros(nb,1,save_intra_annual_n);
    blank_TimeSeries = zeros(1,save_intra_annual_n);
    TimeSlice        = struct(I.OCN_names{1},[]);
    TimeSeries       = struct(I.OCN_names{1},[]);

% BGC Tracers
    for i=1:n_bgc_tracers
        TimeSlice  = setfield(TimeSlice,I.OCN_names{i} ,blank_TimeSlice);
        TimeSeries = setfield(TimeSeries,I.OCN_names{i},blank_TimeSeries);
        if i==n_bgc_tracers
            OUTPUT = setfield(OUTPUT,'TimeSlice' ,TimeSlice);
            OUTPUT = setfield(OUTPUT,'TimeSeries',TimeSeries);
        end
    end
% Carbonate Chemistry
    for i=1:n_carbchem
        TimeSlice  = setfield(TimeSlice,I.ECC_names{i} ,blank_TimeSlice);
        TimeSeries = setfield(TimeSeries,I.ECC_names{i},blank_TimeSeries);
        if i==n_carbchem
            OUTPUT = setfield(OUTPUT,'TimeSlice' ,TimeSlice);
            OUTPUT = setfield(OUTPUT,'TimeSeries',TimeSeries);
        end
    end
% Atmospheric
    for i=1:n_atm
        TimeSeries = setfield(TimeSeries,I.ATM_names{i},blank_TimeSeries);
        if i==n_atm
            OUTPUT = setfield(OUTPUT,'TimeSeries',TimeSeries);
        end
    end
% Particles
    for i=1:n_particles
        TimeSlice  = setfield(TimeSlice,I.SED_names{i} ,blank_TimeSlice);
        TimeSeries = setfield(TimeSeries,I.SED_names{i},blank_TimeSeries);
        if i==n_particles
            OUTPUT = setfield(OUTPUT,'TimeSlice' ,TimeSlice);
            OUTPUT = setfield(OUTPUT,'TimeSeries',TimeSeries);
        end
    end

% EcoEvo OUTPUT (as bulk arrays - too many variables to each have named variable)
    if strcmp(parameters.bgc_pars.uptake_scheme,'eco')
        % Setup Global TimeSlices
        OUTPUT.TimeSlice.PHY            = zeros(parameters.ocn_pars.ni ,parameters.eco_pars.jpmax,save_intra_annual_n);
        % Setup Global Time Series
        OUTPUT.TimeSeries.PHY           = zeros(1,save_intra_annual_n);

        % Setup Global Bioinf TimeSlices
        if parameters.eco_pars.nrgb>0
            OUTPUT.TimeSlice.RGB        = zeros(parameters.ocn_pars.ni ,parameters.eco_pars.jpmax,parameters.eco_pars.nrgb,save_intra_annual_n);
        end
        if parameters.eco_pars.ngenes>0
            OUTPUT.TimeSlice.GENOME     = zeros(parameters.ocn_pars.ni ,parameters.eco_pars.jpmax,parameters.eco_pars.ngenes,save_intra_annual_n);
        end

    end
end

%%

function [ OUTPUT ] = collate_output( yr , TRACERS_t , diagnostics , parameters )
    
    
    % if integration year is in gen_pars.save_years index
    if  ismember(yr,parameters.gen_pars.save_years)

        % reset OUTPUT array
        [ OUTPUT ] = reset_output(parameters);

        gen_pars = parameters.gen_pars;
        ocn_pars = parameters.ocn_pars;
        eco_pars = parameters.eco_pars;
        bgc_pars = parameters.bgc_pars;
        I        = parameters.ind_pars;
        loc_M    = parameters.ocn_pars.M;
        loc_A    = parameters.ocn_pars.A;

        % Unpack data
        PARTICLES_t = diagnostics.PARTICLES;
        ATM_t       = diagnostics.ATM      ;
        ECC_t       = diagnostics.ECC      ;
        if eco_pars.nrgb>0 || eco_pars.ngenes>0
           RGB_t    = diagnostics.RGB   ;
           GENOME_t = diagnostics.GENOME;
        end

        save_timeseries_output = gen_pars.save_timeseries_output;
        save_timeslice_output  = gen_pars.save_timeslice_output;
        save_timeseries_years  = unique(ceil(save_timeseries_output));
        save_timeslice_years   = unique(ceil(save_timeslice_output));
        save_intra_annual_n    = gen_pars.save_intra_annual_n;
        n_dt                   = gen_pars.n_dt;

        write_timeslices = ismember(yr,save_timeslice_years);
        write_timeseries = ismember(yr,save_timeseries_years);
        
        year_chunks = round(linspace(0,n_dt,save_intra_annual_n+1));
        save_ind    = [1+year_chunks(1:end-1);year_chunks(2:end)]';

        % System state at end of year
                                   OUTPUT.TRACERS   =   TRACERS_t(  :,:,end);
        if gen_pars.n_particles>0; OUTPUT.PARTICLES = PARTICLES_t(  :,:,end); end
        if gen_pars.n_atm      >0; OUTPUT.ATM       =       ATM_t(    :,end); end
        if gen_pars.n_carbchem >0; OUTPUT.ECC       =       ECC_t(  :,:,end); end
        if eco_pars.nrgb       >0; OUTPUT.RGB       =       RGB_t(:,:,:,end); end
        if eco_pars.ngenes     >0; OUTPUT.GENOME    =    GENOME_t(:,:,:,end); end

        for i=1:save_intra_annual_n

            i_int = [save_ind(i,1):save_ind(i,2)]; 
            
            % Time-averaged TimeSlices and TimeSeries
            
            % BGC Tracers
            for j=1:gen_pars.n_bgc_tracers
                Data_name = I.OCN_names{j};
                Data      = mean(TRACERS_t(:,j,i_int),3);
                % Time Slices
                if write_timeslices
                    F = getfield(OUTPUT.TimeSlice ,Data_name);                   % Extract TimeSlice Fields from OUTPUT structure
                    F(:,:,i) = Data;                                             % Place Data in extracted Fields
                    OUTPUT.TimeSlice  = setfield(OUTPUT.TimeSlice ,Data_name,F); % Place Fields back in OUTPUT structure
                end
                if write_timeseries
                    % Time Series
                    V = getfield(OUTPUT.TimeSeries,Data_name);                   % Extract TimeSeries Vector from OUTPUT structure
                    V(i) = loc_M'*Data;                                          % Integrate Data and place in extracted Vector
                    OUTPUT.TimeSeries = setfield(OUTPUT.TimeSeries,Data_name,V); % Place Vector back in OUTPUT structure
                end
            end
            
            % Particles
            for j=1:gen_pars.n_particles
                Data_name = I.SED_names{j};
                Data      = sum(PARTICLES_t(:,j,i_int)*parameters.gen_pars.dt*parameters.gen_pars.dt_ratio,3,'omitnan'); % n.b. comes in as mol kg-1 day-1 so integrate with bgc timestep
                % Time Slices
                if write_timeslices
                    F = getfield(OUTPUT.TimeSlice ,Data_name);                   % Extract TimeSlice Fields from OUTPUT structure
                    F(:,:,i) = Data;                                             % Place Data in extracted Fields
                    OUTPUT.TimeSlice  = setfield(OUTPUT.TimeSlice ,Data_name,F); % Place Fields back in OUTPUT structure
                end
                if write_timeseries
                    % Time Series
                    V = getfield(OUTPUT.TimeSeries,Data_name);                   % Extract TimeSeries Vector from OUTPUT structure
                    V(i) = loc_M'*Data;                                          % Integrate Data and place in extracted Vector
                    OUTPUT.TimeSeries = setfield(OUTPUT.TimeSeries,Data_name,V); % Place Vector back in OUTPUT structure
                end
            end
            
            % Atmospheric
            if ismember(yr,save_timeseries_years)
                for j=1:gen_pars.n_atm
                    Data_name = I.ATM_names{j};
                    Data      = mean(ATM_t(j,i_int),2);
                    % No Time Slice data for Scalar Atmosphere
                    % Time Series
                    V = getfield(OUTPUT.TimeSeries,Data_name);                   % Extract TimeSeries Vector from OUTPUT structure
                    V(i) = Data;                                          % Integrate average field
                    OUTPUT.TimeSeries = setfield(OUTPUT.TimeSeries,Data_name,V); % Place Vector back in OUTPUT structure
                end
            end
            
            % Carbonate Chemistry
            for j=1:gen_pars.n_carbchem
                Data_name = I.ECC_names{j};
                Data      = mean(ECC_t(:,j,i_int),3);
                % Time Slices
                if write_timeslices
                    F = getfield(OUTPUT.TimeSlice ,Data_name);                   % Extract TimeSlice Fields from OUTPUT structure
                    F(:,:,i) = Data;                                             % Place Data in extracted Fields
                    OUTPUT.TimeSlice  = setfield(OUTPUT.TimeSlice ,Data_name,F); % Place Fields back in OUTPUT structure
                end
                if write_timeseries
                    % Time Series
                    V = getfield(OUTPUT.TimeSeries,Data_name);                   % Extract TimeSeries Vector from OUTPUT structure
                    V(i) = loc_M'*Data;                                          % Integrate Data and place in extracted Vector
                    OUTPUT.TimeSeries = setfield(OUTPUT.TimeSeries,Data_name,V); % Place Vector back in OUTPUT structure
                end
            end

            % Diagnosed Fluxes
            % for j=1:gen_pars.n_bgc_tracers
            %     Data_name = I.OCN_names{j};
            %     Data      = mean(TRACERS_t(:,j,i_int),3);
            %     % Time Slices
            %     if write_timeslices
            %         F = getfield(OUTPUT.TimeSlice ,Data_name);                   % Extract TimeSlice Fields from OUTPUT structure
            %         F(:,:,i) = Data;                                             % Place Data in extracted Fields
            %         OUTPUT.TimeSlice  = setfield(OUTPUT.TimeSlice ,Data_name,F); % Place Fields back in OUTPUT structure
            %     end
            %     if write_timeseries
            %         % Time Series
            %         V = getfield(OUTPUT.TimeSeries,Data_name);                   % Extract TimeSeries Vector from OUTPUT structure
            %         V(i) = loc_M'*Data;                                          % Integrate Data and place in extracted Vector
            %         OUTPUT.TimeSeries = setfield(OUTPUT.TimeSeries,Data_name,V); % Place Vector back in OUTPUT structure
            %     end
            % end

            % EcoEvo OUTPUT
            if strcmp(bgc_pars.uptake_scheme,'eco')
                iphy   = I.PHY;
                Ib     = parameters.ocn_pars.Ib;
%                 isites = gen_pars.tseries_index;
                Data_name = 'PHY';
                Data      = mean( TRACERS_t(Ib,iphy,i_int) ,3);
                % Time Slices
                if write_timeslices
                    F = getfield(OUTPUT.TimeSlice ,Data_name);                   % Extract TimeSlice Fields from OUTPUT structure
                    F(:,:,i) = Data;                                             % Place Data in extracted Fields
                    OUTPUT.TimeSlice  = setfield(OUTPUT.TimeSlice ,Data_name,F); % Place Fields back in OUTPUT structure
                end

                if write_timeseries
                    % Time Series
                    V = getfield(OUTPUT.TimeSeries,Data_name);                   % Extract TimeSeries Vector from OUTPUT structure
                    V(i) = sum(loc_M(Ib)'*Data);                                    % Integrate Data (space and populations) 
                                                                                 % and place in extracted Vector
                    OUTPUT.TimeSeries = setfield(OUTPUT.TimeSeries,Data_name,V); % Place Vector back in OUTPUT structure
                end

                if eco_pars.nrgb>0 & write_timeslices
                    Data_name = 'RGB';
                    Data      = RGB_t(:,:,:,save_ind(i,1));
                    % Place in OUTPUT.TimeSlice_* structure
                    F = getfield(OUTPUT.TimeSlice ,Data_name);
                    F(:,:,:,i) = Data;
                    OUTPUT.TimeSlice  = setfield(OUTPUT.TimeSlice ,Data_name,F);
                end

                if eco_pars.ngenes>0 & write_timeslices
                    Data_name = 'GENOME';
                    Data      = GENOME_t(:,:,:,save_ind(i,1));
                    % Place in OUTPUT.TimeSlice_* structure
                    F = getfield(OUTPUT.TimeSlice ,Data_name);
                    F(:,:,:,i) = Data;
                    OUTPUT.TimeSlice  = setfield(OUTPUT.TimeSlice ,Data_name,F);
                end


            end
        end
    else
        OUTPUT = [];
    end
end


%%
function write_output ( yr , OUTPUT, parameters , functions)

    % if integration year is in gen_pars.save_years index
    if  ismember(yr,parameters.gen_pars.save_years)

        % netcdf timeslices
        functions.gen_fcns.write_netcdf_output(yr,OUTPUT,parameters);
        
        % text timeseries
        functions.gen_fcns.write_timeseries_output (yr,OUTPUT,parameters);
        
        % MatFile Objects
        functions.gen_fcns.write_MatFiles(yr,OUTPUT,parameters);
        
    end

end

%%
function [ TRACERS, ATM , TRACERS_t , diagnostics , bioinf ] = ...
                   solve_ODEs(ode_solver,yr,TRACERS,ATM,OCEAN,ECC,parameters,functions,forcings,bioinf)

    gen_pars = parameters.gen_pars;
    ocn_pars = parameters.ocn_pars;
    eco_pars = parameters.eco_pars;
    
    % generate vector of time steps for Euler solver (and ode solver outputs)
    tsteps       = (yr-1).*360+[0:gen_pars.dt:360];
    nsteps       = gen_pars.n_dt;             % number of full time steps
    dt           = gen_pars.dt;               % duration of full time step
    n_sub_tsteps = gen_pars.n_sub_tsteps;     % number of sub time steps
    sub_dt       = gen_pars.sub_dt;           % duration of sub time step
    
    % reshape state variables arrays to column vector
    y = [reshape(TRACERS,[],1) ; ATM'];

    % initialise y output and diagnostic arrays
    diagnostics.PARTICLES = zeros(ocn_pars.nb,gen_pars.n_particles,nsteps);
    diagnostics.ATM       = zeros(gen_pars.n_atm,nsteps);
    diagnostics.ECC       = zeros(ocn_pars.nb,gen_pars.n_carbchem,nsteps);
    diagnostics.fluxes.uptake = zeros(ocn_pars.nb,gen_pars.n_bgc_tracers,nsteps);
    if eco_pars.nrgb>0 || eco_pars.ngenes>0
        diagnostics.GENOME    = zeros(ocn_pars.ni,eco_pars.jpmax,eco_pars.ngenes,nsteps);
        diagnostics.RGB       = zeros(ocn_pars.ni,eco_pars.jpmax,eco_pars.nrgb,nsteps);
    end
    
    % use Matlab ODE solver if selected
    if startsWith(gen_pars.integrate_scheme,'ode')
        integrate_scheme = 'ode';
        diags_on = false;
    
        % call solver
        [ ~ , yout ] = ode_solver(...
            @(tspan,y) ... % ODE function
            functions.OMG_fcns.dOMGdt(tspan,y,OCEAN,ECC,parameters,functions,forcings,bioinf,diags_on), ...
            tsteps, ...                      % time steps
            y, ...                           % initial conditions
            gen_pars.odeoptions); % ODE  options
    
        % remove initial condition and transpose
        yout = yout(2:end,:)'; 
        
    else % setup forward euler if ode solver not selected
        integrate_scheme = 'fwd';
        yout = zeros(size(y,1),nsteps);
    end


    

    % Do fwd time stepping (as primary solver or for diagnostics)
    for i=1:nsteps

        diags_on = true;

        % loop n_sub_tsteps (n_sub_tsteps = 1 when using ode solvers)
        for i_sub = 1:n_sub_tsteps

            % get time at sub time step
            t = tsteps(i) + sub_dt*i_sub;

            % get gradients and diagnostics
            [dydt,diagnostics_tmp,bioinf] = functions.OMG_fcns.dOMGdt(t,y,OCEAN,ECC,parameters,functions,forcings,bioinf,diags_on);

            switch integrate_scheme
                case 'ode'
                    % get next time step form ode output
                    y = yout(:,i);

                case {'fwd','newton'}
                    % get next time step by time stepping
                    y = y + dydt*sub_dt;
                    
                    if i_sub == n_sub_tsteps
                        % collate output
                        yout(:,i)=y;
                    end

            end
        end

        % collate output
        diagnostics.PARTICLES(:,:,i) = diagnostics_tmp.PARTICLES;
        diagnostics.ATM        (:,i) = diagnostics_tmp.ATM      ;
        diagnostics.ECC      (:,:,i) = diagnostics_tmp.ECC      ;
        diagnostics.fluxes.uptake   (:,:,i) = diagnostics_tmp.fluxes.uptake   ;
        if eco_pars.nrgb>0 || eco_pars.ngenes>0
            diagnostics.RGB    (:,:,:,i) = bioinf.RGB   ;
            diagnostics.GENOME (:,:,:,i) = bioinf.GENOME;
        end

    end
    
    % remove oceanic tracers from yout (ignoring scalar atmospheric fields)
    TRACERS_t = yout(1:ocn_pars.nb*gen_pars.n_tracers,:);
    
    % Reshape yout array for n_grid_cells x n_tracers x n_days
    TRACERS_t = reshape(TRACERS_t , ocn_pars.nb , gen_pars.n_tracers  , nsteps);

    % get state variables from last timestep
    TRACERS   = TRACERS_t(:,:,end);
    ATM       = yout(ocn_pars.nb*gen_pars.n_tracers+1:end,end)'; % extract ATM scalar(s)
    
end


%%
function [conservation_error , inventory] = check_conserve ( TRACERS , inventory_previous , parameters)

    ocn_pars = parameters.ocn_pars;
    I        = parameters.ind_pars;

    % get index of phosphorus tracers
    switch parameters.bgc_pars.uptake_scheme
        case {'MM','restore'}
            iPhosphorus = [I.PO4 I.DOP];
        case 'eco'
            iPhosphorus = [I.PO4 I.DOP I.PHY];
        otherwise
            iPhosphorus = [I.PO4 I.DOP];
    end
    inventory = sum(ocn_pars.M'*TRACERS(:,iPhosphorus));

    conservation_error = 1 - inventory ./ inventory_previous;

end

%%
function [ parameters ] = report_progress( parameters , varargin )

    if parameters.gen_pars.progress_reporting
    
        if nargin > 1 % full call only with additional parameters
            yr        = varargin{1};
            TRACERS   = varargin{2};
            inventory = varargin{3};
            functions = varargin{4};
    
            switch parameters.bgc_pars.uptake_scheme

                case 'eco'

                    % check for mass conservation in inventory
                    [conservation_error , inventory] = functions.gen_fcns.check_conserve ( TRACERS , inventory , parameters);

                    % record time elapsed and reset
                    elapsed = toc(parameters.gen_pars.timer );
                    parameters.gen_pars.timer  = tic;

                    % calculate nnz_out not-functionally-extinct populations
                    nnz_out = nnz(TRACERS(parameters.ocn_pars.Ib,parameters.ind_pars.PHY)>parameters.eco_pars.functional_extinction);
    
                    % display progress statistics
                    disp(sprintf('%5.0f%13.2f%23.3g%14.0f', [yr elapsed conservation_error.*100 nnz_out]))
    
                    % QUICK PLOT OF SURFACE PO4 AND TOTAL PLANKTON BIOMASS
                    functions.gen_fcns.quick_plot( yr , TRACERS , parameters );

                otherwise

                    %frac_complete = yr/(parameters.gen_pars.runtime+parameters.gen_pars.start_year);
                    frac_complete = (yr-parameters.gen_pars.start_year)/parameters.gen_pars.runtime;
                    if frac_complete>=parameters.gen_pars.run_pc_completed

                        % record time elapsed and reset
                        elapsed = toc(parameters.gen_pars.timer );
                        parameters.gen_pars.timer  = tic;

                        % Display percentage completed
                        disp(sprintf('%5.0f%13.2f%12.0f', [yr elapsed frac_complete*100]))
                        
                        % set next threshold for progress report (minimum +10%)
                        next_output = parameters.gen_pars.run_pc_completed + max([0.1 1./(parameters.gen_pars.runtime+parameters.gen_pars.start_year)]);
                        parameters.gen_pars.run_pc_completed = round(next_output*100)/100; % round to nearest percent (trailing errors from 'max')

                    end
            end
        else

            disp('-----------------------')
            switch parameters.bgc_pars.uptake_scheme
                case {'MM','restore'}
                    disp(' Year|Time Elapsed|% Completed')
                case 'eco'
                    disp(' Year|Time Elapsed|Conservation Error (%)|N Populations')
            end

            % Start progress timer
            parameters.gen_pars.timer  = tic; 
        end
    
    end

end


%%
function [ ] = quick_plot ( yr , TRACERS , parameters )

    ocn_pars = parameters.ocn_pars;
    eco_pars = parameters.eco_pars;
    I        = parameters.ind_pars;
    
    nlat = numel(ocn_pars.lat);
    nlon = numel(ocn_pars.lon);
    
    i = ocn_pars.i;
    j = ocn_pars.j;
    k = ocn_pars.k;
    Ib = ocn_pars.Ib;
    
    PO4 = v2f(TRACERS(Ib,I.PO4),i(Ib),j(Ib),k(Ib)).*1e6; % convert from mol/kg --> mmol/m^3
    PO4 = squeeze(PO4(end,:,:));
    
    clf
    
    % Plot surface PO4
    subplot(311)
    imagesc(PO4)
    axis xy
    colorbar
    caxis([0 2.5])
    axis xy
    title(['Surface PO_4 at end of year ' num2str(yr) ' (snapshot, not annual average)'])
    set(gca,'XTick',[],'YTick',[])
    
    
    ind=0;
    PHYSUM=zeros(size(PO4));
    for ii=1:eco_pars.nT2
        for jj=1:eco_pars.nT1
            ind=ind+1;
            tmp = v2f(TRACERS(Ib,I.PHY(ind)),i(Ib),j(Ib),k(Ib)).*1e6; % convert from mol/kg --> mmol/m^3
            PHYcell{ii,jj} = squeeze(tmp(end,:,:));
            PHYSUM = PHYSUM + squeeze(tmp(end,:,:));
        end
    end
    PHY = cell2mat(PHYcell);
    
    
    
    % Plot surface integrated PHY 
    subplot(312)
    PHYSUM(PHYSUM<=0)=NaN;
    imagesc(log10(PHYSUM))
    colorbar
    caxis([-6 0])
    axis xy
    title(['Surface phytoplankton biomass at end of year ' num2str(yr) ' (snapshot, not annual average)'])
    set(gca,'XTick',[],'YTick',[])
    
    % Plot surface PHY by population
    subplot(313)
    PHY(PHY<=0)=NaN;
    imagesc(log10(PHY))
    colorbar
    caxis([log10(eco_pars.functional_extinction) 0])
    axis xy
    set(gca,'XTick',0.5.*nlon:nlon:nlon*(eco_pars.nT1-0.5),...
            'YTick',0.5.*nlat:nlat:nlat*(eco_pars.nT2-0.5),...
            'TickLength',[0 0],...
            'XTickLabel',num2str(eco_pars.T1,'%5.2f'),...
            'YTickLabel',num2str(eco_pars.T2,'%5.2f'),...
            'XColor','k','YColor','k',...
            'GridColor','none','GridAlpha',1)

    hold on
    xt = 0.5:nlon:nlon*eco_pars.nT1+0.5;
    yt = 0.5:nlat:nlat*eco_pars.nT2+0.5;
    plot([xt(1) xt(end)],[yt;yt],'color','w','LineW',1); % plot horizontal grid lines
    plot([xt;xt],[yt(1) yt(end)],'color','w','LineW',1); % plot vertical grid lines

    xlabel(eco_pars.Tname1)
    ylabel(eco_pars.Tname2)
    grid on
    title(['Surface phytoplankton biomass at end of year ' num2str(yr) ' (snapshot, not annual average)'])
    
    colormap(turbo)
    drawnow
end

%%
function [ gen_pars , eco_pars , I ] = setup_array_indices ( gen_pars , bgc_pars , eco_pars )

% Adding tracers:
%   - add a new if bloack containing
%   gen_pars.n_tracers=gen_pars.n_tracers+1;
%   I.NAME=gen_pars.n_tracers;
%   OCN_names{gen_pars.n_tracers}='NAME';
%   OCN_long_names{gen_pars.n_tracers}='LONG NAME';
%   OCN_units{gen_pars.n_tracers}='mol kg-1';

% Setup Ocean Indices

    gen_pars.n_tracers=0; 

    I.OCN_Indices = '------------------------------------------------';
    % I.PO4=1;
    % I.DOP=2;   
    % OCN_names={'PO4','DOP'};
    % OCN_long_names={'Phosphate','Dissolved Organic Phosphorus'};
    % OCN_units={'mol kg-1','mol kg-1'};


    gen_pars.n_tracers=gen_pars.n_tracers+1;        % PO4
    I.PO4=gen_pars.n_tracers;
    OCN_names{gen_pars.n_tracers}='PO4';
    OCN_long_names{gen_pars.n_tracers}='Phosphate';
    OCN_units{gen_pars.n_tracers}='mol kg-1';
    OCN_to_DOM{gen_pars.n_tracers}='DOP';
    OCN_to_POM{gen_pars.n_tracers}='POP';

    gen_pars.n_tracers=gen_pars.n_tracers+1;        % DOP
    I.DOP=gen_pars.n_tracers;
    OCN_names{gen_pars.n_tracers}='DOP';
    OCN_long_names{gen_pars.n_tracers}='Dissolved Organic Phosphorus';
    OCN_units{gen_pars.n_tracers}='mol kg-1';


    if(bgc_pars.O2_select)        
        gen_pars.n_tracers=gen_pars.n_tracers+1;        % O2
        I.O2=gen_pars.n_tracers;   
        OCN_names{gen_pars.n_tracers}='O2';
        OCN_long_names{gen_pars.n_tracers}='Dissolved Oxygen';
        OCN_units{gen_pars.n_tracers}='mol kg-1';
    end    

    if(bgc_pars.CARBCHEM_select)
        gen_pars.n_tracers=gen_pars.n_tracers+1;        % DIC
        I.DIC=gen_pars.n_tracers;
        OCN_names{gen_pars.n_tracers}='DIC';
        OCN_long_names{gen_pars.n_tracers}='Dissolved Inorganic Carbon';
        OCN_units{gen_pars.n_tracers}='mol kg-1';
        OCN_to_DOM{gen_pars.n_tracers}='DOC';
        OCN_to_POM{gen_pars.n_tracers}='POC';

        gen_pars.n_tracers=gen_pars.n_tracers+1;        % ALK
        I.ALK=gen_pars.n_tracers;
        OCN_names{gen_pars.n_tracers}='ALK';
        OCN_long_names{gen_pars.n_tracers}='Alkalinity';
        OCN_units{gen_pars.n_tracers}='mol kg-1';

        gen_pars.n_tracers=gen_pars.n_tracers+1;        % DOC
        I.DOC=gen_pars.n_tracers;
        OCN_names{gen_pars.n_tracers}='DOC';
        OCN_long_names{gen_pars.n_tracers}='Dissolved Organic Carbon';
        OCN_units{gen_pars.n_tracers}='mol kg-1';

    end

    if(bgc_pars.Fe_cycle)
        gen_pars.n_tracers=gen_pars.n_tracers+1;        % Total Fe
        I.TDFe=gen_pars.n_tracers;
        OCN_names{gen_pars.n_tracers}='TDFe';
        OCN_long_names{gen_pars.n_tracers}='Total Iron';
        OCN_units{gen_pars.n_tracers}='mol kg-1';
        OCN_to_DOM{gen_pars.n_tracers}='DOFe';
        OCN_to_POM{gen_pars.n_tracers}='POFe';

        gen_pars.n_tracers=gen_pars.n_tracers+1;        % Total Ligands
        I.TL=gen_pars.n_tracers;
        OCN_names{gen_pars.n_tracers}='TL';
        OCN_long_names{gen_pars.n_tracers}='Total Ligands';
        OCN_units{gen_pars.n_tracers}='mol kg-1';

        gen_pars.n_tracers=gen_pars.n_tracers+1;        % DOFe
        I.DOFe=gen_pars.n_tracers;
        OCN_names{gen_pars.n_tracers}='DOFe';
        OCN_long_names{gen_pars.n_tracers}='Dissolved Organic Fe';
        OCN_units{gen_pars.n_tracers}='mol kg-1';

    end
    
    % record number of biogeochemical tracers
    gen_pars.n_bgc_tracers=gen_pars.n_tracers; 

    if strcmp(bgc_pars.uptake_scheme,'eco')
        eco_pars.jpmax=eco_pars.nsize.*eco_pars.nTopt.*eco_pars.ntroph;
        gen_pars.n_tracers = gen_pars.n_tracers + eco_pars.jpmax;

        I.PHY=(gen_pars.n_bgc_tracers+1):gen_pars.n_tracers;

        OCN_names{gen_pars.n_bgc_tracers+1}='PHY';
        OCN_long_names{gen_pars.n_bgc_tracers+1}=['Plankton Array [' num2str(gen_pars.n_bgc_tracers+1) ':' num2str(gen_pars.n_tracers) ']'];
        OCN_units{gen_pars.n_bgc_tracers+1}='mol kg-1';
    end

    I.OCN_names      = OCN_names;
    I.OCN_long_names = OCN_long_names;
    I.OCN_units      = OCN_units;
    I.OCN_to_DOM     = OCN_to_DOM;
    I.OCN_to_POM     = OCN_to_POM;

% Setup Particulate Indices

    gen_pars.n_particles=1; 
    
    I.SED_Indices = '------------------------------------------------';
    I.POP=1;
    SED_names={'POP'};
    SED_long_names={'Particulate Organic Phosphorus'};
    SED_units={'mol kg-1'};

    if(bgc_pars.CARBCHEM_select)

        gen_pars.n_particles=gen_pars.n_particles+1;    % POC
        I.POC=gen_pars.n_particles;
        SED_names{gen_pars.n_particles}='POC';
        SED_long_names{gen_pars.n_particles}='Particulate Organic Carbon';
        SED_units{gen_pars.n_particles}='mol kg-1';

        gen_pars.n_particles=gen_pars.n_particles+1;    % CaCO3
        I.CaCO3=gen_pars.n_particles;
        SED_names{gen_pars.n_particles}='CaCO3';
        SED_long_names{gen_pars.n_particles}='Calcium Carbonate';
        SED_units{gen_pars.n_particles}='mol kg-1';
    end

    if bgc_pars.Fe_cycle
        gen_pars.n_particles=gen_pars.n_particles+1;    % detritus
        I.Det=gen_pars.n_particles;
        SED_names{gen_pars.n_particles}='Det';
        SED_long_names{gen_pars.n_particles}='Detritus';
        SED_units{gen_pars.n_particles}='mol kg-1';

        gen_pars.n_particles=gen_pars.n_particles+1;    % PFe
        I.POFe=gen_pars.n_particles;
        SED_names{gen_pars.n_particles}='POFe';
        SED_long_names{gen_pars.n_particles}='Particulate Iron';
        SED_units{gen_pars.n_particles}='mol kg-1';
    end

    if gen_pars.n_particles<gen_pars.n_bgc_tracers
        for n=1:gen_pars.n_bgc_tracers-gen_pars.n_particles
            name=['DUM' num2str(n)];
            gen_pars.n_particles=gen_pars.n_particles+1;    % dummy
            eval(['I.' name '=gen_pars.n_particles;']);
            SED_names{gen_pars.n_particles}=name;
            SED_long_names{gen_pars.n_particles}=name;
            SED_units{gen_pars.n_particles}='mol kg-1';
        end
    end


    gen_pars.n_particles=gen_pars.n_bgc_tracers; % due to mapping matrix multiplication, set this equal to tracers. Will work if SED>=OCN!

    I.SED_names      = SED_names;
    I.SED_long_names = SED_long_names;
    I.SED_units      = SED_units;

% Setup Atmospheric Indices

    gen_pars.n_atm=0; % no default atmospheric scalars
    
    I.ATM_Indices = '------------------------------------------------';

    if(bgc_pars.O2_select)
        gen_pars.n_atm=gen_pars.n_atm+1;                % pO2
        I.pO2=gen_pars.n_atm;
        ATM_names{gen_pars.n_atm}='pO2';
        ATM_OCN{gen_pars.n_atm}=I.O2;
        ATM_long_names{gen_pars.n_atm}='Atmospheric Mixing Ratio of Oxygen';
        ATM_units{gen_pars.n_atm}='mol mol-1';
    end

    if(bgc_pars.CARBCHEM_select)
        gen_pars.n_atm=gen_pars.n_atm+1;                % pCO2
        I.pCO2=gen_pars.n_atm;
        ATM_names{gen_pars.n_atm}='pCO2';
        ATM_OCN{gen_pars.n_atm}=I.DIC;
        ATM_long_names{gen_pars.n_atm}='Atmospheric Mixing Ratio of CO2';
        ATM_units{gen_pars.n_atm}='mol mol-1';
    end

    if gen_pars.n_atm>0
        I.ATM_names      = ATM_names;
        I.ATM_long_names = ATM_long_names;
        I.ATM_units      = ATM_units;
        I.ATM_OCN        = ATM_OCN;
    else
        I = rmfield(I,"ATM_Indices");
    end

% Setup Carbonate System Indices
    
    if(bgc_pars.CARBCHEM_select)
        
        I.ECC_Indices = '------------------------------------------------';

        % carbonate chemistry
        gen_pars.n_carbchem=7;                          
        
        I.CO2=1;
        I.HCO3=2;
        I.CO3=3;
        I.PH=4;
        I.SAT_CA=5;
        I.SAT_AR=6;
        I.H=7;

        I.ECC_names{I.CO2}='CO2';
        I.ECC_long_names{I.CO2}='Dissoled Carbon Dioxide';
        I.ECC_units{I.CO2}='mol kg-1';
        
        I.ECC_names{I.HCO3}='HCO3';
        I.ECC_long_names{I.HCO3}='Bicarbonate';
        I.ECC_units{I.HCO3}='mol kg-1';
        
        I.ECC_names{I.CO3}='CO3';
        I.ECC_long_names{I.CO3}='Carbonate';
        I.ECC_units{I.CO3}='mol kg-1';
        
        I.ECC_names{I.PH}='pH';
        I.ECC_long_names{I.PH}='pH';
        I.ECC_units{I.PH}='n/a';
        
        I.ECC_names{I.SAT_CA}='Omega_Ca';
        I.ECC_long_names{I.SAT_CA}='Calcite Saturation State';
        I.ECC_units{I.SAT_CA}='n/a';
        
        I.ECC_names{I.SAT_AR}='Omega_Ar';
        I.ECC_long_names{I.SAT_AR}='Aragonite Saturation State';
        I.ECC_units{I.SAT_AR}='n/a';
        
        I.ECC_names{I.H}='H_ion';
        I.ECC_long_names{I.H}='Hydrogen Ions';
        I.ECC_units{I.H}='mol kg-1';
        
        I.ECC_K_Indices = '------------------------------------------------';

        I.K1=1;
        I.K2=2;
        I.KP1=3;
        I.KP2=4;
        I.KP3=5;
        I.KSi=6;
        I.Kb=7;
        I.K0=8;
        I.Kw=9;
        I.Ks=10;
        I.KF=11;
        I.Kcal=12;
        I.Karg=13;

    else
        gen_pars.n_carbchem=0;  % no carbonate chemistry
    end

end

%%
function File = GetFullPath(File, Style)
    % GetFullPath - Get absolute canonical path of a file or folder
    % Absolute path names are safer than relative paths, when e.g. a GUI or TIMER
    % callback changes the current directory. Only canonical paths without "." and
    % ".." can be recognized uniquely.
    % Long path names (>259 characters) require a magic initial key "\\?\" to be
    % handled by Windows API functions, e.g. for Matlab's FOPEN, DIR and EXIST.
    %
    % FullName = GetFullPath(Name, Style)
    % INPUT:
    %   Name:  String or cell string, absolute or relative name of a file or
    %          folder. The path need not exist. Unicode strings, UNC paths and long
    %          names are supported.
    %   Style: Style of the output as string, optional, default: 'auto'.
    %          'auto': Add '\\?\' or '\\?\UNC\' for long names on demand.
    %          'lean': Magic string is not added.
    %          'fat':  Magic string is added for short names also.
    %          The Style is ignored when not running under Windows.
    %
    % OUTPUT:
    %   FullName: Absolute canonical path name as string or cell string.
    %          For empty strings the current directory is replied.
    %          '\\?\' or '\\?\UNC' is added on demand.
    %
    % NOTE: The M- and the MEX-version create the same results, the faster MEX
    %   function works under Windows only.
    %   Some functions of the Windows-API still do not support long file names.
    %   E.g. the Recycler and the Windows Explorer fail even with the magic '\\?\'
    %   prefix. Some functions of Matlab accept 260 characters (value of MAX_PATH),
    %   some at 259 already. Don't blame me.
    %   The 'fat' style is useful e.g. when Matlab's DIR command is called for a
    %   folder with les than 260 characters, but together with the file name this
    %   limit is exceeded. Then "dir(GetFullPath([folder, '\*.*], 'fat'))" helps.
    %
    % EXAMPLES:
    %   cd(tempdir);                    % Assumed as 'C:\Temp' here
    %   GetFullPath('File.Ext')         % 'C:\Temp\File.Ext'
    %   GetFullPath('..\File.Ext')      % 'C:\File.Ext'
    %   GetFullPath('..\..\File.Ext')   % 'C:\File.Ext'
    %   GetFullPath('.\File.Ext')       % 'C:\Temp\File.Ext'
    %   GetFullPath('*.txt')            % 'C:\Temp\*.txt'
    %   GetFullPath('..')               % 'C:\'
    %   GetFullPath('..\..\..')         % 'C:\'
    %   GetFullPath('Folder\')          % 'C:\Temp\Folder\'
    %   GetFullPath('D:\A\..\B')        % 'D:\B'
    %   GetFullPath('\\Server\Folder\Sub\..\File.ext')
    %                                   % '\\Server\Folder\File.ext'
    %   GetFullPath({'..', 'new'})      % {'C:\', 'C:\Temp\new'}
    %   GetFullPath('.', 'fat')         % '\\?\C:\Temp\File.Ext'
    %
    % COMPILE:
    %   Automatic: InstallMex GetFullPath.c uTest_GetFullPath
    %   Manual:    mex -O GetFullPath.c
    %   Download:  http://www.n-simon.de/mex
    % Run the unit-test uTest_GetFullPath after compiling.
    %
    % Tested: Matlab 6.5, 7.7, 7.8, 7.13, WinXP/32, Win7/64
    %         Compiler: LCC2.4/3.8, BCC5.5, OWC1.8, MSVC2008/2010
    % Assumed Compatibility: higher Matlab versions
    % Author: Jan Simon, Heidelberg, (C) 2009-2016 matlab.2010(a)n(MINUS)simon.de
    %
    % See also: CD, FULLFILE, FILEPARTS.
    
    % $JRev: R-G V:032 Sum:zBDFj0/m8a0f Date:15-Jan-2013 01:06:12 $
    % $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
    % $UnitTest: uTest_GetFullPath $
    % $File: Tools\GLFile\GetFullPath.m $
    % History:
    % 001: 20-Apr-2010 22:28, Successor of Rel2AbsPath.
    % 010: 27-Jul-2008 21:59, Consider leading separator in M-version also.
    % 011: 24-Jan-2011 12:11, Cell strings, '~File' under linux.
    %      Check of input types in the M-version.
    % 015: 31-Mar-2011 10:48, BUGFIX: Accept [] as input as in the Mex version.
    %      Thanks to Jiro Doke, who found this bug by running the test function for
    %      the M-version.
    % 020: 18-Oct-2011 00:57, BUGFIX: Linux version created bad results.
    %      Thanks to Daniel.
    % 024: 10-Dec-2011 14:00, Care for long names under Windows in M-version.
    %      Improved the unittest function for Linux. Thanks to Paul Sexton.
    % 025: 09-Aug-2012 14:00, In MEX: Paths starting with "\\" can be non-UNC.
    %      The former version treated "\\?\C:\<longpath>\file" as UNC path and
    %      replied "\\?\UNC\?\C:\<longpath>\file".
    % 032: 12-Jan-2013 21:16, 'auto', 'lean' and 'fat' style.
    
    % Initialize: ==================================================================
    % Do the work: =================================================================
    
    % #############################################
    % ### USE THE MUCH FASTER MEX ON WINDOWS!!! ###
    % #############################################
    
    % Difference between M- and Mex-version:
    % - Mex does not work under MacOS/Unix.
    % - Mex calls Windows API function GetFullPath.
    % - Mex is much faster.
    
    % Magix prefix for long Windows names:
    if nargin < 2
        Style = 'auto';
    end
    
    % Handle cell strings:
    % NOTE: It is faster to create a function @cell\GetFullPath.m under Linux, but
    % under Windows this would shadow the fast C-Mex.
    if isa(File, 'cell')
        for iC = 1:numel(File)
            File{iC} = GetFullPath(File{iC}, Style);
        end
        return;
    end
    
    % Check this once only:
    isWIN    = strncmpi(computer, 'PC', 2);
    MAX_PATH = 260;
    
    % Warn once per session (disable this under Linux/MacOS):
    persistent hasDataRead
    if isempty(hasDataRead)
        % Test this once only - there is no relation to the existence of DATAREAD!
        %if isWIN
        %   Show a warning, if the slower Matlab version is used - commented, because
        %   this is not a problem and it might be even useful when the MEX-folder is
        %   not inlcuded in the path yet.
        %   warning('JSimon:GetFullPath:NoMex', ...
        %      ['GetFullPath: Using slow Matlab-version instead of fast Mex.', ...
        %       char(10), 'Compile: InstallMex GetFullPath.c']);
        %end
        
        % DATAREAD is deprecated in 2011b, but still available. In Matlab 6.5, REGEXP
        % does not know the 'split' command, therefore DATAREAD is preferred:
        hasDataRead = ~isempty(which('dataread'));
    end
    
    if isempty(File)  % Accept empty matrix as input:
        if ischar(File) || isnumeric(File)
            File = cd;
            return;
        else
            error(['JSimon:', mfilename, ':BadTypeInput1'], ...
                ['*** ', mfilename, ': Input must be a string or cell string']);
        end
    end
    
    if ischar(File) == 0  % Non-empty inputs must be strings
        error(['JSimon:', mfilename, ':BadTypeInput1'], ...
            ['*** ', mfilename, ': Input must be a string or cell string']);
    end
    
    if isWIN  % Windows: --------------------------------------------------------
        FSep = '\';
        File = strrep(File, '/', FSep);
        
        % Remove the magic key on demand, it is appended finally again:
        if strncmp(File, '\\?\', 4)
            if strncmpi(File, '\\?\UNC\', 8)
                File = ['\', File(7:length(File))];  % Two leading backslashes!
            else
                File = File(5:length(File));
            end
        end
        
        isUNC   = strncmp(File, '\\', 2);
        FileLen = length(File);
        if isUNC == 0                        % File is not a UNC path
            % Leading file separator means relative to current drive or base folder:
            ThePath = cd;
            if File(1) == FSep
                if strncmp(ThePath, '\\', 2)   % Current directory is a UNC path
                    sepInd  = strfind(ThePath, '\');
                    ThePath = ThePath(1:sepInd(4));
                else
                    ThePath = ThePath(1:3);     % Drive letter only
                end
            end
            
            if FileLen < 2 || File(2) ~= ':'  % Does not start with drive letter
                if ThePath(length(ThePath)) ~= FSep
                    if File(1) ~= FSep
                        File = [ThePath, FSep, File];
                    else                        % File starts with separator:
                        File = [ThePath, File];
                    end
                else                           % Current path ends with separator:
                    if File(1) ~= FSep
                        File = [ThePath, File];
                    else                        % File starts with separator:
                        ThePath(length(ThePath)) = [];
                        File = [ThePath, File];
                    end
                end
                
            elseif FileLen == 2 && File(2) == ':'   % "C:" current directory on C!
                % "C:" is the current directory on the C-disk, even if the current
                % directory is on another disk! This was ignored in Matlab 6.5, but
                % modern versions considers this strange behaviour.
                if strncmpi(ThePath, File, 2)
                    File = ThePath;
                else
                    try
                        File = cd(cd(File));
                    catch    % No MException to support Matlab6.5...
                        if exist(File, 'dir')  % No idea what could cause an error then!
                            rethrow(lasterror);
                        else  % Reply "K:\" for not existing disk:
                            File = [File, FSep];
                        end
                    end
                end
            end
        end
        
    else         % Linux, MacOS: ---------------------------------------------------
        FSep = '/';
        File = strrep(File, '\', FSep);
        
        if strcmp(File, '~') || strncmp(File, '~/', 2)  % Home directory:
            HomeDir = getenv('HOME');
            if ~isempty(HomeDir)
                File(1) = [];
                File    = [HomeDir, File];
            end
            
        elseif strncmpi(File, FSep, 1) == 0
            % Append relative path to current folder:
            ThePath = cd;
            if ThePath(length(ThePath)) == FSep
                File = [ThePath, File];
            else
                File = [ThePath, FSep, File];
            end
        end
    end
    
    % Care for "\." and "\.." - no efficient algorithm, but the fast Mex is
    % recommended at all!
    if ~isempty(strfind(File, [FSep, '.']))
        if isWIN
            if strncmp(File, '\\', 2)  % UNC path
                index = strfind(File, '\');
                if length(index) < 4    % UNC path without separator after the folder:
                    return;
                end
                Drive            = File(1:index(4));
                File(1:index(4)) = [];
            else
                Drive     = File(1:3);
                File(1:3) = [];
            end
        else  % Unix, MacOS:
            isUNC   = false;
            Drive   = FSep;
            File(1) = [];
        end
        
        hasTrailFSep = (File(length(File)) == FSep);
        if hasTrailFSep
            File(length(File)) = [];
        end
        
        if hasDataRead
            if isWIN  % Need "\\" as separator:
                C = dataread('string', File, '%s', 'delimiter', '\\');  %#ok<REMFF1>
            else
                C = dataread('string', File, '%s', 'delimiter', FSep);  %#ok<REMFF1>
            end
        else  % Use the slower REGEXP, when DATAREAD is not available anymore:
            C = regexp(File, FSep, 'split');
        end
        
        % Remove '\.\' directly without side effects:
        C(strcmp(C, '.')) = [];
        
        % Remove '\..' with the parent recursively:
        R = 1:length(C);
        for dd = reshape(find(strcmp(C, '..')), 1, [])
            index    = find(R == dd);
            R(index) = [];
            if index > 1
                R(index - 1) = [];
            end
        end
        
        if isempty(R)
            File = Drive;
            if isUNC && ~hasTrailFSep
                File(length(File)) = [];
            end
            
        elseif isWIN
            % If you have CStr2String, use the faster:
            %   File = CStr2String(C(R), FSep, hasTrailFSep);
            File = sprintf('%s\\', C{R});
            if hasTrailFSep
                File = [Drive, File];
            else
                File = [Drive, File(1:length(File) - 1)];
            end
            
        else  % Unix:
            File = [Drive, sprintf('%s/', C{R})];
            if ~hasTrailFSep
                File(length(File)) = [];
            end
        end
    end
    
    % "Very" long names under Windows:
    if isWIN
        if ~ischar(Style)
            error(['JSimon:', mfilename, ':BadTypeInput2'], ...
                ['*** ', mfilename, ': Input must be a string or cell string']);
        end
        
        if (strncmpi(Style, 'a', 1) && length(File) >= MAX_PATH) || ...
                strncmpi(Style, 'f', 1)
            % Do not use [isUNC] here, because this concerns the input, which can
            % '.\File', while the current directory is an UNC path.
            if strncmp(File, '\\', 2)  % UNC path
                File = ['\\?\UNC', File(2:end)];
            else
                File = ['\\?\', File];
            end
        end
    end
end

%%
function [ forcings ] = load_forcing_data(gen_pars,bgc_pars,forcings,I,ocn_pars)


    % set default forcing values for each tracer if not already set by user
    for n=1:gen_pars.n_bgc_tracers
        if ~isfield(bgc_pars,['restore_' I.OCN_names{n} '_val']) eval(['bgc_pars.restore_' I.OCN_names{n} '_val= 1.0 ;']); end
        if ~isfield(bgc_pars,['force_' I.OCN_names{n} '_val']) eval(['bgc_pars.force_' I.OCN_names{n} '_val= 1.0 ;']); end
        if ~isfield(bgc_pars,['restore_' I.OCN_names{n} '_timescale']) eval(['bgc_pars.restore_' I.OCN_names{n} '_timescale= 1.0 ;']); end
    end

    for n=1:gen_pars.n_atm
        if ~isfield(bgc_pars,['restore_' I.ATM_names{n} '_val']) eval(['bgc_pars.restore_' I.ATM_names{n} '_val= 1.0 ;']); end
        if ~isfield(bgc_pars,['force_' I.ATM_names{n} '_val']) eval(['bgc_pars.force_' I.ATM_names{n} '_val= 1.0 ;']); end
        if ~isfield(bgc_pars,['restore_' I.ATM_names{n} '_timescale']) eval(['bgc_pars.restore_' I.ATM_names{n} '_timescale= 1.0 ;']); end
    end

    for n=1:gen_pars.n_particles
        if ~isfield(bgc_pars,['restore_' I.SED_names{n} '_val']) eval(['bgc_pars.restore_' I.SED_names{n} '_val= 1.0 ;']); end
        if ~isfield(bgc_pars,['force_' I.SED_names{n} '_val']) eval(['bgc_pars.force_' I.SED_names{n} '_val= 1.0 ;']); end
        if ~isfield(bgc_pars,['restore_' I.SED_names{n} '_timescale']) eval(['bgc_pars.restore_' I.SED_names{n} '_timescale= 1.0 ;']); end
    end

    forcings.ocn.meta=zeros(gen_pars.n_tracers,4);
    forcings.atm.meta=zeros(gen_pars.n_atm,4);
    forcings.sed.meta=zeros(gen_pars.n_particles,4);
    
    % atmosphere arrays
    forcings.atm.interp=cell(gen_pars.n_atm,1);
    
    % ocean arrays
    forcings.ocn.data=zeros(ocn_pars.nb,gen_pars.n_tracers);
    forcings.ocn.interp=cell(gen_pars.n_tracers,1);

    % sediment arrays
    forcings.sed.data=zeros(ocn_pars.nb,gen_pars.n_particles);
    forcings.sed.interp=cell(gen_pars.n_particles,1);

    
    if(~isempty(gen_pars.forcings_directory))
    
        dir_loc=strcat('../../forcings/',gen_pars.forcings_directory,'/');
        
        % load forcing metadata
        if ~isempty(dir(dir_loc))
            
            dir_list=dir(strcat(dir_loc,'*.dat'));
            dir_n=size(dir_list,1);
            
            % load ocean meta data
            for n=1:gen_pars.n_bgc_tracers
                for nn=1:size(dir_list,1)
                    var_name=I.OCN_names{n};
                    force_name=dir_list(nn).name(1:end-4);
                    if strcmp(var_name,force_name)
                        forcings.ocn.meta(n,:)=load(strcat(dir_loc,dir_list(nn).name));
                    end
                end
            end
            
            % load atm meta data
            for n=1:gen_pars.n_atm
                for nn=1:size(dir_list,1)
                    var_name=I.ATM_names{n};
                    force_name=dir_list(nn).name(1:end-4);
                    if strcmp(var_name,force_name)
                        forcings.atm.meta(n,:)=load(strcat(dir_loc,dir_list(nn).name));
                    end
                end
            end

            % load particle meta data
            for n=1:gen_pars.n_particles
                for nn=1:size(dir_list,1)
                    var_name=I.SED_names{n};
                    force_name=dir_list(nn).name(1:end-4);
                    if strcmp(var_name,force_name)
                        forcings.sed.meta(n,:)=load(strcat(dir_loc,dir_list(nn).name));
                    end
                end
            end
            
            dir_list=dir(strcat(dir_loc,'*.sig'));
            dir_n=size(dir_list,1);
            
            % load ocean forcing data
            for n=1:gen_pars.n_bgc_tracers
                for nn=1:size(dir_list,1)
                    var_name=I.OCN_names{n};
                    force_name=dir_list(nn).name(1:end-4);
                    if strcmp(var_name,force_name)
                        tmp=load(strcat(dir_loc,dir_list(nn).name));
                        %tmp(:,1)=tmp(:,1)-gen_pars.start_year; % scale to starting year
                        tmp(:,1)=tmp(:,1)*360; % years to days
                        if tmp(1,1)~=0
                            tmp=[1.0,0.0;tmp]; % add starting point if not specified
                        else
                            tmp(1,1)=1; % otherwise convert to dt
                        end
                        tmp(find(diff(tmp(:,1))==0)+1,1)=tmp(find(diff(tmp(:,1))==0)+1,1)+1; % set split points
                        forcings.ocn.interp{n}=griddedInterpolant(tmp(:,1),tmp(:,2));
                    end


                    % surface explicit
                    if forcings.ocn.meta(n,4) & strcmp(var_name,force_name)
                        tmp=flipud(load(strcat(dir_loc,force_name,'.sur')));
                        tmp2=zeros(max(ocn_pars.k),36,36); tmp2(1,:,:)=tmp;
                        forcings.ocn.data(:,n)=f2v(tmp2,ocn_pars.i,ocn_pars.j,ocn_pars.rk);
                        if forcings.ocn.meta(n,1)
                            % catch case of restoring to zero!
                            ind=forcings.ocn.data(:,n)==0;
                            forcings.ocn.data(ind,n)=NaN;
                        end
                        %forcings.ocn.data(:,n)=forcings.ocn.data(:,n)./sum(forcings.ocn.data(:,n)); % re-scale so force_val is applied evenly
                    end

                    % uniform
                    if forcings.ocn.meta(n,3)
                        forcings.ocn.data(:,n)=1;
                        forcings.ocn.data(:,n)=forcings.ocn.data(:,n)./sum(forcings.ocn.data(:,n)); % re-scale so force_val is applied evenly
                    end

                end
            end
            
            % load atm forcing data
            for n=1:gen_pars.n_atm
                for nn=1:size(dir_list,1)
                    var_name=I.ATM_names{n};
                    force_name=dir_list(nn).name(1:end-4);
                    if strcmp(var_name,force_name)
                        tmp=load(strcat(dir_loc,dir_list(nn).name));
                        %tmp(:,1)=tmp(:,1)-gen_pars.start_year; % scale to starting year
                        tmp(:,1)=tmp(:,1)*360; % years to timesteps
                        if tmp(1,1)~=0
                         tmp=[1.0,0.0;tmp]; % add starting point if not specified
                        else 
                         tmp(1,1)=1; % otherwise convert to dt
                        end
                        tmp(find(diff(tmp(:,1))==0)+1,1)=tmp(find(diff(tmp(:,1))==0)+1,1)+1; % set split points
                        forcings.atm.interp{n}=griddedInterpolant(tmp(:,1),tmp(:,2));
                    end                    
                end                
            end

             % load particle forcing data
            for n=1:gen_pars.n_particles
                for nn=1:size(dir_list,1)
                    var_name=I.SED_names{n};
                    force_name=dir_list(nn).name(1:end-4);
                    if strcmp(var_name,force_name)
                        tmp=load(strcat(dir_loc,dir_list(nn).name));
                        %tmp(:,1)=tmp(:,1)-gen_pars.start_year; % scale to starting year
                        tmp(:,1)=tmp(:,1)*360; % years to days
                        if tmp(1,1)~=0
                            tmp=[1.0,0.0;tmp]; % add starting point if not specified
                        else
                            tmp(1,1)=1; % otherwise convert to dt
                        end
                        tmp(find(diff(tmp(:,1))==0)+1,1)=tmp(find(diff(tmp(:,1))==0)+1,1)+1; % set split points
                        forcings.sed.interp{n}=griddedInterpolant(tmp(:,1),tmp(:,2));
                    end


                    % surface explicit
                    if forcings.sed.meta(n,4) & strcmp(var_name,force_name)
                        tmp=flipud(load(strcat(dir_loc,force_name,'.sur')));
                        tmp2=zeros(max(ocn_pars.k),36,36); tmp2(1,:,:)=tmp;
                        forcings.sed.data(:,n)=f2v(tmp2,ocn_pars.i,ocn_pars.j,ocn_pars.rk);
                        %forcings.sed.data(:,n)=forcings.sed.data(:,n)./sum(forcings.sed.data(:,n)); % re-scale so force_val is applied evenly
                    end

                    % uniform
                    if forcings.sed.meta(n,3)
                        forcings.sed.data(:,n)=1;
                        forcings.sed.data(:,n)=forcings.sed.data(:,n)./sum(forcings.sed.data(:,n)); % re-scale so force_val is applied evenly
                    end

                end
            end
            
        else
            error('Forcing directory selected but no forcing files found')
        end
    end
        
    % copy over user scaling values
    for n=1:gen_pars.n_bgc_tracers
        forcings.ocn.force_scale(1,n)=getfield(bgc_pars,strcat('force_',I.OCN_names{n},'_val'));
        forcings.ocn.restore_scale(1,n)=getfield(bgc_pars,strcat('restore_',I.OCN_names{n},'_val'));
        forcings.ocn.restore_timescale(1,n)=1./getfield(bgc_pars,strcat('restore_',I.OCN_names{n},'_timescale'))/360; % year -> day-1
    end
    for n=1:gen_pars.n_atm
        forcings.atm.force_scale(1,n)=getfield(bgc_pars,strcat('force_',I.ATM_names{n},'_val'));
        forcings.atm.restore_scale(1,n)=getfield(bgc_pars,strcat('restore_',I.ATM_names{n},'_val'));
        forcings.atm.restore_timescale(1,n)=1./getfield(bgc_pars,strcat('restore_',I.ATM_names{n},'_timescale'))/360; % year -> day-1
    end
    for n=1:gen_pars.n_particles
        forcings.sed.force_scale(1,n)=getfield(bgc_pars,strcat('force_',I.SED_names{n},'_val'));
        forcings.sed.restore_scale(1,n)=getfield(bgc_pars,strcat('restore_',I.SED_names{n},'_val'));
        forcings.sed.restore_timescale(1,n)=1./getfield(bgc_pars,strcat('restore_',I.SED_names{n},'_timescale'))/360; % year -> day-1
    end
        
end
 
 
function [ bgc_pars , gen_pars , eco_pars , gchem_pars ] = check_model_setup ( bgc_pars , gen_pars , eco_pars , gchem_pars , user_pars )
 
    %% check user-defined parameters match default parameters
    for n=1:2:numel(user_pars)

        if strcmp(user_pars{n},'ocean_config')
            continue
        end

        ind=strfind(user_pars{n},'.');
        namestr = user_pars{n}(1:ind-1);
        if numel(ind)>1
            valstr  = user_pars{n}(ind(1)+1:ind(2)-1);
        else
            valstr  = user_pars{n}(ind+1:end);
        end

        if ~eval(['isfield(' namestr ',''' valstr ''')'])
            error('>>> ''%s'' is not a recognised parameter <<<',user_pars{n})
        end

    end

    %% specific checks on certain parameters

    % -> gen_pars.n_dt must be exactly divisible by dt_ratio
    if mod(gen_pars.n_dt,gen_pars.dt_ratio)
        error('>>> gen_pars.dt_ratio must go into gen_pars.n_dt exactly <<<')
    end

    if isinteger(gen_pars.n_sub_tsteps)
        error('>>> gen_pars.n_sub_tsteps must be an integer <<<')
    end

    % -> if euler, cannot have more intra_annual saving than there are timesteps
    if strcmp(gen_pars.integrate_scheme,'dd_euler') & gen_pars.save_intra_annual_n>gen_pars.n_dt
        error('>>> gen_pars.save_intra_annual_n must be not be bigger than the number of timesteps: %i <<<',gen_pars.n_dt)
    end
    
    % -> notify user that model will timestep after newton
    %if strcmp(gen_pars.integrate_scheme,'newton') & gen_pars.runtime>1
    %    warning(['>>> OMG will integrate forward in time for ',runtime,' years after the newton solver <<<'])
    %end

   
 
 
 end
    

%%
 
 
function check_restart_setup ( model_pars , OMG_restart_file , outdir)
 

    load(['../../output/' OMG_restart_file '/Varargin.mat']);

    restart_pars = varargs;

    rstrt_names  = restart_pars(1:2:end-1);
    rstrt_fields = restart_pars(2:2:end  );

    mod_names    = model_pars(1:2:end-1);
    mod_fields   = model_pars(2:2:end  );


    %% check model parameters match default restart parameters and warn if not
    nnn=[];mmm=[];
    for n=1:numel(mod_names)

        % do not compare names of restart files
        if ~strcmp(mod_names{n},'gen_pars.OMG_restart_file')

            

            nn = find(ismember(rstrt_names,mod_names{n}));

            if ~isempty(nn)


                r = rstrt_fields(nn);
                m = mod_fields(n);
                

                r=r{1};
                m=m{1};

                nnn=[nnn;nn];
                mmm=[mmm;n];

                if strcmp(mod_names{n},'gen_pars.save_output_directory')
                    if strcmp(OMG_restart_file,m)
                        errorText = sprintf([['Model output directory (gen_pars.save_output_directory) is same as Restart location.\n'],...
                                            ['The new run will overwrite restart. This is not advised.\n'],...
                                            ['If this is desired, need to add an overwrite option in gen_fcns.check_restart_setup.']]);
                        error(errorText)
                    end
                end

                

                if isstr(r) & isstr(m)
                    if ~strcmp(r,m)
                        warnText = sprintf([[mod_names{n} ' does not match between current run and restart\n'],...
                            ['    restart value = ' r '\n'],...
                            ['      model value = ' m]]);
                        warning(warnText)
                    end
                elseif  isnumeric(r) &  isnumeric(m)
                    if ~isequal(r,m)
                        warnText = sprintf([[mod_names{n} ' does not match between current run and restart\n'],...
                            ['    restart value = ' mat2str(r) '\n'],...
                            ['      model value = ' mat2str(m)]]);
                        warning(warnText)
                    end
                elseif isnumeric(r) &  isstr(m)
                    warnText = sprintf([[mod_names{n} ' does not match between current run and restart\n'],...
                        ['    restart value = ' mat2str(r) '\n'],...
                        ['      model value = ' (m)]]);
                    warning(warnText)
                elseif isstr(r) &  isnumeric(m)
                    warnText = sprintf([[mod_names{n} ' does not match between current run and restart\n'],...
                        ['    restart value = ' (r) '\n'],...
                        ['      model value = ' mat2str(m)]]);
                    warning(warnText)
                end
            end
        end
    end



 
end


%%
function [parameters] = initialise_MatFiles( parameters , outdir , varargin)
 
    if parameters.gen_pars.save2matobj
        
        I        = parameters.ind_pars;
        gen_pars = parameters.gen_pars;
    
        % generate OUTPUT array
        OUTPUT = reset_output(parameters);
    
        %% Initialise .mat file
        % matlab.io.MatFile objects allow you to load and save parts of variables
        % in a MAT-file. Working with part of a variable requires less memory
        % than working with its entire contents.
        % Initialise .mat file
    
        % Extract Restart fields and nested TimeSeries and TimeSlice structures from OUTPUT
        names1 = fieldnames(OUTPUT);                % get first level of OUTPUT field names
        for i=1:numel(names1)                       % loop first level of OUTPUT field names
            field = getfield(OUTPUT,names1{i});   % extract current field as its own entity (structure or field)
            if isstruct(field)                    % if this field is a structure ...
                if contains(names1{i},'TimeSeries')
                    TimeSeries = field; % generate new TimeSeries structure
                elseif contains(names1{i},'TimeSlice')
                    TimeSlices = field; % generate new TimeSeries structure
                else
                    error(['Unknown nested structure ' names1{i} ' in OUTPUT structure'])
                end
            end
        end    
    
        % populate TimeSlices structure
        n_output = numel(gen_pars.save_timeslice_output); % get number of output times
        C        = cell(n_output,1);                      % generate empty cells for each output time
        %     C(:)     = {[]};
        names1   = fieldnames(TimeSlices);                % get all field names
        for i=1:numel(names1)                             % loop all field names
            C(:) = {TimeSlices.(names1{i})(:,:,1)};       % Populate C with empty arrays
            TimeSlices.(names1{i}) = C;                   % place cell arrays within TimeSlices structure
        end
        TimeSlices.Times = gen_pars.save_timeslice_output ; % Append output time data
        TimeSlices.parameters = parameters; % Append parameter structure
    
    
        % populate TimeSeries structure
        n_output = numel(gen_pars.save_timeseries_output); % get number of output times
        C        = zeros(n_output,1);                      % generate empty vector with element for each output time
        names1   = fieldnames(TimeSeries);                 % get all field names
        for i=1:numel(names1)                              % loop all field names
            TimeSeries.(names1{i}) = C;                    % place cell array within TimeSeries structure
        end
        TimeSeries.Times = gen_pars.save_timeseries_output; % Append output time data
        TimeSeries.parameters = parameters; % Append parameter structure


        % Initialise Output MatFiles
        save([outdir '/TimeSlices.mat'], '-struct','TimeSlices','-v7.3');
        save([outdir '/TimeSeries.mat'], '-struct','TimeSeries','-v7.3');

        % Create Matfile objects
        parameters.gen_pars.MatObjTimeSlices = matfile([outdir '/TimeSlices.mat'],'Writable', true);
        parameters.gen_pars.MatObjTimeSeries = matfile([outdir '/TimeSeries.mat'],'Writable', true);
    end
 
end
    

%%
function write_MatFiles( yr , OUTPUT, parameters );
    
    if parameters.gen_pars.save2matobj
        
        I        = parameters.ind_pars;
        gen_pars = parameters.gen_pars;
    
        save_timeslice_output   = gen_pars.save_timeslice_output;
        save_timeseries_output  = gen_pars.save_timeseries_output;
    
        save_timeslice_years    = unique(ceil(save_timeslice_output ));
        save_timeseries_years   = unique(ceil(save_timeseries_output));
    
        % if integration year is in save_timeslice_years
        if  ismember(yr,save_timeslice_years)
            MatObjTimeSlices        = parameters.gen_pars.MatObjTimeSlices;
    
            % get index of MatFile output TimeSlice cells during current year
            out_ind = find(save_timeslice_output<yr & save_timeslice_output>yr-1);
    
            % BGC Tracers
            for j=1:parameters.gen_pars.n_bgc_tracers
                % get data name
                Data_name = I.OCN_names{j};
                % Extract Time Slice from OUTPUT structure
                field = OUTPUT.TimeSlice.(Data_name);
                % Extract Time Slice from OUTPUT structure
                fieldcell = mat2cell(field,size(field,1),size(field,2),ones(1,size(field,3)));
                % rearrange dimensions
                fieldcell = permute(fieldcell,[3 1 2]);
                % Place in MatFile
                MatObjTimeSlices.(Data_name)(out_ind,1)  = fieldcell;
            end
    
            % Particles
            for j=1:parameters.gen_pars.n_particles
                Data_name = I.SED_names{j};
                % Extract Time Slice from OUTPUT structure
                field = OUTPUT.TimeSlice.(Data_name);
                % Extract Time Slice from OUTPUT structure
                fieldcell = mat2cell(field,size(field,1),size(field,2),ones(1,size(field,3)));
                % rearrange dimensions
                fieldcell = permute(fieldcell,[3 1 2]);
                % Place in MatFile
                MatObjTimeSlices.(Data_name)(out_ind,1)  = fieldcell;
            end
    
            % Carbonate Chemistry
            for j=1:parameters.gen_pars.n_carbchem
                Data_name = I.ECC_names{j};
                % Extract Time Slice from OUTPUT structure
                field = OUTPUT.TimeSlice.(Data_name);
                % Extract Time Slice from OUTPUT structure
                fieldcell = mat2cell(field,size(field,1),size(field,2),ones(1,size(field,3)));
                % rearrange dimensions
                fieldcell = permute(fieldcell,[3 1 2]);
                % Place in MatFile
                MatObjTimeSlices.(Data_name)(out_ind,1)  = fieldcell;
            end
    
            % EcoEvo
            if strcmp(parameters.bgc_pars.uptake_scheme,'eco')
                Data_name = 'PHY';
                % Extract Time Slice from OUTPUT structure
                field = OUTPUT.TimeSlice.(Data_name);
                % Extract Time Slice from OUTPUT structure
                fieldcell = mat2cell(field,size(field,1),size(field,2),ones(1,size(field,3)));
                % rearrange dimensions
                fieldcell = permute(fieldcell,[3 1 2]);
                % Place in MatFile
                MatObjTimeSlices.(Data_name)(out_ind,1)  = fieldcell;
            end
        end
    
        % if integration year is in save_timeseries_output
        if  ismember(yr,save_timeseries_years)
            MatObjTimeSeries        = parameters.gen_pars.MatObjTimeSeries;
    
            % get index of MatFile output save_timeseries_output points during current year
            out_ind = find(save_timeseries_output<yr & save_timeseries_output>yr-1);
    
            % Write time integrated TimeSeries
            % BGC Tracers
            for j=1:parameters.gen_pars.n_bgc_tracers
                Data_name = I.OCN_names{j};
                % Extract Time Series data from OUTPUT structure
                field = OUTPUT.TimeSeries.(Data_name);
                % Place in MatFile
                MatObjTimeSeries.(Data_name)(out_ind,1)  = field';
            end
    
            % Particles
            for j=1:parameters.gen_pars.n_particles
                Data_name = I.SED_names{j};
                % Extract Time Series data from OUTPUT structure
                field  = OUTPUT.TimeSeries.(Data_name);
                % Place in MatFile
                MatObjTimeSeries.(Data_name)(out_ind,1)  = field';
            end
    
            % Atmosphere
            for j=1:parameters.gen_pars.n_atm
                Data_name = I.ATM_names{j};
                % Extract Time Series data from OUTPUT structure
                field = OUTPUT.TimeSeries.(Data_name);
                % Place in MatFile
                MatObjTimeSeries.(Data_name)(out_ind,1)  = field';
            end
    
            % Carbonate Chemistry
            for j=1:parameters.gen_pars.n_carbchem
                Data_name = I.ECC_names{j};
                % Extract Time Series data from OUTPUT structure
                field = OUTPUT.TimeSeries.(Data_name);
                % Place in MatFile
                MatObjTimeSeries.(Data_name)(out_ind,1)  = field';
            end
    
            % EcoEvo
            if strcmp(parameters.bgc_pars.uptake_scheme,'eco')
                Data_name = 'PHY';
                % Extract Time Slice from OUTPUT structure
                field = OUTPUT.TimeSeries.(Data_name);
                % Place in MatFile
                MatObjTimeSeries.(Data_name)(out_ind,1)  = field';
            end
    
        end
    end
end
























