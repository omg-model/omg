%% OMG_Initialise

disp(' ')
disp('-----------------------')
disp('OMG!')
disp('-----------------------')

%% Load default parameters and function handles
run('../general/gen_params');
run('../biogeochemistry/bgc_params');
run('../ecoevo/eco_params');
run('../geochemistry/gchem_params');
run('../OMG/omg_params.m');

%% Overwrite default parameters with any passed into the OMG call
if mod(numel(varargin),2)>0
    error('Unequal number of optional arguments')
end
vararg_counter = 0;
i_string=1:2:numel(varargin)-1; % JDW - i_string?

% Check if user defined ocean configuration file
is_ocnconfig=find(strcmp(varargin(i_string),'ocean_config'));
if ~isempty(is_ocnconfig)
    ocean_config_name = varargin{is_ocnconfig+1}; 
    run(['../../config_files/' ocean_config_name '.m']);
    vararg_counter = vararg_counter + 2;
    disp(['Using ''' ocean_config_name ''' ocean configuration.'])
else
    disp('No ocean user configuration file - using default settings')
end

% JDW - overwrite all default options with user-defined
for n=vararg_counter+1:2:numel(varargin)
    eval([varargin{n} '=varargin{n+1};']);
end

% Check model setup
[bgc_pars,gen_pars,eco_pars,gchem_pars]=gen_fcns.check_model_setup(bgc_pars,gen_pars,eco_pars,gchem_pars,varargin);


%% Check matlab has correct number of workers
p = gcp('nocreate');

%% run initialisation scripts
run('../general/initialise_gen.m')
run('../ocean/initialise_ocn.m')
run('../biogeochemistry/initialise_bgc.m')
run('../ecoevo/initialise_eco.m')
run('../geochemistry/initialise_gchem.m')

%% make SUPER-PARAMETER-STRUCTURE-ARRAY!!!
parameters.gen_pars = gen_pars;
parameters.bgc_pars = bgc_pars;
parameters.eco_pars = eco_pars;
parameters.ocn_pars = ocn_pars;
parameters.ind_pars = I;
parameters.gchem_pars = gchem_pars;

%% make SUPER-FUCNTION-HANDLE-ARRAY!!!
functions.gen_fcns=gen_fcns;
functions.bgc_fcns=bgc_fcns;
if strcmp(bgc_pars.uptake_scheme,'eco')
    functions.eco_fcns=eco_fcns;
end
functions.OMG_fcns=OMG_fcns;
functions.gchem_fcns=gchem_fcns;

%% initialise state
run('../OMG/initialise_state.m')


%% Put ode solver options in structural array
% (done here because it needs jpmax)

% JDW - use parameters.eco_pars.jpmax?
% JDW - why set gen_pars not parameters?
% JDW - eco_pars.jpmax not set if other scheme than 'eco' is used.

parameters.gen_pars.odeoptions=odeset('AbsTol', gen_pars.odeoptions.AbsTol,...
    'RelTol'     , gen_pars.odeoptions.RelTol,...
    'MaxStep'    , gen_pars.odeoptions.MaxStep,...
    'RelTol'     , gen_pars.odeoptions.RelTol,...
    'NormControl', gen_pars.odeoptions.NormControl,...
    'Stats'      , gen_pars.odeoptions.Stats);

% non-negative constraint
if strcmp(bgc_pars.uptake_scheme,'eco')    
    parameters.gen_pars.odeoptions.NonNegative     = ones(1,(eco_pars.jpmax+2).*ocn_pars.nb);
end

%% Initialise mass conservation check
[~ , inventory] = functions.gen_fcns.check_conserve ( TRACERS , 1 , parameters);

%% Saving - maybe a better place for this?

OUTPUT = gen_fcns.reset_output(parameters);

% initialise netcdf output
gen_fcns.create_netcdf_output(parameters,OCEAN,functions);
%gen_fcns.create_netcdf_output(strcat(gen_pars.OMG_output_dir,'/omg_fields_netcdf.nc'),numel(gen_pars.save_output(gen_pars.output_type(:,1))),OCEAN.lon,OCEAN.lat,unique(OCEAN.depth),gen_pars.save_output(gen_pars.output_type(:,1)),bgc_pars);
gen_fcns.initialise_timeseries_output(parameters);

%% Initialise ode solver (Matlab ode or froward euler)

ode_solver = str2func(parameters.gen_pars.integrate_scheme);

if startsWith(parameters.gen_pars.integrate_scheme,'ode')
    
    try % quick test of new ode solver function handle
        ode_solver = str2func(parameters.gen_pars.integrate_scheme);
        [~,~] = ode_solver(@(t,y) t, [0 1], 0);
    catch
        error([func2str(ode_solver) ' is not a valid solver.\n%s'],...
            '(Check parameters.gen_pars.integrate_scheme).')
    end
    % Variable time step solvers cannot use dt_ratio~=1;
    parameters.gen_pars.dt_ratio=1; 
    % Variable time step solvers output n_dt times
    parameters.gen_pars.n_sub_tsteps = 1; % (cannot be modified by n_sub_tsteps)

elseif strcmp(parameters.gen_pars.integrate_scheme,'fwd')

elseif strcmp(parameters.gen_pars.integrate_scheme,'newton')

elseif strcmp(parameters.gen_pars.integrate_scheme,'newton_krylov')

else

    error([parameters.gen_pars.integrate_scheme ' is not a valid solver.\n%s'],...
        '(Check parameters.gen_pars.integrate_scheme).')
end

%% Initialise .mat files
if true 
    [parameters] = gen_fcns.initialise_MatFiles(parameters , outdir , varargin);
end

varargs = varargin;
% save input variables from executable script
save([outdir '/Varargin.mat'],'varargs');
save([outdir '/Parameters.mat'],'parameters');

%%
disp('-----------------------')
disp(['OMG! Running with ' strrep(parameters.gen_pars.integrate_scheme,'_',' ') ' integration scheme...'])
if strcmp(parameters.gen_pars.integrate_scheme,'fwd')
    disp(['and fixed integration timestep of ' num2str(gen_pars.dt./gen_pars.n_sub_tsteps) ' days.'])
end

%% tidy memory from initialisation

% from initialise_OMG.m
clear i_string is_ocnconfig is_ecoconfig ocean_config_name vararg_counter i_string p gen_pars n m d dt varargin OMGcns

% from initialise_gen.m
clear runID prefix dirpath genpars I gen_fcns

% from initialise_bgc.m
clear bgc_pars CaCO3_matrix POM_matrix bgc_fcns

% from initialise_ocn.m
clear Ainput A B aa adindx Neg bbindx tstep_factor indneg idx i j tm ocn_pars nms v_index y

% from initialise_eco.m
clear eco_fcns eco_pars

% from initialise_geochemistry.m
clear gchem_fcns gchem_pars
















