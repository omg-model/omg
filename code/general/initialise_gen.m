

%% Dependent parameters
% gen_pars.TM_path    = gen_fcns.GetFullPath(gen_pars.TM_path);              % Convert relative to absolute path
[gen_pars,eco_pars,I] = gen_fcns.setup_array_indices(gen_pars,bgc_pars,eco_pars);     % Set array dimension and indices
gen_pars.runtime      = runtime;
gen_pars.dt           = 360/gen_pars.n_dt;                  % base timestep (days)
gen_pars.sub_dt       = gen_pars.dt./gen_pars.n_sub_tsteps; % duration of sub time step

%% Input/Output Stuff
%% Make new overarching output directory if required
if ~exist('../../output','dir')
    !mkdir ../output
end

% JDW - to add user defined output name

dirpath= gen_fcns.GetFullPath('../../output'); % location of output directories
% use user-defined name if present
if ~isempty(gen_pars.save_output_directory)
    outdir=strcat(dirpath,'/',gen_pars.save_output_directory);
else
    % Make unique directory name for each run
    prefix = 'OMG_output'; % beginning of output directory string
    
    % get date (yyyy-mm-dd)
    [y m d]=ymd(datetime);
    dt=[num2str(y,'%i') '-' num2str(m,'%02i') '-' num2str(d,'%02i')];
    % Check what directories already exist for todays date
    D=dir(strcat(dirpath,'/',prefix,['*' dt '*']));
    if isempty(D)
        runID=1; % if none, begin at 1
    else
        runID=str2num(D(end).name(end-3:end))+1; % otherwise begin one above last existing run
    end
    % create new output directory string
    outdir=strcat(dirpath,'/',prefix,'_',dt,'_Run_',num2str(runID,'%04i'));
end

% save output directory in gen_pars
gen_pars.OMG_output_dir=outdir;

if exist(outdir)
    id =  'MATLAB:RMDIR:RemovedFromPath';
    warning('off',id) % turn off annoying warning message!
    rmdir(outdir,'s') % remove pre-existing directory of same name
    warning('on',id)  % turn annoying warning message back on!  
end
mkdir(outdir) % create new output directory

disp(['Output data will be written to: ' gen_fcns.GetFullPath(outdir)])

%% copy config file to output directory
% JDW - this currently doesn't work, relative path is wrong?

mkdir([outdir '/config_file']) % create new config file archive
if ~isempty(is_ocnconfig)
    copyfile(['../../config_files/' ocean_config_name '.m'],[outdir '/config_file/.']);
end
% if ~isempty(is_ecoconfig)
%     copyfile(['../../config_files/' ecosystem_config_name '.m'],[outdir '/config_file/.']);
% end
% copy code archive to output directory
% FOR BEN - this really screws my PC for some reason, it refuses to
% overwrite the output directory!
%mkdir([outdir '/code_archive']) % create new code archive
%copyfile('../../code',[outdir '/code_archive']);
%if isfield(gen_pars,'runscript')
%    copyfile([gen_pars.runscript '.m'],[outdir '/code_archive/run_script.m']);
%end
% fileattrib([outdir '/code_archive/*'],'-w') % make archive read only
% copy matrix variables to output directory

copyfile([gen_pars.TM_path '/matrix_vars.mat'],[outdir])


gen_pars=gen_fcns.generate_save_indices(gen_pars);


%% Misc
% This is set to 0.1 so that progress reported at each 10% of run
gen_pars.run_pc_completed=0.1;



















