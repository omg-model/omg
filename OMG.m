function [ OUTPUT ] = OMG( runtime , varargin )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Offline Matrix-based GEnIE (OMG) with Matrix Metacommunity Model (MMM)

%% Initialise Matlab path
%OMG(addpath(genpath('.'));
%warning OFF BACKTRACE

%% ------- Initialise Model ------- %%

run('code/OMG/initialise_OMG');

% Start progress timer
parameters = functions.gen_fcns.report_progress( parameters );


%% steady-state solutions
switch parameters.gen_pars.integrate_scheme
    case('newton')
        [TRACERS,ATM]=functions.OMG_fcns.solve_newton(TRACERS,ATM,OCEAN,ECC,parameters,functions,forcings,bioinf,false);
    case('newton_krylov')
        [TRACERS,ATM]=functions.OMG_fcns.solve_newton_krylov(TRACERS,ATM,OCEAN,ECC,parameters,functions,forcings,bioinf,false);
end

%%
%for yr=1:runtime
for    yr=parameters.gen_pars.start_year:runtime+parameters.gen_pars.start_year-1

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SOLVE FOR 1 YEAR  
    [ TRACERS , ATM , TRACERS_t , diagnostics , bioinf ] ...
        = functions.gen_fcns.solve_ODEs(ode_solver,yr,TRACERS,ATM,OCEAN,ECC,parameters,functions,forcings,bioinf);
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% End of year housekeeping

    % collate OUTPUT data (function only evaluated if yr is a save year)
    [OUTPUT] = functions.gen_fcns.collate_output( yr , TRACERS_t , diagnostics , parameters );

    % Write output data to matobj (function only evaluated if yr is a save year)
    functions.gen_fcns.write_output( yr , OUTPUT, parameters , functions );

    % update carbonate chemistry to speed up convergence
    parameters.ocn_pars.T=functions.OMG_fcns.interpolate_var_at_t(360,parameters.ocn_pars,OCEAN.T);
    parameters.ocn_pars.S=functions.OMG_fcns.interpolate_var_at_t(360,parameters.ocn_pars,OCEAN.S);
    [ECC] = functions.gchem_fcns.calc_carbonate_constants(ECC,parameters,functions,true,360);
    [ECC] = functions.gchem_fcns.solve_carbonate_system(ECC,TRACERS,parameters,true);   

    %% report progress through run
    parameters = functions.gen_fcns.report_progress( parameters , yr , TRACERS , inventory , functions );
end


PARTICLES = diagnostics.PARTICLES;
%ATM       = diagnostics.ATM;
%ECC       = diagnostics.ECC;
% GENOME    = diagnostics.genome;
% RGB       = diagnostics.rgb;

%% % Save restart file at end of run
save([outdir '/Restart.mat'],'TRACERS','ATM','PARTICLES','yr');
% save([outdir '/Restart.mat'],'TRACERS','ATM','PARTICLES','ECC','GENOME','RGB','yr');

disp('-----------------------')
disp(['Output data written to ' outdir])
disp('OMG! Finished!')
disp('-----------------------')



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end

