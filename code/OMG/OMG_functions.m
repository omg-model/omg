function [ OMG_fcns ] = OMG_functions ( )

    OMG_fcns.dOMGdt            = @dOMGdt;
    OMG_fcns.OMG_transport     = @OMG_transport;
    OMG_fcns.OMG_POMexport     = @OMG_POMexport;
    OMG_fcns.unpack_params     = @unpack_params;
    OMG_fcns.transportmatrix   = @transportmatrix;
    OMG_fcns.transport_weights = @transport_weights; 
    OMG_fcns.interpolate_var_at_t = @interpolate_var_at_t;
    OMG_fcns.solve_newton      = @solve_newton;
    OMG_fcns.solve_newton_krylov = @solve_newton_krylov;
    OMG_fcns.calc_Jpattern     = @calc_Jpattern;
    OMG_fcns.calc_F            = @calc_F;
    OMG_fcns.nsold             = @nsold;
    OMG_fcns.nsoli             = @nsoli;
    OMG_fcns.calc_preconditioner = @calc_preconditioner;

    % includes local subfunctions for nsold and nsoli functions

end

function [ dCdt_out , diagnostics , bioinf] = dOMGdt ( t , y , OCEAN , ECC , parameters , functions , forcings  , bioinf , diags_on)
    % N.B. 'OUTPUT' will not be passed out if dOMGdt is called by ODE solver
    % Any additional outputs will need to be calculated offline.

    [gen_pars, bgc_pars, eco_pars, ocn_pars, I] = unpack_params(parameters);

    Ib            = ocn_pars.Ib;
    
    % Extract state variables and reshape
    ATM     = y(ocn_pars.nb*gen_pars.n_tracers+1:end)'; % extract ATM scalar(s)
    TRACERS = reshape(y(1:ocn_pars.nb*gen_pars.n_tracers) , ocn_pars.nb , gen_pars.n_tracers ); % reshape incoming TRACER array
    SURFACE = TRACERS(ocn_pars.Ib,:); % Get surface TRACERS
    
    % initialise dXdt arrays
    dCdt      = TRACERS.*0;
    dATMdt    = zeros(1,gen_pars.n_atm);
    PARTICLES = nan(ocn_pars.nb,gen_pars.n_particles);
    
    % get time as day-of-year
    tcyc = rem(t,360);
    
    % ocean transport
    [dCdt_transport] = OMG_transport(tcyc,TRACERS,parameters);
    
    % biogeochemistry can be computed in longer timestep (when using fwd solver)
    % skip this section if tcyc goes into gen_pars.dt_ratio.*gen_pars.sub_dt
    % will only skip if using fwd solver
    if ~strcmp(gen_pars.integrate_scheme,'fwd') || mod(tcyc,gen_pars.dt_ratio.*gen_pars.sub_dt)==0
        % this evaluates the following code unless ...
        %     solver is fwd and timestep is divisible by gen_pars.dt*gen_pars.dt_ratio

        % get ocean state at time t
         parameters.ocn_pars.T=interpolate_var_at_t(tcyc,ocn_pars,OCEAN.T);
         parameters.ocn_pars.S=interpolate_var_at_t(tcyc,ocn_pars,OCEAN.S);
         parameters.ocn_pars.solfor=interpolate_var_at_t(tcyc,ocn_pars,OCEAN.solfor);
         parameters.ocn_pars.seaice=interpolate_var_at_t(tcyc,ocn_pars,OCEAN.seaice);
         parameters.ocn_pars.wspeed=interpolate_var_at_t(tcyc,ocn_pars,OCEAN.wspeed);
         parameters.ocn_pars.MLD=interpolate_var_at_t(tcyc,ocn_pars,OCEAN.MLD);
         parameters.ocn_pars.rho=interpolate_var_at_t(tcyc,ocn_pars,OCEAN.rho);
         parameters.ocn_pars.PAR0=interpolate_var_at_t(tcyc,ocn_pars,OCEAN.PAR0);
         parameters.ocn_pars.atm_temp=interpolate_var_at_t(tcyc,ocn_pars,OCEAN.atm_temp);

        % solve equilibrium carbonate chemistry
        [ECC] = functions.gchem_fcns.calc_carbonate_constants(ECC,parameters,functions,false,tcyc);
        [ECC] = functions.gchem_fcns.solve_carbonate_system(ECC,TRACERS,parameters,false);

        % Atmosphere Forcings
        [dATMdt] = functions.bgc_fcns.atm_forcings ( t , ATM , dATMdt , parameters , forcings );

        % Ocean Forcings
        [dCdt] = functions.bgc_fcns.ocn_forcings ( t , TRACERS , dCdt , parameters , forcings );

        % Particle Forcings
        [PARTICLES] = functions.bgc_fcns.sed_forcings ( t , PARTICLES , dCdt , parameters , forcings );

        % Aeolian Fe input
        [ dCdt ] = functions.bgc_fcns.aeolian_Fe ( dCdt , PARTICLES , parameters );

        % Surface Biological PO4 uptake and OM production
        switch bgc_pars.uptake_scheme
            case 'eco'
                [dCdt(Ib,:),POM_prodn,~,ggr] = functions.eco_fcns.ecosystem(SURFACE,parameters);
                % third output variable is 'invfit' 
            otherwise
                % default PO4-based uptake & export
                [dCdt(Ib,:),POM_prodn  ] = functions.bgc_fcns.SurfaceProd(SURFACE,parameters); 
        end

        % remineralise DOM
        [dCdt] = functions.bgc_fcns.remin_DOM (TRACERS, dCdt , parameters );

        % Implicitly remineralise POM
        [dCdt,PARTICLES] = functions.bgc_fcns.remin_POM ( dCdt , POM_prodn , PARTICLES , parameters );

        % Implicitly remineralise CaCO3 
        [dCdt,PARTICLES] = functions.bgc_fcns.remin_CaCO3(dCdt , POM_prodn , ECC , PARTICLES , parameters); 

        % Air-sea gas exchange
        [ dCdt , dATMdt ] = functions.bgc_fcns.airsea_gas_exchange ( dCdt , dATMdt , TRACERS , ATM , ECC , parameters );

        % Restore ocean atmosphere carbon inventory
        [dCdt] = functions.bgc_fcns.restore_ocnatm_Cinv ( dCdt , TRACERS , ATM , parameters);

        % Bioinformatics (if using forward euler - calculated offline if ode solver used)
        if (eco_pars.nrgb>0 || eco_pars.ngenes>0) && diags_on

            % reshape each genome into 1 dimension
            genome =reshape(bioinf.GENOME,ocn_pars.ni*eco_pars.jpmax,eco_pars.ngenes);
            rgb    =reshape(bioinf.RGB   ,ocn_pars.ni*eco_pars.jpmax,eco_pars.nrgb  );

            % Get surface plankton
            Plankton = TRACERS(ocn_pars.Ib,I.PHY); 
            
            % Add random mutation (where biomass exists)
            [genome,rgb] = functions.eco_fcns.gene_mutate(genome,rgb,Plankton,parameters);

            % Get surface transport matrix
            B_t = transportmatrix(tcyc,ocn_pars,ocn_pars.TMB);

            % Exchange 'molecular clock' tracers ...
            [genome,rgb] = functions.eco_fcns.mol_clock(genome,rgb,Plankton,B_t,ggr,eco_pars.Pmut);

            % reshape back to original
            bioinf.GENOME =reshape(genome,ocn_pars.ni,eco_pars.jpmax,eco_pars.ngenes);
            bioinf.RGB    =reshape(rgb   ,ocn_pars.ni,eco_pars.jpmax,eco_pars.nrgb  );
        end
        
    end
    
    % NET
    % (multiply specific terms with timestep ratio here)
    dCdt_out = ...
            dCdt_transport ...
            + dCdt * parameters.gen_pars.dt_ratio ...
            ;
        
    dATMdt_out = ...
        dATMdt * parameters.gen_pars.dt_ratio ...
        ;
        
    % reshape to output    
    dCdt_out = [reshape( dCdt_out , [] , 1 ) ; dATMdt_out'];

    if diags_on
        diagnostics.PARTICLES = PARTICLES;
        diagnostics.ATM       = ATM      ;
        diagnostics.ECC       = ECC.state;
    else
        diagnostics = [];
    end

end
%%

function [dCdt] = OMG_transport(t,TRACERS,parameters)
    
    [gen_pars, bgc_pars, eco_pars, ocn_pars, I] = unpack_params(parameters);

    % ocean transport
    A = ocn_pars.TMA;
    B = ocn_pars.TMB;
    
    dCdt = TRACERS.*0;
   
    if strcmp(bgc_pars.uptake_scheme,'eco')
        % Get time-interpolated full matrix
        A_t = transportmatrix(t,ocn_pars,A);
        switch(eco_pars.plankton_advect)
            case '3D'
                dCdt  = A_t * TRACERS;
                % Plankton biomass fluxes below the surface are diverted to POM
                dCdt(ocn_pars.Ii,I.POM) = dCdt(ocn_pars.Ii,I.POM) + sum(dCdt(ocn_pars.Ii,I.PHY),2);
                dCdt(ocn_pars.Ii,I.PHY) = 0;                
            case '2D'
                dCdt(:,[I.PO4 I.DOP])   = A_t * TRACERS(:,[I.PO4 I.DOP]);
                % Get time-interpolated surface matrix
                B_t = transportmatrix(t,ocn_pars,B);
                dCdt(ocn_pars.Ib,I.PHY) = B_t * TRACERS(ocn_pars.Ib,I.PHY);
            case 'none'
                dCdt(:,[I.PO4 I.DOP])  = A_t * TRACERS(:,[I.PO4 I.DOP]);
        end
    else
        %dCdt  = A_t * TRACERS;
        % calculate temporal weighting
        [TM1, TM2, w1, w2] = transport_weights(t,ocn_pars,A);
        dCdt = (TM1*TRACERS*w1) + (TM2*TRACERS*w2);
    end
    
end
%%

function [gen_pars, bgc_pars, eco_pars, ocn_pars, I] = unpack_params(parameters)

    gen_pars = parameters.gen_pars;
    bgc_pars = parameters.bgc_pars;
    ocn_pars = parameters.ocn_pars;
    I        = parameters.ind_pars;
    if strcmp(parameters.bgc_pars.uptake_scheme,'eco')
        eco_pars = parameters.eco_pars;
    else
        eco_pars.nrgb   = 0;
        eco_pars.ngenes = 0;
    end

end



%%
function [TM_t] = transportmatrix(t,ocn_pars,TM)

    if ~iscell(TM)
            TM_t = TM;
    else
        % calculate temporal weighting
        [TM1 TM2 w1 w2] = transport_weights(t,ocn_pars,TM);
        
        % calculate weighted average (linear interpolation)
%         output = (TM1*input).*w1 + (TM2*input).*w2;
        TM_t = TM1.*w1 + TM2.*w2;
    end
end


%%
function [TM1 TM2 w_before w_after] = transport_weights(t,ocn_pars,TM)
    
    % pad time data and TMs
    TIME = ocn_pars.timesteps;
    TIMEcyclic = [TIME(end)-360 TIME' TIME(1)+360]';
    TMcyclic   = cat(1,TM(end),TM,TM(1));
    
    % find matrices before and after t
    i_before = max(find(TIMEcyclic<=t));
    i_after  = min(find(TIMEcyclic>t));
    
    % get times before and after t
    t_before = TIMEcyclic(i_before);
    t_after  = TIMEcyclic(i_after);
    
    % calculate weights
    w_before = (t-t_before)./(t_after-t_before);
    w_after  = (t_after -t)./(t_after-t_before);
    
    % get time resolved transport matrices
    TM1 = TMcyclic{i_before};
    TM2 = TMcyclic{i_after};
end

%%
function [ interped ] = interpolate_var_at_t(t,ocn_pars,data)

    % pad time
    TIME = ocn_pars.timesteps;
    TIMEcyclic = [TIME(end)-360 TIME' TIME(1)+360]';
    i=[max(numel(TIME)) 1:numel(TIME) 1]';
    
    % find times before and after t
    i_before = find(TIMEcyclic<=t,1,'last');
    i_after  = find(TIMEcyclic>t,1);
    
    % get times before and after t
    t_before = TIMEcyclic(i_before);
    t_after  = TIMEcyclic(i_after);
    
    % calculate weights
    w_before = (t_after -t)./(t_after-t_before);
    w_after  = (t-t_before)./(t_after-t_before);

    interped = data(:,i(i_before))*w_before + data(:,i(i_after))*w_after;

end

%%
function [ TRACERS, ATM ] = solve_newton (TRACERS,ATM,OCEAN,ECC,parameters,functions,forcings,bioinf,diags_on)

    % calculate a sparse pattern for the Jacobian
    Jpattern=calc_Jpattern(parameters);

    % create state vector
    x = [reshape(TRACERS,[],1) ; ATM'];
    pars_in={TRACERS,ATM,OCEAN,ECC,parameters,functions,forcings,bioinf,diags_on};

    % call nsold solver
    %parms=[40, 1, 0, 0]; % straight newton
    parms=[40,1000,0.5,0]; % recommended
    tol=[1e-10,1e-10];
    f=1;
    [xsol, it_hist, ierr] = nsold(x,f,tol,parms,pars_in);

    % output
    if ierr==0
        ATM     = xsol(parameters.ocn_pars.nb*parameters.gen_pars.n_tracers+1:end)';
        TRACERS = reshape(xsol(1:parameters.ocn_pars.nb*parameters.gen_pars.n_tracers) , parameters.ocn_pars.nb , parameters.gen_pars.n_tracers );
    elseif ierr==1
        disp(['termination criterion not met after ' num2str(parms(1)) ' iterations'])
        it_hist
    elseif ierr==2
        disp('failure in the line search.')
        it_hist
    end

end

%%
function [ TRACERS, ATM ] = solve_newton_krylov (TRACERS,ATM,OCEAN,ECC,parameters,functions,forcings,bioinf,diags_on)

    % calculate a sparse pattern for the Jacobian
    %Jpattern=calc_Jpattern(parameters);

    % create state vector
    xinit = [reshape(TRACERS,[],1) ; ATM'];
    parameters.gen_pars.left_precond=true; 
    parameters.gen_pars.update_precond=false;
    pars_in={TRACERS,ATM,OCEAN,ECC,parameters,functions,forcings,bioinf,diags_on};

    % calculate a left preconditioner 
    if parameters.gen_pars.left_precond
        %disp('calculating left preconditioner...')
        % timestepping info
        nsteps   = parameters.gen_pars.n_dt; 
        dt       = parameters.gen_pars.dt;
        tsteps = (1-1).*360+[0:parameters.gen_pars.dt:360];
    
        % get solution at end of period T, e.g., 1 year: x=G(x);
        x=xinit;
        for i=0:nsteps
            t = tsteps(i+1); 
            x=x+functions.OMG_fcns.dOMGdt(t,x,OCEAN,ECC,parameters,functions,forcings,bioinf,diags_on)*dt;
        end

        [parameters] = calc_preconditioner(x,pars_in);
    end
        
    % wrap up OMG variables into one to pass around nsoli subfcns
    pars_in={TRACERS,ATM,OCEAN,ECC,parameters,functions,forcings,bioinf,diags_on};

    % call nsoli solver
    parms = [40,40,0.9,1,20]; % recommended
    tol=[1e-10,1e-10];
    f=1;
    [xsol, it_hist, ierr] = nsoli(xinit,f,tol, parms,pars_in);

    % output
    if ierr==0
        ATM     = xsol(parameters.ocn_pars.nb*parameters.gen_pars.n_tracers+1:end)';
        TRACERS = reshape(xsol(1:parameters.ocn_pars.nb*parameters.gen_pars.n_tracers) , parameters.ocn_pars.nb , parameters.gen_pars.n_tracers );
        disp(['solution found after ' num2str(it_hist(end,2)) ' calls to the model function'])
    elseif ierr==1
        disp(['termination criterion not met after ' num2str(parms(1)) ' iterations'])
        it_hist
    elseif ierr==2
        disp('failure in the line search.')
        it_hist
    end

end


%%
function [Jpattern] = calc_Jpattern (parameters)

    tmp=sparse(parameters.ocn_pars.nb,parameters.ocn_pars.nb);
    Jpattern=sparse((parameters.ocn_pars.nb*parameters.gen_pars.n_tracers+parameters.gen_pars.n_atm),...
        (parameters.ocn_pars.nb*parameters.gen_pars.n_tracers+parameters.gen_pars.n_atm));
    JTM=sparse(parameters.ocn_pars.nb,parameters.ocn_pars.nb);
    Jgasex=sparse(parameters.ocn_pars.nb,1);
    
    % lazy implementation that is not exact but hopefully not far off
    for n=1:parameters.ocn_pars.n_A
        JTM=JTM+parameters.ocn_pars.TMA{n};
    end
    
    Jgasex(parameters.ocn_pars.Ib,1)=1; 
    
    
    tmp=spones(JTM+parameters.bgc_pars.POM_matrix);
    %tmp(:,parameters.ocn_pars.nb+1:end)=Jgasex;
    %tmp(parameters.ocn_pars.nb+1:end,:)=Jgasex';
    
    Jpattern(1:parameters.ocn_pars.nb*parameters.gen_pars.n_tracers,1:parameters.ocn_pars.nb*parameters.gen_pars.n_tracers)=repmat(tmp,parameters.gen_pars.n_tracers,parameters.gen_pars.n_tracers);
    if parameters.gen_pars.n_atm>0
        Jpattern(end,1:parameters.ocn_pars.nb*parameters.gen_pars.n_tracers)=repmat(Jgasex',1,parameters.gen_pars.n_tracers);
        Jpattern(1:parameters.ocn_pars.nb*parameters.gen_pars.n_tracers,end)=repmat(Jgasex,parameters.gen_pars.n_tracers,1);
        Jpattern((parameters.ocn_pars.nb*parameters.gen_pars.n_tracers)+parameters.gen_pars.n_atm,(parameters.ocn_pars.nb*parameters.gen_pars.n_tracers)+parameters.gen_pars.n_atm)=1;
    end

end

%%
function [F,Jacobian,parameters]=calc_F(x,pars_in)

    % JDW: need to consider scaling

    %TRACERS=pars_in{1};
    %ATM=pars_in{2};
    OCEAN=pars_in{3};
    ECC=pars_in{4};
    parameters=pars_in{5};
    functions=pars_in{6};
    forcings=pars_in{7};
    bioinf=pars_in{8};
    diags_on=pars_in{9};

    switch(parameters.gen_pars.integrate_scheme)
        case('newton')
            % F(x) = dxdt = 0
            F=functions.OMG_fcns.dOMGdt(1,x,OCEAN,ECC,parameters,functions,forcings,bioinf,diags_on);
        case('newton_krylov')
            % G(x) = F(x,t+T)
            % left-preconditioning using preconditioner defintion of
            % Khatiwala (2008) Ocean Modelling

            % timestepping info
            nsteps   = parameters.gen_pars.n_dt; 
            dt       = parameters.gen_pars.dt;
            tsteps = (1-1).*360+[0:parameters.gen_pars.dt:360];

            % get solution at end of period T, e.g., 1 year: x=G(x);
            xinit=x;
            for i=0:nsteps
                t = tsteps(i+1); 
                x=x+functions.OMG_fcns.dOMGdt(t,x,OCEAN,ECC,parameters,functions,forcings,bioinf,diags_on)*dt;
            end
            
            % F(x)
            if parameters.gen_pars.left_precond
                %disp(['calculating preconditioned model period ...'])
                z=(parameters.gen_pars.preconditioner.Q * ...
                    (parameters.gen_pars.preconditioner.U \ ...
                    (parameters.gen_pars.preconditioner.L \ ...
                    (parameters.gen_pars.preconditioner.P * ...
                    (parameters.gen_pars.preconditioner.R \ (x-xinit)))))); % (x-xint) = G(x)
                F=parameters.gen_pars.preconditioner.PK*z; % preconditioned G
            else
                %disp(['calculating model period ...'])
                F=x-xinit;
            end

    end

    % Optional Jacobian output using numjac() for nsold.m
    if nargout>1

        FAC=[];
        G=[];
        t=1;

        Jpattern=calc_Jpattern(parameters);
        [Jacobian,FAC,G] = numjac(@(t,x) functions.OMG_fcns.dOMGdt(t,x,OCEAN,ECC,parameters,functions,forcings,bioinf,diags_on),t,x,F,ones(numel(x),1)*1e-10,FAC,false,Jpattern,G);

    end

end

%%
function [parameters] = calc_preconditioner(x,pars_in)

    parameters=pars_in{5};


        %TRACERS=pars_in{1};
        %ATM=pars_in{2};
        OCEAN=pars_in{3};
        ECC=pars_in{4};
        
        functions=pars_in{6};
        forcings=pars_in{7};
        bioinf=pars_in{8};
        diags_on=pars_in{9};
        
        % timestepping info
        nsteps   = parameters.gen_pars.n_dt; 
        dt       = parameters.gen_pars.dt;

        FAC=[];
        G=[];
        t=180;
        
        f=functions.OMG_fcns.dOMGdt(t,x,OCEAN,ECC,parameters,functions,forcings,bioinf,diags_on);
        
        Jpattern=calc_Jpattern(parameters);
        [f_prime,FAC,G] = numjac(@(t,x) functions.OMG_fcns.dOMGdt(t,x,OCEAN,ECC,parameters,functions,forcings,bioinf,diags_on),t,x,f,ones(numel(x),1)*1e-10,FAC,false,Jpattern,G);

        f_prime=(dt*nsteps)*f_prime; % ∆t*f' where ∆t=360 days
        [L,U,P,Q,R]=lu(f_prime);

        I=speye(size(f_prime));
        PK=I-f_prime;

        parameters.gen_pars.preconditioner.L=L;
        parameters.gen_pars.preconditioner.U=U;
        parameters.gen_pars.preconditioner.P=P;
        parameters.gen_pars.preconditioner.Q=Q;
        parameters.gen_pars.preconditioner.R=R;
        parameters.gen_pars.preconditioner.PK=PK;


end

%%
function [sol, it_hist, ierr, x_hist] = nsold(x,f,tol,parms,pars_in)
% NSOLD  Newton-Armijo nonlinear solver
%
% Factor Jacobians with Gaussian Elimination
%
% Hybrid of Newton, Shamanskii, Chord
%
% C. T. Kelley, April 1, 2003.
%
% This code comes with no guarantee or warranty of any kind.
%
% function [sol, it_hist, ierr, x_hist] = nsold(x,f,tol,parms)
%
% inputs:
%        initial iterate = x
%        function = f
%        tol = [atol, rtol] relative/absolute
%                           error tolerances
%        parms = [maxit, isham, rsham, jdiff, nl, nu]
%        maxit = maxmium number of iterations
%                default = 40
%        isham, rsham: The Jacobian matrix is
%                computed and factored after isham
%                updates of x or whenever the ratio
%                of successive l2 norms of the
%                 nonlinear residual exceeds rsham.
%
%            isham = -1, rsham = .5 is the default
%            isham =  1, rsham = 0 is Newton's method,
%            isham = -1, rsham = 1 is the chord method,
%            isham =  m, rsham = 1 is the Shamanskii method with
%                        m steps per Jacobian evaluation
%
%                       The Jacobian is computed and factored
%                       whenever the stepsize
%                       is reduced in the line search.
%
%       jdiff = 1: compute Jacobians with forward differences
%       jdiff = 0: a call to f will provide analytic Jacobians
%                         using the syntax [function,jacobian] = f(x)
%                 defaults = [40, 1000, .5, 1]
%
%       nl, nu: lower and upper bandwidths of a banded Jacobian.
%               If you include nl and nu in the parameter list,
%               the Jacobian will be evaluated with a banded differencing
%               scheme and stored as a sparse matrix.
%
%
%
% output:
%    sol = solution
%    it_hist = array of iteration history, useful for tables and plots
%                The two columns are the residual norm and
%                number of step size reductions done in the line search.
%
%        ierr = 0 upon successful termination
%        ierr = 1 if after maxit iterations
%             the termination criterion is not satsified
%        ierr = 2 failure in the line search. The iteration
%             is terminated if too many steplength reductions
%             are taken.
%
%    x_hist = matrix of the entire interation history.
%             The columns are the nonlinear iterates. This
%             is useful for making movies, for example, but
%             can consume way too much storage. This is an
%             OPTIONAL argument. Storage is only allocated
%             if x_hist is in the output argument list.
%
%
% internal parameter:
%       debug = turns on/off iteration statistics display as
%               the iteration progresses
%
% Here is an example. The example computes pi as a root of sin(x)
% with Newton's method, forward difference derivatives,
% and plots the iteration history. Note that x_hist is not in
% the output list.
%
%
%  x = 3; tol = [1.d-6, 1.d-6]; params = [40, 1, 0];
%  [result, errs, ierr] = nsold(x, 'sin', tol, params);
%  result
%  semilogy(errs)
%
%
% Set the debug parameter, 1 turns display on, otherwise off.
%
debug = 0;
%
% Initialize it_hist, ierr, and set the iteration parameters.
%
ierr = 0;
maxarm = 20;
maxit = 40;
isham = -1;
rsham = .5;
jdiff = 1;
iband = 0;
if nargin >= 4 & length(parms) ~= 0
    maxit = parms(1); isham = parms(2); rsham = parms(3); 
        if length(parms) >= 4
            jdiff = parms(4);
        end
        if length(parms) >= 6
            nl = parms(5); nu = parms(6);
            iband = 1;
        end
    end
rtol = tol(2); atol = tol(1);
it_hist = [];
n = length(x);
if nargout == 4, x_hist = x; end
fnrm = 1;
itc = 0;
%
% evaluate f at the initial iterate
% compute the stop tolerance
%
%f0 = feval(f,x); jdw
f0 = calc_F(x,pars_in);
fnrm = norm(f0);
it_hist = [fnrm,0];
fnrmo = 1;
itsham = isham;
stop_tol = atol+rtol*fnrm;
%
% main iteration loop
%
while(fnrm > stop_tol & itc < maxit)
%
% keep track of the ratio (rat = fnrm/frnmo)
% of successive residual norms and 
% the iteration counter (itc)
%
    rat = fnrm/fnrmo;
    outstat(itc+1, :) = [itc fnrm rat];
    fnrmo = fnrm; 
    itc = itc+1;
%
% evaluate and factor the Jacobian
% on the first iteration, every isham iterates, or
% if the ratio of successive residual norm is too large
%
    if(itc == 1 | rat > rsham | itsham == 0 | armflag == 1)
        itsham = isham;
    jac_age = -1;
    if jdiff == 1 
            if iband == 0
            [l, u] = diffjac(x,f,f0);
            else
            jacb = bandjac(f,x,f0,nl,nu); 
            %[l,u] = lu(jacb); jdw
            [L,U,P,Q,R]=lu(jac);
            end
        else
            %[fv,jac] = feval(f,x); jdw
            [fv,jac] = calc_F(x,pars_in);
        %[l,u] = lu(jac); jdw
        [L,U,P,Q,R]=lu(jac);
        end
    end
    itsham = itsham-1;
%
% compute the Newton direction 
%
    %tmp = -l\f0;
    %direction = u\tmp; jdw
    direction=  -(Q * (U \ (L \ (P * (R \ f0)))));
%
% Add one to the age of the Jacobian after the factors have been
% used in a solve. A fresh Jacobian has an age of -1 at birth.
%
    jac_age = jac_age+1;
    xold = x; fold = f0;
    [step,iarm,x,f0,armflag] = armijo(direction,x,f0,f,maxarm,pars_in);
%
% If the line search fails and the Jacobian is old, update it.
% If the Jacobian is fresh; you're dead.
%
    if armflag == 1  
       if jac_age > 0
          sol = xold;
          x = xold; f0 = fold;    
          disp('Armijo failure; recompute Jacobian.');
       else
          disp('Complete Armijo failure.');
          sol = xold;
          ierr = 2;
          return
       end
    end
    fnrm = norm(f0);
    it_hist = [it_hist',[fnrm,iarm]']';
    if nargout == 4, x_hist = [x_hist,x]; end
    rat = fnrm/fnrmo;
    if debug == 1, disp([itc fnrm rat]); end
    outstat(itc+1, :) = [itc fnrm rat];
% end while
end
sol = x;
if debug == 1, disp(outstat); end
%
% on failure, set the error flag
%
if fnrm > stop_tol, ierr = 1; end
%
end

function [sol, it_hist, ierr, x_hist] = nsoli(x,f,tol, parms,pars_in)
% NSOLI  Newton-Krylov solver, globally convergent 
%        solver for f(x) = 0
%
% Inexact-Newton-Armijo iteration
%
% Eisenstat-Walker forcing term
%
% Parabolic line search via three point interpolation.
%
% C. T. Kelley, April 27, 2001
%
% This code comes with no guarantee or warranty of any kind.
%
% function [sol, it_hist, ierr, x_hist] = nsoli(x,f,tol,parms)
%
% inputs:
%        initial iterate = x
%        function = f
%        tol = [atol, rtol] relative/absolute
%            error tolerances for the nonlinear iteration
%        parms = [maxit, maxitl, etamax, lmeth, restart_limit]
%            maxit = maxmium number of nonlinear iterations
%                default = 40
%            maxitl = maximum number of inner iterations before restart
%                in GMRES(m), m = maxitl 
%                default = 40
%                
%                For iterative methods other than GMRES(m) maxitl
%                is the upper bound on linear iterations.
%
%            |etamax| = Maximum error tolerance for residual in inner
%                iteration. The inner iteration terminates
%                when the relative linear residual is
%                smaller than eta*| F(x_c) |. eta is determined
%                by the modified Eisenstat-Walker formula if etamax > 0.
%                If etamax < 0, then eta = |etamax| for the entire
%                iteration.
%                default: etamax = .9
%
%            lmeth = choice of linear iterative method
%                    1 (GMRES), 2 GMRES(m), 
%                    3 (BICGSTAB), 4 (TFQMR)
%                 default = 1 (GMRES, no restarts)
%
%            restart_limit = max number of restarts for GMRES if
%                    lmeth = 2
%                  default = 20
%
% output:
%        sol = solution
%        it_hist(maxit,3) = l2 norms of nonlinear residuals
%            for the iteration, number of function evaluations,
%            and number of steplength reductions
%        ierr = 0 upon successful termination
%        ierr = 1 if after maxit iterations
%             the termination criterion is not satsified
%        ierr = 2 failure in the line search. The iteration
%             is terminated if too many steplength reductions
%             are taken.
%
%    x_hist = matrix of the entire interation history.
%             The columns are the nonlinear iterates. This
%             is useful for making movies, for example, but
%             can consume way too much storage. This is an
%             OPTIONAL argument. Storage is only allocated
%             if x_hist is in the output argument list.
%
%
%
% internal parameters:
%       debug = turns on/off iteration statistics display as
%               the iteration progresses
%
%       alpha = 1.d-4, parameter to measure sufficient decrease
%
%       sigma0 = .1, sigma1 = .5, safeguarding bounds for the linesearch
%
%       maxarm = 20, maximum number of steplength reductions before
%                    failure is reported
%
%
%
%
% Set the debug parameter; 1 turns display on, otherwise off.
%
debug = 0;
%
% Set internal parameters.
%
alpha = 1.d-4; sigma0 = .1; sigma1 = .5; maxarm = 20; gamma = .9;
%
% Initialize it_hist, ierr, x_hist, and set the default values of
% those iteration parameters which are optional inputs.
%
ierr = 0; maxit = 40; lmaxit = 40; etamax = .9; it_histx = zeros(maxit,3);
lmeth = 1; restart_limit = 20;
if nargout == 4, x_hist = x; end
%
% Initialize parameters for the iterative methods.
% Check for optional inputs.
%
gmparms = [abs(etamax), lmaxit];
if nargin == 5 % jdw: fudge!
    maxit = parms(1); lmaxit = parms(2); etamax = parms(3);
    it_histx = zeros(maxit,3);
    gmparms = [abs(etamax), lmaxit];
    if length(parms)>= 4
       lmeth = parms(4);
    end
    if length(parms) == 5
       gmparms = [abs(etamax), lmaxit, parms(5), 1];
    end
end
%
rtol = tol(2); atol = tol(1); n = length(x); fnrm = 1; itc = 0;
%
% Evaluate f at the initial iterate,and
% compute the stop tolerance.
%
%f0 = feval(f,x);
[f0,~,parameters] = calc_F(x,pars_in);

fnrm = norm(f0);
it_histx(itc+1,1) = fnrm; it_histx(itc+1,2) = 0; it_histx(itc+1,3) = 0;
fnrmo = 1;
stop_tol = atol + rtol*fnrm;
outstat(itc+1, :) = [itc fnrm 0 0 0];
%
% main iteration loop
%
while(fnrm > stop_tol & itc < maxit)
%
% Keep track of the ratio (rat = fnrm/frnmo)
% of successive residual norms and 
% the iteration counter (itc).
%
    % JDW: update preconditioner in outer iteration only?
%     if itc>0 & parameters.gen_pars.update_precond
%         disp('updating preconditioner...')
%         parameters.gen_pars.left_precond=false;
%         pars_in{5}=parameters;
%         [tmp_x,~,parameters]=calc_F(x,pars_in);
%         [parameters] = calc_preconditioner(tmp_x,pars_in);
%         parameters.gen_pars.left_precond=true;
%         pars_in{5}=parameters;
%     end
    
    rat = fnrm/fnrmo;
    fnrmo = fnrm; 
    itc = itc+1;
    [step, errstep, inner_it_count,inner_f_evals] = ...
         dkrylov(f0, f, x, gmparms, lmeth, pars_in);

%
%   The line search starts here.
%
    xold = x;
    lambda = 1; lamm = 1; lamc = lambda; iarm = 0;
    xt = x + lambda*step;
    %ft = feval(f,xt); jdw
    ft = calc_F(xt,pars_in);
    nft = norm(ft); nf0 = norm(f0); ff0 = nf0*nf0; ffc = nft*nft; ffm = nft*nft;
    while nft >= (1 - alpha*lambda) * nf0;
%
%   Apply the three point parabolic model.
%
        if iarm == 0
            lambda = sigma1*lambda; 
        else
            lambda = parab3p(lamc, lamm, ff0, ffc, ffm); 
        end
%
% Update x; keep the books on lambda.
%
        xt = x+lambda*step;
        lamm = lamc;
        lamc = lambda;
%
% Keep the books on the function norms.
%
        %ft = feval(f,xt); jdw
        ft = calc_F(xt,pars_in);
        nft = norm(ft);
        ffm = ffc;
        ffc = nft*nft;
        iarm = iarm+1;
        if iarm > maxarm
            disp(' Armijo failure, too many reductions ');
            ierr = 2;
            disp(outstat)
            it_hist = it_histx(1:itc+1,:);
        if nargout == 4, x_hist = [x_hist,x]; end
            sol = xold;
            return;
        end
    end
    x = xt;
    f0 = ft;
%
%   End of line search.
%
    if nargout == 4, x_hist = [x_hist,x]; end
    fnrm = norm(f0);
    it_histx(itc+1,1) = fnrm; 
%
%   How many function evaluations did this iteration require?
%
    it_histx(itc+1,2) = it_histx(itc,2)+inner_f_evals+iarm+1;
    if itc == 1, it_histx(itc+1,2) = it_histx(itc+1,2)+1; end;
    it_histx(itc+1,3) = iarm;
%
    rat = fnrm/fnrmo;
%
%   Adjust eta as per Eisenstat-Walker.
%
    if etamax > 0
        etaold = gmparms(1);
        etanew = gamma*rat*rat;
        if gamma*etaold*etaold > .1
            etanew = max(etanew,gamma*etaold*etaold);
        end
        gmparms(1) = min([etanew,etamax]);
        gmparms(1) = max(gmparms(1),.5*stop_tol/fnrm);
    end
%
    outstat(itc+1, :) = [itc fnrm inner_it_count rat iarm];
%
end
sol = x;
it_hist = it_histx(1:itc+1,:);
if debug == 1
    disp(outstat)
    it_hist = it_histx(1:itc+1,:);
end
%
% on failure, set the error flag
%
if fnrm > stop_tol, ierr = 1; end
end

%%
function [step,iarm,xp,fp,armflag] = armijo(direction,x,f0,f,maxarm,pars_in)
iarm = 0;
sigma1 = .5;
alpha = 1.d-4;
armflag = 0;
xp = x; fp = f0; 
%
    xold = x;
    lambda = 1; lamm = 1; lamc = lambda; iarm = 0;
    step = lambda*direction;
    xt = x + step;
    %ft = feval(f,xt); jdw
    ft = calc_F(xt,pars_in);
    nft = norm(ft); nf0 = norm(f0); ff0 = nf0*nf0; ffc = nft*nft; ffm = nft*nft;
    while nft >= (1 - alpha*lambda) * nf0;
%
%   Apply the three point parabolic model.
%
        if iarm == 0
            lambda = sigma1*lambda;
        else
            lambda = parab3p(lamc, lamm, ff0, ffc, ffm);
        end
%
% Update x; keep the books on lambda.
%
        step = lambda*direction;
        xt = x + step;
        lamm = lamc;
        lamc = lambda;
%
% Keep the books on the function norms.
%
        %ft = feval(f,xt); jdw
        ft = calc_F(xt,pars_in);
        nft = norm(ft);
        ffm = ffc;
        ffc = nft*nft;
        iarm = iarm+1;
        if iarm > maxarm
            disp(' Armijo failure, too many reductions ');
            armflag = 1;
            sol = xold;
            return;
        end
    end
    xp = xt; fp = ft;
%
%   end of line search
end

function lambdap = parab3p(lambdac, lambdam, ff0, ffc, ffm)
% Apply three-point safeguarded parabolic model for a line search.
%
% C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
% function lambdap = parab3p(lambdac, lambdam, ff0, ffc, ffm)
%
% input:
%       lambdac = current steplength
%       lambdam = previous steplength
%       ff0 = value of \| F(x_c) \|^2
%       ffc = value of \| F(x_c + \lambdac d) \|^2
%       ffm = value of \| F(x_c + \lambdam d) \|^2
%
% output:
%       lambdap = new value of lambda given parabolic model
%
% internal parameters:
%       sigma0 = .1, sigma1 = .5, safeguarding bounds for the linesearch
%

%
% Set internal parameters.
%
sigma0 = .1; sigma1 = .5;
%
% Compute coefficients of interpolation polynomial.
%
% p(lambda) = ff0 + (c1 lambda + c2 lambda^2)/d1
%
% d1 = (lambdac - lambdam)*lambdac*lambdam < 0
%      so, if c2 > 0 we have negative curvature and default to
%      lambdap = sigam1 * lambda.
%
c2 = lambdam*(ffc-ff0)-lambdac*(ffm-ff0);
if c2 >= 0
    lambdap = sigma1*lambdac; return
end
c1 = lambdac*lambdac*(ffm-ff0)-lambdam*lambdam*(ffc-ff0);
lambdap = -c1*.5/c2;
if lambdap < sigma0*lambdac, lambdap = sigma0*lambdac; end
if lambdap > sigma1*lambdac, lambdap = sigma1*lambdac; end

end

%%
function [step, errstep, total_iters, f_evals] = dkrylov(f0, f, x, params, lmeth , pars_in)
% Krylov linear equation solver for use in nsoli
%
% C. T. Kelley, April 1, 2003
%
%
% This code comes with no guarantee or warranty of any kind.
%
% function [step, errstep, total_iters, f_evals] 
%                              = dkrylov(f0, f, x, params, lmeth)
%
%
% Input:  f0 = function at current point
%         f = nonlinear function
%              the format for f is  function fx = f(x)
%              Note that for Newton-GMRES we incorporate any
%              preconditioning into the function routine.
%         x = current point
%         params = vector to control iteration
%              params(1) = relative residual reduction factor
%              params(2) = max number of iterations
%              params(3) = max number of restarts for GMRES(m)
%              params(4) (Optional) = reorthogonalization method in GMRES
%                   1 -- Brown/Hindmarsh condition (default)
%                   2 -- Never reorthogonalize (not recommended)
%                   3 -- Always reorthogonalize (not cheap!)
%
%         lmeth = method choice
%              1 GMRES without restarts (default)
%              2 GMRES(m), m = params(2) and the maximum number
%                   of restarts is params(3) 
%              3 Bi-CGSTAB
%              4 TFQMR
%
% Output: x = solution
%         errstep = vector of residual norms for the history of
%                 the iteration
%         total_iters = number of iterations
%
%

%
% initialization
%
lmaxit = params(2);
restart_limit = 20;
if length(params) >= 3
    restart_limit = params(3);
end
if lmeth == 1, restart_limit = 0; end
if length(params) == 3
%
% default reorthogonalization
%
     gmparms = [params(1), params(2), 1];
elseif length(params) == 4
%
% reorthogonalization method is params(4)
%
     gmparms = [params(1), params(2), params(4)];
else
     gmparms = [params(1), params(2)];
end
%
% linear iterative methods
%
if lmeth == 1 | lmeth == 2  % GMRES or GMRES(m) 
%
% compute the step using a GMRES routine especially designed
% for this purpose
%
    [step, errstep, total_iters] = dgmres(f0, f, x, gmparms,pars_in);
    kinn = 0;
%
%   restart at most restart_limit times
%
    while total_iters == lmaxit & ...
          errstep(total_iters) > gmparms(1)*norm(f0) & ...
          kinn < restart_limit
        kinn = kinn+1;
        [step, errstep, total_iters] = dgmres(f0, f, x, gmparms,step,pars_in);
    end
    total_iters = total_iters+kinn*lmaxit;
    f_evals = total_iters+kinn;
%
% Bi-CGSTAB
%
elseif lmeth == 3
    [step, errstep, total_iters] = dcgstab(f0, f, x, gmparms,pars_in);
    f_evals = 2*total_iters;
%
% TFQMR
%
elseif lmeth == 4 
    [step, errstep, total_iters] = dtfqmr(f0, f, x, gmparms,pars_in);
    f_evals = 2*total_iters;
else
    error(' lmeth error in fdkrylov')
end

end

%%
function z = dirder(x,w,f,f0,pars_in)
% Finite difference directional derivative
% Approximate f'(x) w
% 
% C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
% function z = dirder(x,w,f,f0)
%
% inputs:
%           x, w = point and direction
%           f = function
%           f0 = f(x), in nonlinear iterations
%                f(x) has usually been computed
%                before the call to dirder

%
% Use a hardwired difference increment.
%
epsnew = 1.d-7;
%
n = length(x);
%
% scale the step
%
if norm(w) == 0
    z = zeros(n,1);
return
end
%
% Now scale the difference increment.
%
xs=(x'*w)/norm(w);
if xs ~= 0.d0
     epsnew=epsnew*max(abs(xs),1.d0)*sign(xs);
end
epsnew=epsnew/norm(w);
%
% del and f1 could share the same space if storage
% is more important than clarity.
%
del = x+epsnew*w;
%f1 = feval(f,del); jdw
f1 = calc_F(del,pars_in);
z = (f1 - f0)/epsnew;

end

%%
function [x, error, total_iters] = dgmres(f0, f, xc, params, xinit,pars_in)
% GMRES linear equation solver for use in Newton-GMRES solver
%
% C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
% function [x, error, total_iters] = dgmres(f0, f, xc, params, xinit)
%
%
% Input:  f0 = function at current point
%         f = nonlinear function
%              the format for f is  function fx = f(x)
%              Note that for Newton-GMRES we incorporate any 
%              preconditioning into the function routine.
%         xc = current point
%         params = two dimensional vector to control iteration
%              params(1) = relative residual reduction factor
%              params(2) = max number of iterations
%            params(3) (Optional) = reorthogonalization method
%                   1 -- Brown/Hindmarsh condition (default)
%                   2 -- Never reorthogonalize (not recommended)
%                   3 -- Always reorthogonalize (not cheap!)
%
%         xinit = initial iterate. xinit = 0 is the default. This
%              is a reasonable choice unless restarted GMRES
%              will be used as the linear solver.
%
% Output: x = solution
%         error = vector of residual norms for the history of
%                 the iteration
%         total_iters = number of iterations
%
% Requires givapp.m, dirder.m 

%
% initialization
%
errtol = params(1);
kmax = params(2);
reorth = 1;
if length(params) == 3
    reorth = params(3);
end
%
% The right side of the linear equation for the step is -f0. 
%
b = -f0;
n = length(b);
%
% Use zero vector as initial iterate for Newton step unless
% the calling routine has a better idea (useful for GMRES(m)).
%
x = zeros(n,1);
r = b;
if nargin == 6 % jdw ==5 to ==6 to account for extra pars_in argument
    x = xinit;
    r = -dirder(xc, x, f, f0,pars_in)-f0;
else
    pars_in = xinit; % jdw if <6 then pars_in will be the last argument, enetered under xinit
end
%
%
h = zeros(kmax);
v = zeros(n,kmax);
c = zeros(kmax+1,1);
s = zeros(kmax+1,1);
rho = norm(r);
g = rho*eye(kmax+1,1);
errtol = errtol*norm(b);
error = [];
%
% Test for termination on entry.
%
error = [error,rho];
total_iters = 0;
if(rho < errtol) 
%   disp(' early termination ')
return
end
%
%
v(:,1) = r/rho;
beta = rho;
k = 0;
%
% GMRES iteration
%
while((rho > errtol) & (k < kmax))
    k = k+1;
%
%   Call directional derivative function.
%
    v(:,k+1) = dirder(xc, v(:,k), f, f0,pars_in);
    normav = norm(v(:,k+1));
%
%   Modified Gram-Schmidt
%
    for j = 1:k
        h(j,k) = v(:,j)'*v(:,k+1);
        v(:,k+1) = v(:,k+1)-h(j,k)*v(:,j);
    end
    h(k+1,k) = norm(v(:,k+1));
    normav2 = h(k+1,k);
%
%   Reorthogonalize?
%
if  (reorth == 1 & normav + .001*normav2 == normav) | reorth ==  3
    for j = 1:k
        hr = v(:,j)'*v(:,k+1);
	h(j,k) = h(j,k)+hr;
        v(:,k+1) = v(:,k+1)-hr*v(:,j);
    end
    h(k+1,k) = norm(v(:,k+1));
end
%
%   Watch out for happy breakdown.
%
    if(h(k+1,k) ~= 0)
    v(:,k+1) = v(:,k+1)/h(k+1,k);
    end
%
%   Form and store the information for the new Givens rotation.
%
    if k > 1
        h(1:k,k) = givapp(c(1:k-1),s(1:k-1),h(1:k,k),k-1);
    end
%
%   Don't divide by zero if solution has  been found.
%
    nu = norm(h(k:k+1,k));
    if nu ~= 0
%        c(k) = h(k,k)/nu;
        c(k) = conj(h(k,k)/nu);
        s(k) = -h(k+1,k)/nu;
        h(k,k) = c(k)*h(k,k)-s(k)*h(k+1,k);
        h(k+1,k) = 0;
        g(k:k+1) = givapp(c(k),s(k),g(k:k+1),1);
    end
%
%   Update the residual norm.
%
    rho = abs(g(k+1));
    error = [error,rho];
%
%   end of the main while loop
%
end
%
% At this point either k > kmax or rho < errtol.
% It's time to compute x and leave.
%
y = h(1:k,1:k)\g(1:k);
total_iters = k;
x = x + v(1:n,1:k)*y;

end


%%
function vrot = givapp(c,s,vin,k)
%  Apply a sequence of k Givens rotations, used within gmres codes.
% 
%  C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
%  function vrot = givapp(c, s, vin, k)
%
vrot = vin;
for i = 1:k
    w1 = c(i)*vrot(i)-s(i)*vrot(i+1);
%
%   Here's a modest change that makes the code work in complex
%   arithmetic. Thanks to Howard Elman for this.
%
%    w2 = s(i)*vrot(i)+c(i)*vrot(i+1);
    w2 = s(i)*vrot(i)+conj(c(i))*vrot(i+1);
    vrot(i:i+1) = [w1,w2];
end

end
%%
function [x, error, total_iters] = ...
                     dcgstab(f0, f, xc, params, xinit,pars_in)
% Forward difference Bi-CGSTAB solver for use in nsoli
%
% C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
% function [x, error, total_iters]
%                    = dcgstab(f0, f, xc, params, xinit)
%
% Input:  f0 = function at current point
%         f = nonlinear function
%              the format for f is  function fx = f(x)
%              Note that for Newton-GMRES we incorporate any
%              preconditioning into the function routine.
%         xc = current point
%         params = two dimensional vector to control iteration
%              params(1) = relative residual reduction factor
%              params(2) = max number of iterations
%
%         xinit = initial iterate. xinit = 0 is the default. This
%              is a reasonable choice unless restarts are needed.
%
%
% Output: x = solution
%         error = vector of residual norms for the history of
%                 the iteration
%         total_iters = number of iterations
%
% Requires: dirder.m
%

%
% initialization
%
b = -f0; 
n = length(b); errtol = params(1)*norm(b); kmax = params(2); error = []; 
rho = zeros(kmax+1,1);
%
% Use zero vector as initial iterate for Newton step unless
% the calling routine has a better idea (useful for GMRES(m)).
%
x = zeros(n,1);
r = b;
if nargin == 6 % jdw ==5 to ==6
    x = xinit;
    r = -dirder(xc, x, f, f0,pars_in)-f0;
else
    pars_in=xinit;
end
%
hatr0 = r;
k = 0; rho(1) = 1; alpha = 1; omega = 1;
v = zeros(n,1); p = zeros(n,1); rho(2) = hatr0'*r;
zeta = norm(r); error = [error,zeta];
%
% Bi-CGSTAB iteration
%
while((zeta > errtol) & (k < kmax))
    k = k+1;
    if omega == 0
       error('Bi-CGSTAB breakdown, omega = 0');
    end
    beta = (rho(k+1)/rho(k))*(alpha/omega);
    p = r+beta*(p - omega*v);
    v = dirder(xc,p,f,f0,pars_in);
    tau = hatr0'*v;
    if tau == 0
        error('Bi-CGSTAB breakdown, tau = 0');
    end 
    alpha = rho(k+1)/tau;
    s = r-alpha*v; 
    t = dirder(xc,s,f,f0,pars_in);
    tau = t'*t;
    if tau == 0
       error('Bi-CGSTAB breakdown, t = 0');
    end
    omega = t'*s/tau; 
    rho(k+2) = -omega*(hatr0'*t);
    x = x+alpha*p+omega*s;
    r = s-omega*t;
    zeta = norm(r);
    total_iters = k;
    error = [error, zeta];
end

end

%%
function [x, error, total_iters] = ...
                     dtfqmr(f0, f, xc, params, xinit,pars_in)
% Forward difference TFQMR solver for use in nsoli
%
% C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
% function [x, error, total_iters]
%                    = dtfqmr(f0, f, xc, params, xinit)
%
%
%
% Input:  f0 = function at current point
%         f = nonlinear function
%              the format for f is  function fx = f(x)
%              Note that for Newton-GMRES we incorporate any
%              preconditioning into the function routine.
%         xc = current point
%         params = two dimensional vector to control iteration
%              params(1) = relative residual reduction factor
%              params(2) = max number of iterations
%
%         xinit = initial iterate. xinit = 0 is the default. This
%              is a reasonable choice unless restarts are needed.
%
%
% Output: x = solution
%         error = vector of residual norms for the history of
%                 the iteration
%         total_iters = number of iterations
%
% Requires: dirder.m
%

%
% initialization
%
b = -f0;
n = length(b); errtol = params(1)*norm(b); kmax = params(2); error = []; 
x = zeros(n,1);
r = b;
if nargin == 6 % jdw ==5 to ==6
    x = xinit;
    r = -dirder(xc, x, f, f0,pars_in)-f0;
else
    pars_in=xinit;
end
%
u = zeros(n,2); y = zeros(n,2); w = r; y(:,1) = r; 
k = 0; d = zeros(n,1); 
v = dirder(xc, y(:,1),f,f0,pars_in);
u(:,1) = v;
theta = 0; eta = 0; tau = norm(r); error = [error,tau];
rho = tau*tau;
%
% TFQMR iteration
%
while( k < kmax)
    k = k+1;
    sigma = r'*v;
%
    if sigma == 0
        error('TFQMR breakdown, sigma = 0')
    end
%
    alpha = rho/sigma;
%
% 
%
    for j = 1:2
%
%   Compute y2 and u2 only if you have to
%
        if j == 2 
            y(:,2) = y(:,1)-alpha*v;
            u(:,2) = dirder(xc, y(:,2),f,f0,pars_in);
        end
        m = 2*k-2+j;
        w = w-alpha*u(:,j);
        d = y(:,j)+(theta*theta*eta/alpha)*d;
        theta = norm(w)/tau; c = 1/sqrt(1+theta*theta);
        tau = tau*theta*c; eta = c*c*alpha;
        x = x+eta*d;
%
%   Try to terminate the iteration at each pass through the loop
%
        if tau*sqrt(m+1) <=  errtol
            error = [error, tau];
            total_iters = k;
            return
        end
    end
%
%
%
    if rho == 0
        error('TFQMR breakdown, rho = 0')
    end
%
    rhon = r'*w; beta = rhon/rho; rho = rhon;
    y(:,1) = w + beta*y(:,2);
    u(:,1) = dirder(xc, y(:,1),f,f0,pars_in);
    v = u(:,1)+beta*(u(:,2)+beta*v);
    error = [error, tau];
    total_iters = k;
end

end
%
%

