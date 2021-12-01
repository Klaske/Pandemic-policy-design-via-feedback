classdef Anderson_COVID_SEEIR
    %
    % Matlab implementation of SEEIR model Anderson et al. for BC
    %
    
   properties
      prms              % Model parameters
      settings          % Simulation settings
      cntr              % Control settings
      Res               % Placeholder for simulation results
   end
   
   methods
      function obj = Anderson_COVID_SEEIR(varargin)
          %
          % Anderson_COVID_SEEIR(NonDefault)
          % 
          % Initialize simulation class obj = Anderson_COVID_SEEIR(NonDefault)
          % 
          % Initialized with default values according to "
          % Estimating the impact of COVID-19 control measures using a 
          % Bayesian model of physical distancing", S. C. Anderson et al. 
          % medrxiv April 17, 2020
          %
          % Any settings, parameters of controller settings other than default 
          % can be defined in NonDefault (optional). (Similar to Update_Settings, see below)
          % 
          % Inputs:
          % NonDefault - optional - structure of the form update.{prms, settings, cntr}  
          %
          N = 5100000;  % Population
          i0 = 8;       % Seeding pandemic
          % SEEIR parameters
          params.D    = 5; 
          params.k1   = 0.2;
          params.k2   = 1;
          params.N    = N;
          params.q    = 0.05;
          params.ud    = 0.1;
          params.R0     = 3; % This is a bit higher than what Anderson et al. found for the average
          params.ur     = 0.02;
          % Parameters for variant (no variant by default) 
          params.R0v    = 0; 
          params.Iv_t   = [];
          
          % Simulation settings (for ODE solver) 
          settings.maxstep = 0.5;
          settings.Tspan = [0:365]; % time in days
          
          % Social distancing parameters according to Anderson et al. 
          cntr.f        = 1-0.78; % 
          % t0 = Feb 1st 2020
          % Gradual ramp up of f with 
          % t1 = March 15 and t2 = March 22. 
          cntr.ton      = 29+15; 
          cntr.ton2     = 29+15+7;
          cntr.toff     = 29+15+7+10*30; % 10 months after t2
             
          % Parameters for observed case model (not used) 
          cntr.ws       = (wblpdf([0:1:45], 9.85, 1.73));
          cntr.k2       = params.k2;
          cntr.psi      = [0.1*ones(1, 29+14+1), 0.3*ones(1,settings.Tspan(end)-(29+14))]; % +1 as the simulation starts at zero
          
          % Feedback controller settings (note that the controller
          % parameters are not initialized and need to be defined by an
          % update prior to simulation          
          cntr.h        = 1;
          cntr.t_ini_fix= 7;
          cntr.Tt       = 20; 
          cntr.I0       = .5
          
          cntr.setting  = 'PID';
          cntr.PIDsetting = 'NonLinPID';  % Used in delayed implementation
          cntr.sp       = 100;
          
          % SEEIR model does not include explicit delays. The delay is
          % implemented in the (discrete) controller, and not as part of
          % the ODE model.
          cntr.delayI   = 14;
          cntr.delayE   = 6;
          % No delay unless specified. No vaccination unless specified.
          % Initialize at zero
          cntr.d        = zeros(settings.Tspan(end)+1, 1); 
          cntr.p_vac    = zeros(settings.Tspan(end)+1, 1); 
          
          % Initial conditions for pandemic simulation (x0)
          fsi = params.ud/(params.ud+params.ur);
          nsi = 1-fsi;
          settings.x0 = [nsi*(N-i0), 0.4*nsi*i0, 0.1*nsi*i0, 0.5*nsi*i0, 0, ...
              fsi*(N-i0), 0.4*fsi*i0, 0.1*fsi*i0, 0.5*fsi*i0, 0 , zeros(1, 8)];

          % Save the parameters in the class
          obj.prms = params; obj.cntr = cntr; obj.settings = settings;
          
          % Allow for changed settings with input. This only replaces the
          % fieldnames present in the NonDefault structure
          if ~isempty(varargin)
              fn = fieldnames(varargin{1});
              for n = fn'
                  prmn_fn = fieldnames(varargin{1}.(n{1}));
                  for p_n = prmn_fn'
                    obj.(n{1}).(p_n{1}) = varargin{1}.(n{1}).(p_n{1});
                  end
              end
          end
      end
      
      function obj = Update_Settings(obj, update)
          %
          % obj = Update_Settings(obj, update)
          % 
          % Update settings, this function updates only the fields in the
          % update structure. 
          %
          % Input: 
          % update  - structure of the form update.{prms, settings, cntr}
          %
          fn = fieldnames(update);
          for n = fn'
              prmn_fn = fieldnames(update.(n{1}));
              for p_n = prmn_fn'
                obj.(n{1}).(p_n{1}) = update.(n{1}).(p_n{1});
              end
          end
          if obj.settings.Tspan(end)> length(obj.cntr.d)
              obj.cntr.d = [obj.cntr.d; zeros(obj.settings.Tspan(end)-length(obj.cntr.d)+1, 1)];
              obj.cntr.p_vac = [obj.cntr.p_vac; zeros(obj.settings.Tspan(end)-length(obj.cntr.p_vac)+1, 1)];
          end
          obj.cntr.k2 = obj.prms.k2; % Ensuring values used in simulations are the same.
      end
      
      function [obj, T, X, f] = Simulate(obj, model, varargin)
          % 
          % obj = Simulate(obj, model, {disturbance} )
          % 
          % Simulates COVID model according to settings, prms and cntr. 
          % 
          % Inputs:
          % model:  str  - Results saved in 
          %                 obj.Res.(model)
          % disturbance (optional) : This input is used as the disturbance signal, and cntr.d as 
          %                          saved in the class is replaced by this
          %                          input. Note that this affects future
          %                          simulations too. 
          %
          
          if nargin > 2
              % update the disturbance cntr.d
              obj.cntr.d(varargin{1})=varargin{2};
          end
          
          % ODE solver settings 
          opts = odeset('MaxStep',obj.settings.maxstep, 'InitialStep', obj.settings.maxstep);
          
          % Controller code runs outside of ODE solver, (discrete) and is
          % called by ODE solver. The memory needs to be cleared (reset)
          % prior to simulation. 
          clear COVID_Control_PID
          % Run simulation
          [T, X] = ode45(@(t, x) Anderson_COVID(t, x, obj.prms, obj.cntr),obj.settings.Tspan, obj.settings.x0, opts);
          % Save results
          obj.Res.(model).T = T;obj.Res.(model).X = X;
          
          % Generate controller output, disturbance and reported cases
          clear COVID_Control_PID
          f = zeros(size(obj.Res.(model).T)); d = zeros(size(obj.Res.(model).T)); y_rep = zeros(size(obj.Res.(model).T));p_vac = zeros(size(obj.Res.(model).T));
          for idx = 1:length(obj.Res.(model).T)
              [f(idx), d(idx), y_rep(idx), p_vac(idx)] = COVID_Control_PID(obj.Res.(model).T(idx), obj.Res.(model).X(idx, :), obj.cntr);
          end
          
          % Save results
          obj.Res.(model).f             = f;
          obj.Res.(model).y_reported    = y_rep;
          obj.Res.(model).d             = d;
          obj.Res.(model).p_vac         = p_vac;
      end
      
      function obj = SMCController(obj, k, tau2, delay, tau_c, varargin)
          %
          % obj = SMCController(obj, k, tau2, delay, tau_c)
          % 
          % Tune controller using SMC rules (Skogestad) for integral system
          % with lag and delay
          %
          % Inputs: 
          % k       - system gain
          % tau2    - system lag
          % delay   - system delay
          % tau_c   - desired closed-loop time constant
          % PI(optional) - if PI (varargin) is not empty, PI tuning rules
          %                are used rather than PID, i.e. kd = 0; 
          %   
          % Output:
          % obj     - Calculated controller parameters are saved in
          %         obj.cntr
          %
          if nargin > 5 && ~isempty(varargin{1})
              % Make it a PI controller instead
              delay = delay +tau2; 
              tau2 = 0;
          end
          Kc = 1/k*1/(tau_c + delay); tauI = 4*(tau_c +delay); tauD = tau2;
          obj.cntr.kp = Kc*(1+ tauD/tauI);
          Ti = tauI*(1+ tauD/tauI); Td = tauD/ (1+ tauD/tauI) ;
          obj.cntr.ki = obj.cntr.kp/Ti;
          obj.cntr.kd = obj.cntr.kp*Td;
      end
     
   end
end


function dxdt = Anderson_COVID(t, x, params, cntr)
    %
    % SEEIQR model with social distancing
    %
    S = x(1); E1 = x(2); E2 = x(3); I = x(4); Q = x(5);

    Sd = x(6); E1d = x(7); E2d = x(8); Id = x(9); Qd = x(10);
    
    E1v = x(11); E2v = x(12); Iv = x(13); Qv = x(14);

    E1dv = x(15); E2dv = x(16); Idv = x(17); Qdv = x(18);

    R0 = params.R0;D = params.D;k1 = params.k1;k2 = params.k2;
    N = params.N;q = params.q;ur = params.ur; r = params.ud;
    R0v = params.R0v;
     
    [f, d, ~, p_vac] = COVID_Control_PID(t, x, cntr);
    if isfield(params, 'u_offset')
        f = f + params.u_offset(floor(t)+1); 
    end

    % Initialization (seed) of variant
    if floor(t) == params.Iv_t
        Iv_in = params.Iv_in; 
    else 
        Iv_in = 0; 
    end
    % Both variants come out of susceptible class (i.e. assumes no re-infection)
    dSdt = -(R0/(D+1/k2))*(I+E2 + f*(Id+E2d))*S/N - (R0v/(D+1/k2))*(Iv+E2v + f*(Idv+E2dv))*S/N - r*S + ur*Sd -p_vac*N*S/(S+Sd); % Reduce susceptible by vaccination, d_vac given as daily percentage (low numbers)
    dE1dt = (R0/(D+1/k2))*(I+E2 + f*(Id+E2d))*S/N - k1*E1 -r*E1 + ur*E1d + d; % add disturbance, new cases from external source 
    dE2dt = k1*E1 -k2*E2 -r*E2 + ur*E2d;
    dIdt = k2*E2 - q*I -  I/D - r*I+ ur*Id;
    dQdt = q*I - Q/D -r*Q + ur*Qd;
    %dRdt = I/D + Q/D -r*R+ur*Rd;
    % Add variant with higher R0
    dE1vdt = (R0v/(D+1/k2))*(Iv+E2v + f*(Idv+E2dv))*S/N - k1*E1v -r*E1v + ur*E1dv +Iv_in; 
    dE2vdt = k1*E1v -k2*E2v -r*E2v + ur*E2dv;
    dIvdt = k2*E2v - q*Iv -  Iv/D - r*Iv+ ur*Idv; 
    dQvdt = q*Iv - Qv/D -r*Qv + ur*Qdv;
    
    dSddt = -(f*R0/(D+1/k2))*(I+E2 + f*(Id+E2d))*Sd/N - (f*R0v/(D+1/k2))*(Iv+E2v + f*(Idv+E2dv))*Sd/N+ r*S -ur*Sd -p_vac*N*Sd/(S+Sd); % Reduce susceptible by vaccination, d given as daily percentage (low numbers)
    dE1ddt = (f*R0/(D+1/k2))*(I+E2 + f*(Id+E2d))*Sd/N - k1*E1d + r*E1 - ur*E1d;
    dE2ddt = k1*E1d - k2*E2d + r*E2 - ur*E2d;
    dIddt = k2*E2d - q*Id-  Id/D + r*I - ur*Id;
    dQddt = q*Id - Qd/D +r*Q - ur*Qd;
    %dRddt = Id/D+Qd/D +r*R - ur*Rd;
    dE1dvdt = (f*R0v/(D+1/k2))*(Iv+E2v + f*(Idv+E2dv))*Sd/N - k1*E1dv + r*E1v - ur*E1dv;
    dE2dvdt = k1*E1dv - k2*E2dv + r*E2v - ur*E2dv;
    dIdvdt = k2*E2dv - q*Idv-  Idv/D + r*Iv - ur*Idv;
    dQdvdt = q*Idv - Qdv/D +r*Qv - ur*Qdv;
    dxdt = [dSdt dE1dt dE2dt dIdt dQdt dSddt dE1ddt dE2ddt dIddt dQddt dE1vdt dE2vdt dIvdt dQvdt dE1dvdt dE2dvdt dIdvdt dQdvdt]'; % Excluded R
end


