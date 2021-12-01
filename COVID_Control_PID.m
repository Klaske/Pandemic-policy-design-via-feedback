function [f, d, y_rep, p_vac] = COVID_Control_PID(t, x, cntr)
    %
    % [f, d, y_rep, p_vac] = COVID_Control_PID(t, x, cntr)
    % 
    %   Calculate control action according to controller details in cntr
    % 
    %   Inputs: 
    %   t, x    - time and state from ODE definition
    %   cntr    - structure with controller settings and parameters
    % 
    %   Outputs: 
    %   f       - control action
    %   d       - disturbance on E1, import of new cases
    %   y_rep   - y used for feedback
    %   p_vac   - daily percentage of vaccinations 
    
    persistent per_var;
    % initialize
    if any([t == 0, isempty(per_var)])
        per_var.f_trigger = 'off';
        per_var.y_memory =  zeros(max(max(cntr.delayI, cntr.delayE), length(cntr.ws))+2, 1);
        per_var.y_rep_memory =  zeros(length(cntr.ws)+2, 1);
        per_var.E2_memory =  sum(x([3, 8]))*ones(length(cntr.ws)+2, 1);
        per_var.t_trigger = 0;
        per_var.f_kmin1 = 0;
        per_var.I = cntr.I0;
        per_var.f_delayed = cntr.f(1); 
        per_var.t_update_last = 0;
        per_var.dE1 = true;
        if isfield(cntr, 'filter')
            per_var.filt_u = zeros(size(cntr.filter.num{:},2)-1,1);
            per_var.filt_y = zeros(size(cntr.filter.den{:},2)-1,1);
        end
    end
    
    % Choose policy according to cntr.PIDsetting
    if strcmp(cntr.PIDsetting, 'SD')
        %
        % Social distancing with effectiveness fixed at cntr.f, starting at
        % ton, stopping at toff: 
        %
        f = SD(t, cntr); 
        [~, y_rep, per_var] = Prep_Signals(t, x, cntr, per_var);
    elseif strcmp(cntr.PIDsetting, 'LS') 
        % 
        % Lightswitch approach, implemented according to cntr.setting
        % 
        [f, per_var, y_rep] = SD_lightswitch(t, x, cntr, per_var);
    elseif strcmp(cntr.PIDsetting, 'Log_quantized')
        % Quantized after PID
        %
        [f, per_var, y_rep] = PID_delayed(t, x, cntr, per_var);
        % Define quantization scale (i.e. levels of f that can be achieved)
        Qs = [0.05:0.15:1];
        [~, idx] = min(abs(f-Qs) );
        f = Qs(idx);
    else
        %
        % Implementation of PID control. If cntr.setting
        % is 'E_delayed' if 'I_delayed', feedback from the delayed
        % measurement is used. 
        %
        [f, per_var, y_rep] = PID_delayed(t, x, cntr, per_var);    
    end
    d = cntr.d(floor(t)+1); % Disturbance (updates once a day)
    p_vac = cntr.p_vac(floor(t)+1); % percentage vaccinated (per day)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = SD(t, cntr)
    % f = SD(t, cntr)
    %
    % Social distancing starting at t1, ramping up in linear fashion until
    % t2, to remain constant after t2. Arbitrary vectors of f can be passed
    % on too.
    %
    % Inputs: 
    % cntr.f        = value of f to use in simulation. If a vector is passed (length(cntr.f)>1),
    %               cntr.f should be a vector with a value of f for each
    %               second of simulation in settings.Tspan
    % cntr.ton      = Time that social distancing starts 
    % cntr.ton2     = Time where social distancing reaches it's full value 
    % cntr.toff     = Time social distancing is turned off
    % cntr.f_off    = (optional) value f reaches after cntr.toff
    %
    if length(cntr.f)>1
        f = cntr.f(floor(t)+1); 
    else
        tmp_f = [1:-(1-cntr.f)/(cntr.ton2-cntr.ton):cntr.f cntr.f];
        if t >= cntr.ton && t < cntr.ton2
            f = tmp_f(floor(t)-cntr.ton+2);
        elseif t >= cntr.ton2 && t < cntr.toff
            f = cntr.f;
        elseif isfield(cntr, 'f_off') && t >= cntr.toff
            f = cntr.f_off; 
        else
            f = 1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f, per_var, y_rep] = SD_lightswitch(t, x, cntr, per_var)
    % f = SD_lightswitch(t, cntr)
    %
    % Social distancing lightswitch approach, switching between 
    % effectiveness f and 30% of f. Switch on and off when at 80 and 30
    % percent of capacity. Capacity measures can be defined through
    % cntr.setting. ('I', 'I_delayed', 'E2_delayed', 'I_reported')
    %
    [y, y_rep, per_var] = Prep_Signals(t, x, cntr, per_var);
    
    if t < cntr.ton
        f = 1;
    elseif t <= cntr.ton + cntr.t_ini_fix;
        f = cntr.f;
    else  
        if y > 0.8 * cntr.threshold
            f = cntr.f; 
        elseif y < 0.2 * cntr.threshold
            f = 1-(1-cntr.f)*.3; 
        else
            f = per_var.f_kmin1;
        end
    end
    per_var.f_kmin1 = f;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [f, per_var, y_rep] = PID_delayed(t, x, cntr, per_var)
    % [f, per_var] = PID_delayed(t, x, cntr, per_var)
    % 
    % Implementation of PID control. 
    %
    % Choose the correct signal and delay for feedback, as defined in cntr.
    [y, y_rep, per_var, delay] = Prep_Signals(t, x, cntr, per_var);

    % Controller parameters
    kp = cntr.kp;
    ki = cntr.ki;
    kd = cntr.kd;

    % Discretize controller for the defined sampling time cntr.h
    ki = ki*cntr.h;
    kd = kd/cntr.h;
    aw = cntr.h/cntr.Tt; 
    
    % Setpoint
    if isfield(per_var, 'sp')
        sp = per_var.sp;
    elseif length(cntr.sp) > 1
        %
        % To change the setpoint during the simulation, cntr.sp needs
        % to be defined as a vector with the same length as Tspan. The
        % setpoint can be changed per day.
        %
        sp = cntr.sp(floor(t)+1);
    else
        sp = cntr.sp;
    end
    if length(cntr.f) > 1
        f = cntr.f(floor(t)+1);
    else
        f = cntr.f;
    end
    
    % Control action
    if t < cntr.ton
        % Initial period of uncontrolled growth
        f = 1;     
    elseif t >= cntr.ton && t < cntr.ton2
        % Gradual effect of initial lockdown (Mar-April 2020)
        tmp_f = [1:-(1-f)/(cntr.ton2-cntr.ton):f f];
        f = tmp_f(floor(t)-cntr.ton+2);
    elseif t >= cntr.ton2 && t < cntr.ton + cntr.t_ini_fix
        % Constant of predefined (open-loop) restrictions
        f = f;
    else  
        % Start feedback
        if strcmp(per_var.f_trigger,  'off')
            % Initialize, the controller starts one day later
            f = f;
            per_var.f_kmin1 = f;
            per_var.f_trigger = 'on';
            per_var.t_trigger = t;
            per_var.t_update_last = t;
        else
            if t - per_var.t_trigger >= cntr.h
                % This code will be called whenever the ODE solver is
                % called. Update control action only after the sampling
                % interval has passed. 
                f = per_var.I- kp*(y-sp) -  kd*(y-per_var.y_memory(end-delay-1));
                
                % Update the integrator and the memory for trigger times
                per_var.I = per_var.I- ki*(y-sp);
                per_var.t_trigger = t;
                per_var.t_update_last = t;
                % Constrain f T
                if isfield(cntr, 'MaxF')
                    fmax = cntr.MaxF;
                else
                    fmax = .8;
                end
                
                % Anti-windup on the integrator
                if f >fmax;
                    per_var.I = per_var.I + aw*(fmax-f);
                    f = fmax;
                elseif f < 0
                    per_var.I = per_var.I + aw*(0-f);
                    f = 0;
                end
                per_var.f_kmin1 = f;
            end
            % Apply the control action
            f = per_var.f_kmin1;
        end 
    end

    % Constrain f (This should be unnecessary)
    if f >1;
        f = 1;
    elseif f < 0
        f = 0;
    end

    
end

function [y, y_rep, per_var, delay] = Prep_Signals(t, x, cntr, per_var)
    % Define feedback signal
    if strcmp(cntr.setting, 'I_delayed')
        if length(x)> 12
            y_fb = sum(x([4, 9, 13, 17]));
        else
            y_fb = sum(x([4, 9]));
        end
        delay = cntr.delayI;   
    elseif strcmp(cntr.setting, 'E_delayed')
        if length(x)> 12
            y_fb = sum(x([3, 8, 12, 16]));
        else
            y_fb = sum(x([3, 8]));
        end
        delay = cntr.delayE;
    else % use I, not delayed
        y_fb = sum(x([4, 9]));
        delay = 0;
    end
    if any([strcmp(cntr.PIDsetting, 'Log_PI'), strcmp(cntr.PIDsetting, 'Log_quantized')])
        y_fb = log(y_fb);
    end
    
    if (t - per_var.t_update_last)>=1
        % add noise
        if isfield(cntr, 'noise_var')
            y_fb = y_fb+cntr.noise_var*randn(1,1);
                        
            if ~isfinite(y_fb)
                y_fb = 0;
            end

        end
        % Filter the measurement
        if isfield(cntr, 'filter')
            if ~isfinite(y_fb)
                y_fb = 0;
            end
            per_var.filt_u = [y_fb ;per_var.filt_u(1:end-1)];
            y_fb = cntr.filter.num{1}(2:end)*per_var.filt_u -cntr.filter.den{1}(2:end)*per_var.filt_y;
            per_var.filt_y = [y_fb ;per_var.filt_y(1:end-1)];
        end
        % Update the memory
        per_var.y_memory = [per_var.y_memory(2:end); y_fb];
        per_var.E2_memory = [per_var.E2_memory(2:end); sum(x([3, 8]))];
        per_var.t_update_last = t;
    end
    
    % Return the delayed feedback variable
    y = per_var.y_memory(end-delay); 
    y_rep =y;
end
