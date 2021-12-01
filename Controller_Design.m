function [kp, ki, kd] = Controller_Design(varargin)
    %% Controller is tuned using an approximate model identified (estimated) from 
    % a series of step responses on the nonlinear model. The approximate model
    % describes the relation between f (input) and log(y), where y =
    % Infections (I_T). The delay can be defined (varargin), the nominal model assumes feedback is delayed by 14 days.
    % This response can be approximated by an integral system with lag
    % and delay. 
    %
    % The approximate model is then used to design a PI controller using
    % Skogestad's SMC rules for an integral system with lag and delay.

    %% ID from series of step responses
    % This method is motivated by model reduction methods using a simulated
    % experiment and system identification (estimation of the lower order
    % model). An advantage of this model is that frequency weighting on the
    % model fit can be achieved either by filtering or by the choice of
    % excitation signal. In this case, rather than using a random excitation,
    % we choose a stepwise excitation to allow for validation of the
    % approximation compared to the nonlinear model. 

    % Excitation signal 
    rndnr = [1 .2 .4 .6 .34 .12 .98 .65 .45 .58 .32 .76 .56 .28];
    f_ID = ones(1, 40)*rndnr(1); 
    for ij = 2:length(rndnr)
    f_ID = [f_ID rndnr(ij)*ones(1,40)];
    end

    %% Simulate Anderson et al model with this input sequence
    % Instance of simulation class
    CV = Anderson_COVID_SEEIR(); 
    % Update some settings
    cntr.setting = 'SD'; cntr.PIDsetting = 'SD';
    cntr.toff = 51+31+15;
    cntr.f_off = CV.cntr.f*2.5; 
    cntr.f = f_ID;
    settings.Tspan = 0:360;

    update.cntr = cntr;
    update.settings = settings;
    CV = CV.Update_Settings(update);

    % Simulate nominal model
    [CV, Tnom, Xnom, fnom] = CV.Simulate('test');

    %% Use (non-delayed) data for estimation of approximate model (identification)
    % Response of infected compartments to excitation signal f_ID
    Inom = sum(Xnom( :, [4 9]), 2);
    offset = 1;

    % Output is log of I
    y = log(Inom);

    % System is know to contain integral action. To get an accurate model and
    % accurate integrator, identify the model using the derivative and add the
    % integrator afterwards. Force start at zero to avoid initial condition
    % problems. This is a constant that can be added again later. 
    yd = diff(y)-max(diff(y));

    % Define data in format for SysID toolbox
    dat = iddata(yd(20:end), f_ID(19:359)'-offset, 1);

    % Can force zero initial conditions, not necessary in this case as
    % we forced everything to start at zero. Note that these offsets don't
    % affect the dynamics, only the initial conditions. 
    % opt_k = oeOptions('InitialCondition','zero');
    % Identify model using output error approach
    Moe = oe(dat,[1, 1, 1]);

    %% Validate model 
    % Generate response for validation
    y_m = lsim(Moe, f_ID-offset); 
    % Add initial condition to simulation
    y_m_I = lsim(Moe*tf(1, [1 -1], 1), f_ID-offset)+lsim(tf(max(diff(y)), [1 -1], 1), ones(size(f_ID)))+.885; 

    % Validate model
    figure
    subplot(211)
    plot(y(20:end), 'r', 'linewidth', 2); hold on; plot(y_m_I(20:360), 'b', 'linewidth', 2);grid on
    xlabel('Time [days]'); ylabel('ln(I)')
    subplot(212)
    plot(yd(20:end)+max(diff(y)), 'r', 'linewidth', 2); hold on; plot(y_m(20:360)+max(diff(y)), 'b', 'linewidth', 2);grid on
    xlabel('Time [days]'); ylabel('dln(I) / dt')
    set(gcf, 'color', [1 1 1])

    %% Model parameters, SMC method uses continuous time model parameters
    Mcontinuous = d2c(Moe); 
    k = Mcontinuous.b(1)/Mcontinuous.f(2); 
    tau2 = 1/Mcontinuous.f(2);

    % Design variable 
    if isempty(varargin)
        tau_c = 30; % (conservative design robust to uncertainty, robustness analysis not included in paper or this material)
        delay = 14;
    else
        tau_c = varargin{1}; 
        delay = varargin{2};
    end
    CV = CV.SMCController( k, tau2, delay, tau_c, 'PI');

    kp = CV.cntr.kp; ki = CV.cntr.ki; kd = CV.cntr.kd; 
end