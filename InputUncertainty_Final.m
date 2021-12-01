close all; clear all; clc
%% Simulate with PI controller 

cntr.setting = 'I_delayed'; % The ODE model does not contain any delays. Delays are implemented in the controller. 
cntr.h = 1; % If the objective (desired response) is faster, the sampling time needs to correspond. 
% For comparison and simplicity we use daily updates for all examples here.
cntr.PIDsetting = 'Log_PI'; 

% The simulations don't include scaling, only simulate x. 
% Setpoint for 20 hospitalizations assuming 15% is hospitalized. More
% recent estimates are ~7.3% hospitalizations, for which this setpoint
% would correspond to ~ 10 (9.733) hospitalizations. Number of infections as reported are
% the same
cntr.sp = log([1/.15*20*ones(366*2, 1)]);   

CV = Anderson_COVID_SEEIR(); 
cntr.delayI = 2; 

% Controller settings, initial values
cntr.I0       = .548;% To achieve steady state prior to disturbance
cntr.ton = 1; 
cntr.ton2 = 1;  
settings.Tspan = [0:400];
%
%   The controlled system is a hybrid system. In the implementation it is not
% implemented as a hybrid system and accurate sampling of the controller is
% not ensured. The timing of the update of the control action depends on
% the simulation time for the ODE. For the purpose of closed-loop control
% and illustration of principles, this is sufficient. This allows the code
% to maintain a structure commonly used in epidemiology, in ODE form,
% rather than using a hybrid toolbox or simulink. 
%   In this simulation, an open-loop simulation is used to match the
% initial conditions of the Proportional controller to the PI response.
% With maxstep = 0.2 this fit is reasonably accurate. 
% 
settings.maxstep = 0.2; 
cntr.f = ones(size(settings.Tspan))*1; 

% To achieve steady state prior to disturbance
cntr.t_ini_fix = 27;
update.cntr = cntr; update.settings = settings;
CV = CV.Update_Settings(update);

% Design controller 
[cntr.kp, cntr.ki, cntr.kd] = Controller_Design(15, 2); 
update.cntr = cntr; 
CV = CV.Update_Settings(update);

% Introduce input uncertainty
CV.prms.u_offset = [zeros(1,200) .075*ones(1,280)];

[CV, T_PI, X_PI, f_PI] = CV.Simulate('test');
% Result (I --> Infections at each time t)
I_PI =  sum(X_PI(:, [4 9]),2);

% Simulate quantized PI controller
CV.cntr.PIDsetting = 'Log_quantized';%
[CV, T_PIQ, X_PIQ, f_PIQ] = CV.Simulate('test');
% Result (I --> Infections at each time t)
I_PIQ =  sum(X_PIQ(:, [4 9]),2);

%% Simulate proportional controller
% Response up to change similar, using open-loop simulation
cntr.I0       = f_PI(200);
cntr.f(1:204) = f_PI(1:204); 
cntr.t_ini_fix = 199;
% Use same proportional gain, take out integral action
cntr.ki = 0; 
cntr.Tt = inf;
update.cntr = cntr; 
CV = CV.Update_Settings(update);

[CV, T_P, X_P, f_P] = CV.Simulate('test');
% Result (I --> Infections at each time t)
I_P =  sum(X_P(:, [4 9]),2);

% Simulate quantized proportional controller
cntr.I0       = f_PIQ(200);
cntr.f(1:204) = f_PIQ(1:204); 
cntr.PIDsetting = 'Log_quantized';
update.cntr = cntr; 
CV = CV.Update_Settings(update);

[CV, T_PQ, X_PQ, f_PQ] = CV.Simulate('test');
% Result (I --> Infections at each time t)
I_PQ =  sum(X_PQ(:, [4 9]),2);

%% Plot results
figure

subplot(221)
plot([I_P], 'color', [0.4940    0.1840    0.5560], 'linewidth', 2); hold on
plot([I_PI], 'color', [0.9290 0.6940 0.1250], 'linewidth', 2); hold on
xlim([150 400]); ylim([100 450]); grid on
xticks(150:50:400); xticklabels({'0', '50','100','150','200', '250'}); 
xlabel('Time [days]')
ylabel('Number of infections (I_T)')

plot([200 200], [0 500], 'color',  [ .3 .3 .3], 'linestyle', '--')

subplot(223)
plot(f_P, 'color', [0.4940    0.1840    0.5560], 'linewidth', 2); hold on
plot(f_PI, 'color', [0.9290 0.6940 0.1250], 'linewidth', 2); hold on
ylabel('Intervention u(t)')
xlim([150 400]); ylim([0 1]); grid on
xticks(150:50:400); xticklabels({'0', '50','100','150','200', '250'}); 
xlabel('Time [days]')
plot([200 200], [0 1], 'color',   [ .3 .3 .3], 'linestyle', '--')

subplot(222)
plot([I_PIQ], 'color', [0.9290 0.6940 0.1250], 'linewidth', 2); hold on
plot([I_PQ], 'color', [0.4940    0.1840    0.5560], 'linewidth', 2)
xlim([150 400]); ylim([100 450]); grid on
xticks(150:50:400); xticklabels({'0', '50','100','150','200', '250'}); 
xlabel('Time [days]')
ylabel('Number of infections (I_T)')
plot([200 200], [0 500], 'color',  [ .3 .3 .3], 'linestyle', '--')

subplot(224)
plot(f_PIQ, 'color', [0.9290 0.6940 0.1250], 'linewidth', 2); hold on
plot(f_PQ, 'color', [0.4940    0.1840    0.5560], 'linewidth', 2)
ylabel('Intervention u(t)')
xlim([150 400]); ylim([0 1]); grid on
xticks(150:50:400); xticklabels({'0', '50','100','150','200', '250'}); 
xlabel('Time [days]')
plot([200 200], [0 1], 'color',  [ .3 .3 .3], 'linestyle', '--')
set(gcf, 'color', [1 1 1]); 
PI_AUC = sum(I_PIQ(150:400))
P_AUC = sum(I_PQ(150:400))

% Paper states that \int I PIQ  is 85% higher than \int I PQ. This is
% unclear, it is > 7x higher, or \int I PQ ~ 15% of \int I PIQ.
Results= {['\int I(t) PI  = ' num2str(floor(sum(I_PI(150:400)-I_PI(150))))] ;['\int I(t) P  = ' num2str(floor(sum(I_P(150:400)-I_P(150))))]; ...
	['\int I PIQ  = ' num2str(floor(sum(I_PIQ(150:400)-I_PIQ(150))))] ;['\int I PQ = ' num2str(floor(sum(I_PQ(150:400)-I_PQ(150))))];...
    ['\int v(t) PI  = ' num2str(round(sum(f_PI(150:400)-f_PI(150)), 1))] ;['\int v(t) P  = ' num2str(round(sum(f_P(150:400)-f_P(150)), 1))]; ...
	['\int v PIQ  = ' num2str(round(sum(f_PIQ(150:400)-f_PIQ(150)), 1))] ;['\int v PQ = ' num2str(round(sum(f_PQ(150:400)-f_PIQ(150)), 1))]}

