%% Properties of the system
% This file simulates the SEEIR model and generates the results presented in
% Fig. 1. 
%
%
%% Find steady state f_0 
cntr.setting = 'I_delayed'; % In ODE model this is different as the delay is in the controller. 
cntr.h = 1; % If the objective (desired response) is faster, the sampling time needs to correspond. 
% For comparison and simplicity we use daily updates for all examples here.
cntr.PIDsetting =  'Log_PI'; 

% The simulations don't include scaling, only simulate x. 
% Setpoint for 400 
cntr.sp = log([400*ones(366*2, 1)]);   

CV = Anderson_COVID_SEEIR(); 

% Design controller, as for robust design with 14 days delay
[cntr.kp, cntr.ki, cntr.kd] = Controller_Design(30, 14); 
cntr.I0       = .55;% To achieve close to steady state

cntr.ton = 1; 
cntr.ton2 = 1;  
settings.Tspan = [0:2*365];
cntr.f = ones(size(settings.Tspan))*1; 
cntr.t_ini_fix = 34; % Settings to achieve steady state before disturbance, for comparison 
update.cntr = cntr; update.settings = settings;
CV = CV.Update_Settings(update);

% Simulate
[CV, Tcl_log, Xcl_log, fcl_log] = CV.Simulate('test');

x0 = Xcl_log(end, :); 
f_0 = fcl_log(end);
%% Simulations
% Define the three inputs for scenarios
input0 = zeros(1,125); 
input1 = [zeros(1,9) .24*ones(1,25) -.12*ones(1,50) zeros(1,41)];
input2 = [zeros(1,9) -.12*ones(1,50) .24*ones(1,25) zeros(1,41)];

CV.cntr.PIDsetting = 'SD'; % implement as feedforward social distancing. 
% set x0 
CV.settings.x0 = x0
CV.settings.Tspan = 1:124
% Constant u = f_0
CV.cntr.f = f_0; 
% Simulate
[CV, T_constant, X_constant, f_constant] = CV.Simulate('constant');
I_constant =  sum(X_constant(:, [4 9]),2);
% f = f_0 + input 1
CV.cntr.f = f_0+input1; 
% Simulate
[CV, T_input1, X_input1, f_input1] = CV.Simulate('input1');
I_input1 =  sum(X_input1(:, [4 9]),2);
% f = f_0 + input 2
CV.cntr.f = f_0+input2; 
% Simulate
[CV, T_input2, X_input2, f_input2] = CV.Simulate('input2');
I_input2 =  sum(X_input2(:, [4 9]),2);


%% Calculate total cases and interventions
Results_totalcases = {['\int I(t)  = ' num2str(floor(sum(I_constant)))] ;['\int I(t)  = ' num2str(floor(sum(I_input1)))]; ['\int I(t)  = ' num2str(floor(sum(I_input2)))]}
Results_totalinterventions = {['\int u_e(t)  = ' num2str(round(sum(f_constant-f_0), 1))] ;['\int u_e(t)  = ' num2str(round(sum(f_input1-f_0),1)) ];['\int u_e(t)  = ' num2str(round(sum(f_input2-f_0),1))]}
%% Results
figure('position', [400 400 800 250])

subplot(321); plot(I_constant, 'color', [0 0.4470 0.7410], 'linewidth', 2); hold on; grid on;  ylabel('Number of infections I_T'); ylim([0 2000])
subplot(323); plot(I_input1, 'color', [0.8500 0.3250 0.0980], 'linewidth', 2) ; hold on; grid on;  ylabel('Number of infections I_T'); ylim([0 2000])
subplot(325); plot(I_input2, 'color', [0.9290 0.6940 0.1250], 'linewidth', 2);  hold on; grid on;  ylabel('Number of infections I_T'); ylim([0 2000])
xlabel('Time [days]');

subplot(322); plot(f_constant, 'color', [0 0.4470 0.7410], 'linewidth', 2); hold on; grid on; ylabel('Intervention u(t)'); ylim([0 1])
subplot(324); plot(f_input1, 'color', [0.8500 0.3250 0.0980], 'linewidth', 2); hold on; grid on; ylabel('Intervention u(t)'); ylim([0 1])
subplot(326); plot(f_input2, 'color', [0.9290 0.6940 0.1250], 'linewidth', 2); hold on; grid on; ylabel('Intervention u(t)'); ylim([0 1])
xlabel('Time [days]');
set(gcf, 'color', [1 1 1])

