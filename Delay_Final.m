% Combined delay and sampling figure 
% Start with delay
close all; clear all; clc
%% Illustrate impact of time delays using nominal model SEEIR. 

cntr.setting = 'I_delayed'; % The delay is in the controller. 
cntr.h = 1; % Set sampling time
cntr.PIDsetting =  'Log_PI'; 

% The simulations don't include scaling, only simulate x. 
% Setpoint for 20 hospitalizations assuming 15% is hospitalized. More
% recent estimates are ~7.3% hospitalizations, for which this setpoint
% would correspond to ~ 10 (9.733) hospitalizations. Number of infections as reported are
% the same. 
cntr.sp = log([1/.15*20*ones(366*2, 1)]);   
% Disturbance of 100 exposed on day 450 (in steady-state starting point)
cntr.d = zeros(450, 1); cntr.d(450) = 100; 

CV = Anderson_COVID_SEEIR(); 

% Design controller, as for robust design with 14 days delay
[cntr.kp, cntr.ki, cntr.kd] = Controller_Design(30, 14); 
cntr.I0       = .55;% To achieve steady state prior to disturbance

cntr.ton = 1; 
cntr.ton2 = 1;  
settings.Tspan = [0:2*365];
cntr.f = ones(size(settings.Tspan))*1; 
cntr.t_ini_fix = 24; % Settings to achieve steady state before disturbance, for comparison 
update.cntr = cntr; update.settings = settings;
CV = CV.Update_Settings(update);

% Simulate
[CV, Tcl_log, Xcl_log, fcl_log] = CV.Simulate('test');
% Result 
I_outcl_log =  sum(Xcl_log(:, [4 9]),2);

% Reduce time delay 
CV.cntr.delayI = 2; 
cntr.t_ini_fix = 26;% To achieve steady state prior to disturbance
update.cntr = cntr; 
CV = CV.Update_Settings(update);

% Design a controller for the reduced delay and simulate 
[cntr.kp, cntr.ki, cntr.kd] = Controller_Design(15, 2); 
cntr.t_ini_fix = 24;% To achieve steady state prior to disturbance
update.cntr = cntr; 
CV = CV.Update_Settings(update);
% Simulate
[CV, Tcl_log_design_sd, Xcl_log_design_sd, fcl_log_design_sd] = CV.Simulate('test');
% Result 
I_outcl_log_design_sd =  sum(Xcl_log_design_sd( :, [4 9]),2);


%% Plot Results Infectious (no delay)
figure(22)
subplot(221)
plot([I_outcl_log], 'color', [0 0.4470 0.7410], 'linewidth', 2); hold on
plot([I_outcl_log_design_sd], 'color', [0.9290 0.6940 0.1250], 'linewidth', 2)
xlim([450 550]); ylim([ 100 200]); grid on
xticks(430:20:550); xticklabels({'0', '0', '20','40','60','80', '100'}); 
xlabel('Time [days]')
ylabel('Number of infections (I_T)')

subplot(223)
plot(fcl_log, 'color', [0 0.4470 0.7410], 'linewidth', 2); hold on
plot(fcl_log_design_sd, 'color', [0.9290 0.6940 0.1250], 'linewidth', 2)
ylabel('Intervention u(t)')
xlim([450 550]); ylim([0.45 .6]); grid on
xticks(430:20:550); xticklabels({'0','0', '20','40','60','80', '100'}); 
xlabel('Time [days]')
set(gcf, 'color', [1 1 1]); 

Results= {['\int I(t) slow  = ' num2str(floor(sum(I_outcl_log(450:550))))] ;['\int I(t) faster  = ' num2str(floor(sum(I_outcl_log_design_sd(450:550))))]; ...
	['\int slow u_e(t)  = ' num2str(round(sum(fcl_log(450:550)-fcl_log(450)), 1))] ;['\int faster u_e(t)  = ' num2str(round(sum(fcl_log_design_sd(450:550)-fcl_log(450)), 1))]}

Results_outbreak= {['\int I(t) slow  = ' num2str(floor(sum(I_outcl_log(450:550)-I_outcl_log(450))))] ;['\int I(t) faster  = ' num2str(floor(sum(I_outcl_log_design_sd(450:550)-I_outcl_log_design_sd(450))))]; ...
	['\int slow u_e(t)  = ' num2str(round(sum(fcl_log(450:550)-fcl_log(450)), 1))] ;['\int faster u_e(t)  = ' num2str(round(sum(fcl_log_design_sd(450:550)-fcl_log(450)), 1))]}
%% Illustrate impact of sampling time using nominal model SEEIR. 
cntr.setting = 'I_delayed'; % Delay is in the controller
colored = 2
for h = [1 7 14 28] %loop for different sampling times
    cntr.h = h; 
    cntr.PIDsetting =  'Log_PI'; 

    % The simulations don't include scaling, only simulate x. 
    % Setpoint for 20 hospitalizations assuming 15% is hospitalized. More
    % recent estimates are ~7.3% hospitalizations, for which this setpoint
    % would correspond to ~ 10 (9.733) hospitalizations. Number of infections as reported are
    % the same.
    cntr.sp = log([1/.15*20*ones(366*2, 1)]);   
    % Disturbance of 100 exposed on day 450 (in steady-state starting point)
    cntr.d = zeros(470, 1); cntr.d(450) = 100; 

    CV = Anderson_COVID_SEEIR(); 

    % Design controller, as for robust design with 14 days delay
    [cntr.kp, cntr.ki, cntr.kd] = Controller_Design(30, 14); 
    cntr.I0       = .55;% To achieve steady state prior to disturbance

    cntr.ton = 1; 
    cntr.ton2 = 1;  
    settings.Tspan = [0:2*365];
    cntr.f = ones(size(settings.Tspan))*1; 
    cntr.t_ini_fix = 25-cntr.h; % Settings to achieve steady state before disturbance, for comparison 
    update.cntr = cntr; update.settings = settings;
    CV = CV.Update_Settings(update);

    % Simulate
    [CV, Tcl_log, Xcl_log, fcl_log] = CV.Simulate('test');
    % Result 
    I_outcl_log =  sum(Xcl_log(:, [4 9]),2);

    % Reduce time delay 
    CV.cntr.delayI = 2; 
    update.cntr = cntr; 
    CV = CV.Update_Settings(update);

    % Design a controller for the reduced delay and simulate 
    [cntr.kp, cntr.ki, cntr.kd] = Controller_Design(15, 2); 
    cntr.t_ini_fix = 25-cntr.h;% To achieve steady state prior to disturbance
    update.cntr = cntr; 
    CV = CV.Update_Settings(update);
    % Simulate
    [CV, Tcl_log_design_sd, Xcl_log_design_sd, fcl_log_design_sd] = CV.Simulate('test');
    % Result 
    I_outcl_log_design_sd =  sum(Xcl_log_design_sd( :, [4 9]),2);


%% Plot Results Infectious (no delay) only fast design
if colored == 1
        if h == 1
            clr = [0.9290 0.6940 0.1250]
            wdth = 2
            stle = '-'
        elseif h == 7 
            clr = [0.8500 0.3250 0.0980]
            wdth = 2
            stle = '-'
        elseif h == 14
            clr = [0.4940 0.1840 0.5560]
            wdth = 2
            stle = '-'
        else
            clr = [0.4660 0.6740 0.1880]
            wdth = 2
            stle = '-'
        end
elseif colored == 2
        if h == 1
            clr = [0.9290 0.6940 0.1250]
            wdth = 2
            stle = '-'
        elseif h == 7 
            clr = [0.9290 0.6940 0.1250]
            wdth = 2
            stle = ':'
        elseif h == 14
            clr = [0.9290 0.6940 0.1250]
            wdth = 2
            stle = '-.'
        else
            clr = [0.9290 0.6940 0.1250]
            wdth = 2
            stle = '--'
        end
end
        figure(22)
        subplot(222)
        
        plot([I_outcl_log_design_sd], 'color', clr, 'linewidth', wdth, 'linestyle', stle); hold on
        xlim([450 650]); ylim([ 100 200]); grid on
        xticks(450:40:650); xticklabels({'0', '40','80',  '120', '160','200'}); 
        xlabel('Time [days]')
        ylabel('Number of infections (I_T)')

        subplot(224)
        
        stairs(fcl_log_design_sd, 'color', clr, 'linewidth', wdth, 'linestyle', stle); hold on
        ylabel('Intervention u(t)')
        xlim([450 650]); ylim([0.45 .6]); grid on
        xticks(450:40:650); xticklabels({'0', '40','80',  '120', '160','200'}); 
        xlabel('Time [days]')
        set(gcf, 'color', [1 1 1]);


end