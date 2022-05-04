close all; clear all;clc
%% Simulate appearance of new variant, vaccination and combination of both with PI controller 

cntr.setting = 'I_delayed'; % The ODE model does not contain any delays. Delays are implemented in the controller. 
cntr.h = 1; % For illustration, daily updates are used. 
cntr.PIDsetting = 'Log_PI'; 
cntr.MaxF = 1; % Allow full opening after vaccination
% The simulations don't include scaling, only simulate x. 
% Setpoint for 20 hospitalizations assuming 15% is hospitalized. More
% recent estimates are ~7.3% hospitalizations, for which this setpoint
% would correspond to ~ 10 (9.733) hospitalizations. Number of infections as reported are
% the same
cntr.sp = log([1/.15*20*ones(366*3, 1)]);   

CV = Anderson_COVID_SEEIR(); 
cntr.delayI = 2; 

% Controller parameters as default, simulate initial BC peak in April 2020.
% Design controller 
[cntr.kp, cntr.ki, cntr.kd] = Controller_Design(15, 2); 
update.settings.Tspan = [0:365*3];
update.cntr = cntr; 
CV = CV.Update_Settings(update);

%% Variant
CV.cntr.p_vac = zeros(size(CV.settings.Tspan)); 
CV.prms.Iv_in = 5; 
%CV.prms.R0v = CV.prms.R0*1.4; % 140% , UK variant
%CV.prms.Iv_t = 326; % Starting in December (UK variant) 
CV.prms.R0v = CV.prms.R0*2.5; % 250% , Delta 
%CV.prms.Iv_t = 326 + 55; % Starting in Feb 2021 (P1 variant) 
CV.prms.Iv_t = 326 + 55 +60 -15; % Starting mid March 2021 (Delta variant) 
[CV, T_Var, X_Var, f_Var] = CV.Simulate('test');
% Result (I --> Infections at each time t)
I_V1 = sum(X_Var(:, [4 9]),2); 
I_V2 = sum(X_Var(:, [13 17]),2); 
I_Var =  I_V1+I_V2;

%% Vaccination
CV.cntr.p_vac = zeros(size(CV.settings.Tspan)); 
importvaccinationrate
% t0 = Feb 1st 2020
tstart = days(BCVaccinationprogress.week_end(3)-datetime(2020, 2,1)-7)
% create vector with daily vacination rate
vacs = repelem(BCVaccinationprogress.weeklyrate(3:end)/7,7);
%
CV.cntr.p_vac(tstart+21:tstart+21+length(vacs)-1) = vacs/100; 
CV.prms.Iv_in = 0; 

[CV, T_Vac, X_Vac, f_Vac] = CV.Simulate('test');
% Result (I --> Infections at each time t)
I_Vac =  sum(X_Vac(:, [4 9 13 17]),2);

%% Vaccination and variant 
CV.prms.Iv_in = 5; 
[CV, T_Vac_Var, X_Vac_Var, f_Vac_Var] = CV.Simulate('test');
% Result (I --> Infections at each time t)
I_Vac_V1 = sum(X_Vac_Var(:, [4 9]),2); 
I_Vac_V2 = sum(X_Vac_Var(:, [13 17]),2); 
I_Vac_Var = I_Vac_V1 + I_Vac_V2;

%% Plot results V4
figure
l_orange = [0.9290, 0.6940, 0.1250]; 
d_orange = [0.8500 0.3250 0.0980]; 
l_blue = [0.3010, 0.7450, 0.9330];
t_var = CV.prms.Iv_t % Introduction of variant
t_vac = tstart %start of vaccinations
subplot(231)
% plot_shaded(1:length(I_V1)+2, [zeros(2, 1); I_V1]', zeros(length(I_V1)+2, 1)', d_orange, 2.5); hold on
% plot_shaded(1:length(I_V2)+2, [zeros(2, 1); I_V2]', zeros(length(I_V2)+2, 1)', d_orange, 1.5)
plot([I_V1], 'color',min(d_orange*2.5, [1 1 1]), 'linewidth', 2); hold on
plot([ I_V2],'color', min(d_orange*1.5, [1 1 1]), 'linewidth', 2)

plot([ I_Var], 'color', d_orange, 'linewidth', 2); hold on
plot([t_var t_var], [0 1000], 'color', l_blue, 'linestyle', '--')

xlim([0 2*360]); ylim([0 1000]); grid on
xlabel('Time [days]')
ylabel('Number of infections (I_T)')
title('Variant (2.5 times more contagious)')
subplot(234)
plot(f_Var, 'color', d_orange, 'linewidth', 2); hold on
plot([t_var t_var], [0 1], 'color', l_blue, 'linestyle', '--')
ylabel('Level of activity u(t)')
xlim([0 2*360]); ylim([0 1]); grid on
xlabel('Time [days]')
set(gcf, 'color', [1 1 1]); 

subplot(232)
%plot_shaded(1:length(I_Vac)+2, [zeros(2, 1); I_Vac]', zeros(length(I_Vac)+2, 1)', d_orange, 2.5); hold on
plot([ I_Vac], 'color', d_orange, 'linewidth', 2); hold on
xlim([0 2*360]); ylim([0 1000]); grid on
plot([t_vac t_vac], [0 1000], 'color', l_blue*.5, 'linestyle', '--')
xlabel('Time [days]')
ylabel('Number of infections (I_T)')
title('Vaccinations')

subplot(235)
plot(f_Vac, 'color', d_orange, 'linewidth', 2); hold on
ylabel('Level of activity u(t)')
plot([t_vac t_vac], [0 1], 'color', l_blue*.5, 'linestyle', '--')
xlim([0 2*360]); ylim([0 1]); grid on
xlabel('Time [days]')
set(gcf, 'color', [1 1 1]); 

subplot(233)
%plot_shaded(1:length(I_Vac_V1)+2, [zeros(2, 1); I_Vac_V1]', zeros(length(I_Vac_V1)+2, 1)', d_orange, 2.5); hold on
%plot_shaded(1:length(I_Vac_V2)+2, [zeros(2, 1); I_Vac_V2]', zeros(length(I_Vac_V2)+2, 1)', d_orange, 1.5)
plot([ I_Vac_V1], 'color', min(d_orange*2.5, [1 1 1]), 'linewidth', 2); hold on
plot([ I_Vac_V2], 'color', min(d_orange*1.5, [1 1 1]), 'linewidth', 2); hold on
plot([I_Vac_Var], 'color', d_orange, 'linewidth', 2); hold on
xlim([0 2*360]); ylim([0 1000]); grid on
plot([t_var t_var], [0 1000], 'color', l_blue, 'linestyle', '--')
plot([t_vac t_vac], [0 1000], 'color', l_blue*.5, 'linestyle', '--')
xlabel('Time [days]')
ylabel('Number of infections (I_T)')
title('Variant and vaccinations') 

subplot(236)
plot(f_Vac_Var, 'color', d_orange, 'linewidth', 2); hold on
ylabel('Level of activity u(t)')
xlim([0 2*360]); ylim([0 1]); grid on
plot([t_var t_var], [0 1], 'color', l_blue, 'linestyle', '--')
plot([t_vac t_vac], [0 1], 'color', l_blue*.5, 'linestyle', '--')
xlabel('Time [days]')
set(gcf, 'color', [1 1 1]); 

