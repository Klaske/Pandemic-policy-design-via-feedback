%% Robustness illustrated with Monte Carlo simulation of 400 realizations
% For Anderson et al. model for BC. Feedback from hospitalizations( I with 14 day delay)
clear all; close all

nr_sims = 400;

cntr.setting = 'I_delayed'; % Note that the delay is in the controller. 
cntr.h = 14; 
cntr.PIDsetting =  'Log_PI'; % PI controller with feedback from log(x)
% Design controller
[cntr.kp, cntr.ki, cntr.kd] = Controller_Design(30, 14); 
cntr.I0       = 1-.78/3*2; % 2/3 of maximum as estimated for BC in March/ April 2020. 

% The simulations don't include scaling, only simulate x. 
% Setpoint for 20 hospitalizations assuming 7.3% is hospitalized. 
cntr.sp = log([1/.073*20*ones(30*7, 1) ;1/.073*20*ones(30*7,1)]);   

CV = Anderson_COVID_SEEIR(); 

cntr.t_ini_fix = 28+31+30-CV.cntr.ton-cntr.h+1; % Set the start date for closed-loop start of May.  
update.cntr = cntr; update.settings = settings;
CV = CV.Update_Settings(update);

% Monte carlo randomization
% Fix seed for reproducibility
s = rng(15);
% 
prm_names = {'D', 'k1', 'k2', 'q', 'R0'}
default_prms = [5 .2 1 .05 3];
W = logspace(-5, 2)
RDN_Split = randn(nr_sims,6);
stdvs = default_prms/5; % These are chosen as exaggerated uncertainty to illustrate the robustness principle.  

%% Simulate
for it = 1:nr_sims
    prm_values(it, :) = max(default_prms + stdvs.*RDN_Split(it, 1:5), [.01 .01 .01 .01 .01]);
    for nm = 1:length(prm_names)
        CV.prms.(prm_names{nm}) = prm_values(it, nm); 
    end
    CV.cntr.delayI = ceil(14 + 4*RDN_Split(it, 6));
    [CV, Tcl_log(:,it), Xcl_log(it, :, :), fcl_log(:,it)] = CV.Simulate('test');
end
% Result (hospitalizations, 7.3% of I)
I_outcl_log =  sum(Xcl_log(:, :, [4 9]), 3)*.073;

%% Open loop response
CV_OL = CV; 
CV_OL.cntr.PIDsetting = 'SD';
CV_OL.Res = [];
CV_OL.cntr.f = [ones(1, 44), 0.22*ones(1, 90-44), .61*ones(1, 366-90)];
CV_OL.cntr.f(1:70) = fcl_log(1:70, 1); 

for it = 1:nr_sims
    for nm = 1:length(prm_names)
        CV_OL.prms.(prm_names{nm}) = prm_values(it, nm);
    end
    CV_OL.cntr.delayI = ceil(14 + 4*RDN_Split(it, 6));
    [CV_OL, Tol_log(:,it), Xol_log(it, :, :), fol_log(:,it)] = CV_OL.Simulate('test');
end
% Result (hospitalizations, 7.3% of I)
I_outol_log =  sum(Xol_log(:, :, [4 9]), 3)*.073;   


 %% Plot results, plotting I instead of hospitalizations for consistency of presentation
figure
subplot(221)
x = [1:size(I_outol_log(:, :)', 1)]+datetime(2020,2,1); 
ys = prctile(I_outol_log(:, :)/.073, [0 25 50 75 100], 1)';
blue = [ 0    0.4470    0.7410]; 
orange = [0.8500    0.3250    0.0980]; 
% 0-100 percentile
plot_shaded(x, ys(:,1)', ys(:,5)', orange, 2); hold on; 
% 25-75 percentile
plot_shaded(x, ys(:,2)', ys(:,4)', orange, 1.3); hold on; 
plot(x, ys(:,3),'color', orange/1.5, 'linewidth', 2); hold on
ylim([0 25000]); xlims_ = xlim; 

subplot(223)
x2 = [1:size(I_outol_log(:, :)', 1)]+datetime(2020,2,1); 
ys2 = prctile(fol_log(:, :), [0 25 50 75 100], 2); 
% 0-100 percentile
plot_shaded(x2, ys2(:,1)', ys2(:,5)', orange, 2); hold on; 
% 25-75 percentile
plot_shaded(x2, ys2(:,2)', ys2(:,4)', orange, 1.3); hold on; 
plot(x2, ys2(:,3),'color', orange/1.5, 'linewidth', 2); hold on
ylim([0 1]); xlim(xlims_)

ys_cl = prctile(I_outcl_log(:, :)/.073, [0 25 50 75 100], 1)';
ys2_cl = prctile(fcl_log(:, :), [0 25 50 75 100], 2); 

subplot(222)
plot_shaded(x, ys_cl(:,1)', ys_cl(:,5)', blue, 2); hold on; 
% 25-75 percentile
plot_shaded(x, ys_cl(:,2)', ys_cl(:,4)', blue, 1.3); hold on; 
plot(x, ys_cl(:,3),'color', blue/1.5, 'linewidth', 2); hold on
ylim([0 25000]); xlim(xlims_)

subplot(224)
plot_shaded(x2, ys2_cl(:,1)', ys2_cl(:,5)', blue, 2); hold on; 
% 25-75 percentile
plot_shaded(x2, ys2_cl(:,2)', ys2_cl(:,4)', blue, 1.3); hold on; 
plot(x2, ys2_cl(:,3),'color', blue/1.5, 'linewidth', 2); hold on
ylim([0 1]); xlim(xlims_)

subplot(221); grid on; 
 ylabel('Number of Infections (I_T)');
 
subplot(223); grid on;
 ylabel('Level of activity u(t)');
subplot(222); grid on; 
 ylabel('Number of Infections (I_T)');

subplot(224); grid on;
 ylabel('Level of activity u(t)');
 set(gcf, 'color', [1 1 1]); 
 
%% This figure showing percentiles does not tell the whole story, as the
% different models peak at different times. 
% outputs I_outol_log and I_outcl_log are 7.3% of infections
% (hosptializations): nr_models_over_1000 => over 1000/.073 =~13699 cases
% nr_models_cl_over_40 => over 40/0.073 =~ 548 cases

nr_models_over_1000 = length(find((max(I_outol_log(:,5*30:end)')>1000)))
perc_models_over_1000 = length(find((max(I_outol_log(:,5*30:end)')>1000)))/nr_sims*100
perc_models_over_500 = length(find((max(I_outol_log(:,5*30:end)')>500)))/nr_sims*100
nr_models_over_500 = length(find((max(I_outol_log(:,5*30:end)')>500)))

nr_models_cl_over_40 = length(find((max(I_outcl_log(:,5*30:end)')>40)))
perc_models_cl_over_40 = length(find((max(I_outcl_log(:,5*30:end)')>40)))/nr_sims*100