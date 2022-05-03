# Pandemic-policy-design-via-feedback
Code to generate results in "Effective pandemic policy design through feedback does not need accurate predictions" K.van Heusden et al. 

This code was developed in Matlab 2019b and has not been tested in any other version. 

This code is based on a Matlab implementation of the model by Anderson et al., which is available here
https://github.com/carolinecolijn/distancing-impact-covid19.git
Alterations to the original model have been made to incorporate feedback and additional scenarios, as detailed in the paper. 

main.m generates all figures in the paper: 

% To generate figure 1, run SystemPropertiesSEEIR
% To generate figure 3, run Uncertainty_Final
    % This includes Monte Carlo simulations and may take a while to run
% To generate figure 4, run Delay_Final
% To generate figure 5, run InputUncertainty_Final
% To generate figure 6, run VaccinVariant_BCrate
% To generate figure 7 in the appendix, run Step_Illustration

Simulations use a Matlab implementation of SEEIR model Anderson et al. for BC. 
in Anderson_COVID_SEEIR.m
- Implemented as ODE, equivalent to published version
- Using ode45 solver

Time-varying signals are included in the differential equation by calling the function
[f, d, ~, p_vac] = COVID_Control_PID(t, x, cntr); in the definition of the differential equations. 
This is required for feedback control. 

The controller COVID_Control_PID is implemented as a discrete controller with
sampling time of h (in days). 

Note that this combination is a hybrid system. We chose to implement it in ODE form
 to allow the code to maintain a structure commonly used in epidemiology and 
easier to read/ follow for experts in epidemiological modeling, 
rather than using a hybrid toolbox or simulink. 
- In this implementation, accurate sampling of the controller is
not ensured. The timing of the update of the control action depends on
the simulation time for the ODE. For the purpose of illustration of principles 
of feedback control, this is sufficient and can be made arbirarily accurate by limiting
the step size of the solver. 
- The ODE solver may step back in time, which can cause issues in a controller with 
integral action when called and coded in this form if not taken into account.
To avoid these issues, the output of the controller (COVID_Control_PID.m) 
depends on the (absolute) simulation time t. 

Controller_Design.m
- estimates an approximate model from a series of step responses on the nonlinear model
- determines PI(D) controller tuning using Skogestad's SMC method. 
