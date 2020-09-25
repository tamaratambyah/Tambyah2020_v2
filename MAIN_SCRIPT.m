% Filename: MAIN_SCRIPT.m
% Author: Tamara A. Tambyah, September 2020
% Queensland University of Technology, Brisbane, Australia 
% This function generates Figure 3 (b,d,f,h) and Figure 6 (b,d,f,h) in the 
% manuscript: "A free boundary mechanobiological model of epithelial tissues"
% 
% WARNING: The continuum simulation is computationally expensive and takes
% some time to run.
% To test the continuum simulation and to see example outputs,
% change t_end = 5 and trecord = [1 3 5 t_end+1] in set_parameters.m. This
% simulation takes approximately 10 mintues. 
% To generate Figure 3 (b,d,f,h) and Figure 6 (b,d,f,h), run the code as
% is. This simulation will take a couple of hours. 

clear all
close all 

%% RHO GPTASE
% generates Figure 3 (b,d,f,h)

[L,t_end,trecord] = set_parameters(0); % get all parameters 
run_discrete_simulation(L,0,t_end); % run discrete simulation 
run_continuum_simulation(L,t_end,trecord,0); % run continuum simulation

% compare discrete and continuum solutions for density and chemical
load('discrete_solution_gtpase')
load('continuum_solution_gtpase')
compare_discrete_continuum_solutions(x_midpoints,discrete_density,q_all,q_IC,trecord,s_all,s_IC,'q')
compare_discrete_continuum_solutions(x_midpoints,discrete_C1,C1_all,C1_IC,trecord,s_all,s_IC,'C1')

%% TURING PATTERN
% generates Figure 6 (b,d)

[L,t_end,trecord] = set_parameters(1); % get source term parameters 
run_discrete_simulation(L,1,t_end); % run discrete simulation 
run_continuum_simulation(L,t_end,trecord,1); % run continuum simulation

% compare discrete and ccontinuum solutions for density and chemical 
load('discrete_solution_turing')
load('continuum_solution_turing')
compare_discrete_continuum_solutions(x_midpoints,discrete_density,q_all,q_IC,trecord,s_all,s_IC,'q')
compare_discrete_continuum_solutions(x_midpoints,discrete_C1,C1_all,C1_IC,trecord,s_all,s_IC,'C1')
compare_discrete_continuum_solutions(x_midpoints,discrete_C2,C2_all,C2_IC,trecord,s_all,s_IC,'C2')

