function [L,t_end,trecord] = set_parameters(multi)
% FUNCTION SET_PARAMETERS
% returns all source term and numerical solver parameters.

global N eta
global b gamma n p G_T l0 phi G_h beta
global diff diff_2 n1 n2 n3 n4
global dz dt max_iters err_tol z nodes tol

N = 20; % number of cells
eta = 1; % mobility coefficient 
L = 10; % length of tissue

%%% SOURCE TERM PARAMETERS
if multi == 1
    %%% parameters for turing pattens, Figure 6 (b,d,f,h).
    diff = 0.5;
    diff_2 = 5;
    
    % source term parameters from Equation (3.4,3.5)
    n1 = 0.1; n2 = 1; n3 = 0.5; n4 = 1;
    
    t_end = 100; % final time
    trecord = [10 20 40 90 t_end+1]; % times to record solution at

else
    %%% parameters for gtpases, Figure 3 (b,d,f,h)
    diff = 1;
    diff_2 = 1;
    
    % source term parameters from Equation (3.1)
    b = 0.2; gamma = 1.5; beta = 0.3; n = 4; G_T = 2;
    p = 4; l0 = 1; phi = 0.65; G_h = 0.4;
    
    t_end = 600; % final time
    trecord = [60 120 200 350 430 t_end+1]; % times to record solution at
end

%%% NEWTON SOLVER PARAMETERS 
dz = 0.001; % spatial step
dt = 1e-3; % temporal step
max_iters = 1000; % max iterations for Newton solver
err_tol = 1e-8; % error tolerance for Newton solver

z = 0:dz:1; % define fixed domain
z = z';

nodes = ceil(1/dz) + 1; % number of nodes in x direction
tol = 1e-3; % tolerance for time record


end