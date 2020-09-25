function run_discrete_simulation(L,multi,t_end)
% FUNCTION RUN_DISCRETE_SIMULATION 
% returns the solution to Equations (2.1-2.11)
fprintf('\nRunning discrete simulation\n')

global N 
%% Inital conditions 
x_IC = linspace(0,L,N+1)'; % uniform initial condition for position 
C1_IC = ones(N,1); % uniform initial condition for chemical 
C2_IC = 0.5.*ones(N,1);

%% solve ode

if multi == 1
    [t, x_sol] = ode15s(@(t,x) ode_func_turing(t,x), [0:0.1:t_end], [x_IC; C1_IC; C2_IC]);
else
    [t, x_sol] = ode15s(@(t,x) ode_func_gtpase(t,x), [0:0.1:t_end], [x_IC; C1_IC; C2_IC]);
end

% extract positions and C solutions 
x_plot = x_sol(:,1:N+1);
discrete_C1 = x_sol(:,N+2:2*N+1);
discrete_C2 = x_sol(:,2*N+2:end);
%% plotting variables

discrete_density = zeros(length(t),N);
for i = 1:N
    discrete_density(:,i) = 1./((x_plot(:,i+1) - x_plot(:,i)));
end

t_plot = zeros(size(x_plot));
for i = 1:size(t_plot,2)
   t_plot(:,i) =  t;
end

x_midpoints = zeros(size(x_plot,1),N);
for i = 1:size(x_plot,1)
   x_midpoints(i,:) = get_transport_rates_and_intervals(x_plot(i,:));
end

%% plot characteristics
plot_characteristics(t_plot,x_plot,discrete_C1,N,'C1')
colormap cool

plot_characteristics(t_plot,x_plot,discrete_density,N,'q')

if multi == 1
    plot_characteristics(t_plot,x_plot,discrete_C2,N,'C2')
    colormap winter
    save('discrete_solution_turing');
else
    save('discrete_solution_gtpase');
end




end
