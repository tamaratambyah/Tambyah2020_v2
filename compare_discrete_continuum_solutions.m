function compare_discrete_continuum_solutions(x_midpoints,discrete_density,cont_density,q_IC,trecord,s_all,L,axis_name)
% FUNCTION COMPARE_DISCRETE_CONTINUUM_SOLUTIONS
% returns a figure comparing discrete and continuum solutions 

global z 
figure
hold on

%%% PLOT CONTINUUM SOLUTIONS
plot(z.*L,q_IC, 'LineWidth', 1.2, 'DisplayName', '0'); % initial condition 

% solutions in trecord
for i = 1:length(trecord)-1
    plot(z.*s_all(:,i),cont_density(:,i), 'LineWidth', 1.2, 'DisplayName', num2str(trecord(i)))
end

ax = gca;
ax.ColorOrderIndex = 1;

%%% PLOT DISCRETE SOLUTIONS
plot(x_midpoints(1,:),discrete_density(1,:), 'o','LineWidth', 1, 'DisplayName', [num2str(0)]) % initial condition

% solutions in trecord
for i = 1:length(trecord)-1
    plot(x_midpoints(trecord(i)/0.1+1,:), discrete_density(trecord(i)/0.1+1,:), 'o','LineWidth', 1, 'DisplayName', num2str(trecord(i)))
end


xlabel('\it x','Interpreter', 'latex')
ylabel(['\it ' axis_name],'Interpreter', 'latex')
title(['Discrete-continuum comparison: \it ' axis_name],'Interpreter', 'latex')
set(gca,'FontName', 'Times New Roman')  % Set it to times
set(gca,'FontSize', 16)
box on

%%% CONSTRUCT LEGEND
time_leg = [0 trecord];

for i = 1:length(time_leg)-1
    legendInfo{i} = ['{\it t }= ' num2str(time_leg(i))] ;
end

legend(legendInfo,'FontSize', 16,'FontName', 'Times New Roman')


end