function plot_characteristics(t_plot,x_plot,variable,N,name)
% FUNCTION PLOT_CHARACTERISTICS
% returns the plot of characteristics from the discrete model


figure; 
hold on 
% plot cell boundaries
for i = 1:N+1
    plot3(x_plot(:,i), t_plot(:,i), 20*ones(size(x_plot(:,i))), ...
        'k', 'LineWidth', 2)
end

% plot concentration of chemical 
for i = 1:N
    surf(x_plot(:,i:i+1), t_plot(:,i:i+1),[variable(:,i) variable(:,i)], ...
        'EdgeColor', 'none')
    view(2)
end

h = colorbar;
ylabel(h, ['\it ' name],'Interpreter', 'latex')
xlabel('\it x','Interpreter', 'latex')
ylabel('\it t','Interpreter', 'latex')
title(['Discrete characteristics: \it ' name],'Interpreter', 'latex')

set(gca,'FontName', 'Times New Roman')  % Set it to times
set(gca,'FontSize', 16)
box on
xlim([0 inf])


end
