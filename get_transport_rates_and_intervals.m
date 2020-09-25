function [x_midpoints,lengths,T_L,T_R,H_L,H_R] = get_transport_rates_and_intervals(x)
% FUNCTION GET_TRANSPORT_RATES_AND_INTERVALS 
% returns the transport rates from Equations (2.4-2.7), and computes the
% voronoi partition

global  N diff diff_2

%%% COMPUTE VORONOI PARTITION. Refer to Supplementary Material Section 2ai
a = ones(N);
a(1) = 0;
b = ones(N);
c = zeros(N);
d = 2.*x(1:end-1);
d(1) = (x(2)+x(1))/2;

x_midpoints = tridi(a,b,c,d);

%%% COMPUTE CELL LENGTHS 
lengths = zeros(N,1);
for i = 1:N
   lengths(i) = x(i+1)-x(i); 
end

%%% COMPUTE LEFT AND RIGHT TRANSPORT RATES FOR CHEMICAL 1
T_R = zeros(N,1); % initalise right transport rates
T_L = zeros(N,1); % intialise left transport rates 

% internal cells from Equation (2.5,2.6)
for i = 2:N-1
    den = (x_midpoints(i)-x_midpoints(i-1))+(x_midpoints(i+1)-x_midpoints(i));
    T_L(i) = 2*diff/((x_midpoints(i)-x_midpoints(i-1))*den);
    T_R(i) =  2*diff/((x_midpoints(i+1)-x_midpoints(i))*den);
    
end

% cell 1 from Equation(2.7)
den = 2*x_midpoints(1)+(x_midpoints(2)-x_midpoints(1));
T_R(1) = 2*diff/((x_midpoints(2)-x_midpoints(1))*den);

% cell N from Equation (2.8)
den = (x_midpoints(end)-x_midpoints(end-1)) + 2*(x(N+1) - x_midpoints(end));
T_L(end) = 2*diff/((x_midpoints(end)-x_midpoints(end-1))*den) ; 


%%% COMPUTE LEFT AND RIGHT TRANSPORT RATES of CHEMICAL 2
H_R = zeros(N,1); % initalise right transport rates
H_L = zeros(N,1); % intialise left transport rates 

% internal cells from Equation (2.5,2.6)
for i = 2:N-1
    den = (x_midpoints(i)-x_midpoints(i-1))+(x_midpoints(i+1)-x_midpoints(i));
    H_L(i) = 2*diff_2/((x_midpoints(i)-x_midpoints(i-1))*den);
    H_R(i) =  2*diff_2/((x_midpoints(i+1)-x_midpoints(i))*den);
    
end

% cell 1 from Equation(2.7)
den = 2*x_midpoints(1)+(x_midpoints(2)-x_midpoints(1));
H_R(1) = 2*diff_2/((x_midpoints(2)-x_midpoints(1))*den);

% cell N from Equation (2.8)
den = (x_midpoints(end)-x_midpoints(end-1)) + 2*(x(N+1) - x_midpoints(end));
H_L(end) = 2*diff_2/((x_midpoints(end)-x_midpoints(end-1))*den) ; 
end