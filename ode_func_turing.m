function f = ode_func_turing(t,x)
% FUNCTION ODE_FUNC_TURING
% returns the ode function for the discrete model for Case study 3, Equations (2.1-2.10) 

global N n1 n2 n3 n4
global eta 

[x_midpoints,lengths,T_L,T_R,H_L,H_R] = get_transport_rates_and_intervals(x(1:N+1));

% compute resting length
a = ones(N,1); k = ones(N,1);

% initalise output
f = zeros(3*N+1,1);


%%% ODES FOR CELL BOUNDARIES
f(1) = 0; % fixed boundary at x = 0 from Equation (2.4)
for i = 2:N
    f(i) = 1/eta*(k(i)*(x(i+1) - x(i) - a(i)) - k(i-1)*(x(i) - x(i-1) - a(i-1)));
    % linear force law in Equation (2.3) 
end
f(N+1) = -1/eta*k(N)*(x(N+1)-x(N)-a(N)); % free boundary at x = L(t) from Equation (2.4)

%%% ODES FOR CHEMICAL 1
% internal cells from Equation (2.10)
for i = N+3:2*N
    j = i - (N+1);
    h = i + N;
    f(i) = (n1 - n2*x(i) + n3*x(i)^2*x(h)) .....
    + 1/lengths(j)*(T_R(j-1)*x(i-1)*lengths(j-1) ...
        - (T_R(j) + T_L(j))*x(i)*lengths(j) ...
        + T_L(j+1)*x(i+1)*lengths(j+1)) ...
        -  x(i)*( f(j+1)-f(j))/lengths(j);
end

% boundary cell 1 from Equation (2.9)
i = N+2;
j = i - (N+1);
h = i + N;
f(i) = (n1 - n2*x(i) + n3*x(i)^2*x(h)).....
+ 1/lengths(1)*(T_L(2)*x(i+1)*lengths(j+1) -T_R(1)*x(i)*lengths(j)) ...
    - x(i)*( f(j+1)-f(j))/lengths(j);

% boundary cell N from Equation (2.11)
i = 2*N+1;
j = i - (N+1);
h = i + N;
f(i) = (n1 - n2*x(i) + n3*x(i)^2*x(h))...
+  1/lengths(end)*(T_R(j-1)*x(i-1)*lengths(j-1)- T_L(j)*x(i)*lengths(j)) ...
    - x(i)*( f(j+1)-f(j))/lengths(j);


%%% ODES FOR CHEMICAL 2
% internal cells from Equation (2.10)
for i = 2*N+3:3*N
    j = i - 2*N -1; % x
    h = i - N; % G
    f(i) = (n4 - n3*x(h)^2*x(i)) .....
        + 1/lengths(j)*(H_R(j-1)*x(i-1)*lengths(j-1) ...
        - (H_R(j) + H_L(j))*x(i)*lengths(j) ...
        + H_L(j+1)*x(i+1)*lengths(j+1)) ...
        -  x(i)*( f(j+1)-f(j))/lengths(j);
end

% boundary cell 1 from Equation (2.9)
i = 2*N+2;
j = i - 2*N - 1;
h = i - N; 
f(i) = (n4 - n3*x(h)^2*x(i) )....
   + 1/lengths(1)*(H_L(2)*x(i+1)*lengths(j+1) -H_R(1)*x(i)*lengths(j)) ...
    - x(i)*( f(j+1)-f(j))/lengths(j);

% boundary cell N from Equation (2.11)
i = 3*N+1;
j = i - 2*N - 1;
h = i - N; 
f(i) = (n4 - n3*x(h)^2*x(i) )....
    +  1/lengths(end)*(H_R(j-1)*x(i-1)*lengths(j-1)- H_L(j)*x(i)*lengths(j)) ...
    - x(i)*( f(j+1)-f(j))/lengths(j);





end


