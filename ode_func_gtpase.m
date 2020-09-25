function f = ode_func_gtpase(t,x)
% FUNCTION ODE_FUNC_GTPASE
% returns the ode function for the discrete model for Case study 2, Equations (2.1-2.10) 

global N 
global eta 
global b gamma n G_T beta l0 G_h p phi

feedback = @(x2,x1,a) x2-x1-a;

[x_midpoints,lengths,T_L,T_R] = get_transport_rates_and_intervals(x(1:N+1));

C = x(N+2:end); % extract chemical

% compute resting length
a =  l0 - phi.*C.^p./(G_h.^p + C.^p);
k = 1 + 0.05.*C;
  
% initalise output
f = zeros(3*N+1,1);

%%% ODES FOR CELL BOUNDARIES
f(1) = 0; % fixed boundary at x = 0 from Equation (2.4)
for i = 2:N
    f(i) = 1/eta*(k(i)*(x(i+1) - x(i) - a(i)) - k(i-1)*(x(i) - x(i-1) - a(i-1)));
    % linear force law in Equation (2.3) 
end
f(N+1) = -1/eta*k(N)*(x(N+1)-x(N)-a(N)); % free boundary at x = L(t) from Equation (2.4)

%%% ODES FOR CHEMICAL
% internal cells from Equation (2.10)
for i = N+3:2*N
    j = i - (N+1);
    f(i) = (b +beta*(feedback(x(j+1),x(j),a(j))) + gamma*x(i)^n/(1+x(i)^n))*(G_T - x(i)) - x(i)...
        + 1/lengths(j)*(T_R(j-1)*x(i-1)*lengths(j-1) ...
        - (T_R(j) + T_L(j))*x(i)*lengths(j) ...
        + T_L(j+1)*x(i+1)*lengths(j+1)) ...
        -  x(i)*( f(j+1)-f(j))/lengths(j);
end

% boundary cell 1 from Equation (2.9)
i = N+2;
j = i - (N+1);
f(i) = (b +beta*(feedback(x(j+1),x(j),a(j))) + gamma*x(i)^n/(1+x(i)^n))*(G_T - x(i)) - x(i)...
   + 1/lengths(1)*(T_L(2)*x(i+1)*lengths(j+1) -T_R(1)*x(i)*lengths(j)) ...
    - x(i)*( f(j+1)-f(j))/lengths(j);

% boundary cell N from Equation (2.11)
i = 2*N+1;
j = i - (N+1);
f(i) = (b +beta*(feedback(x(j+1),x(j),a(j))) + gamma*x(i)^n/(1+x(i)^n))*(G_T - x(i)) - x(i)...
    +  1/lengths(end)*(T_R(j-1)*x(i-1)*lengths(j-1)- T_L(j)*x(i)*lengths(j)) ...
    - x(i)*( f(j+1)-f(j))/lengths(j);


end
