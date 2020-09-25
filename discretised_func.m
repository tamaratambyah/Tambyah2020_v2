function g = discretised_func(q,q_old,k,a,nodes,dt,s_old,s)
% FUNCTION DISCRETISED_FUNC 
% returns the discretised verision of Equation (2.32)

global dz z eta
g = zeros(nodes,1); % initialise output

r = dt/(eta*dz^2*s_old^2); % constant 
ds_dt = (s-s_old)/dt; % change in s

% boundary condition at x = 0, Equation (2.34)
g(1) = k(2)*(1/q(2)-a(2)) - k(1)*(1/q(1)-a(1)); 

% internal nodes, Equation (2.32)
for i = 2:length(z)-1
    C = z(i)*dt/(2*s_old*dz);
    g(i) = - q(i) + q_old(i) ...
        -r*(k(i-1)*(1/q(i-1)-a(i-1)) - 2*k(i)*(1/q(i)-a(i)) + k(i+1)*(1/q(i+1)-a(i+1)))...
        + C*ds_dt*(q(i+1)-q(i-1)) ;
end

% boundary condition at x = L(t), Equation (2.35)
i = nodes;
g(i) = 1/(2*q(i)*s_old*dz)*(k(i)*(1/q(i)-a(i)) - k(i-1)*(1/q(i-1)-a(i-1))) ...
    + k(i)*(1/q(i)-a(i)) ...
    + (k(i)/(8*q(i))*(1/q(i)*(a(i)-a(i-1)) - 1/q(i-1)*(a(i)-a(i-1)) ) ...
    + 1/(4*q(i)^2)*(a(i)-a(i-1))*(k(i)-k(i-1)) ...
    + a(i)/(8*q(i))*(1/q(i)*(k(i)-k(i-1)) - 1/q(i-1)*(k(i)-k(i-1)) ) )*1/(s_old^2*dz^2);




end