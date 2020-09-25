function [L,D,U] = jacobian(q,k,nodes,dt,a,s_old,s)
% FUNCTION JACOBIAN 
% returns the jacobian relating in Equation (2.32) in tridiagonal form 

global dz z eta

%%% INTIALISE OUTPUT
L = zeros(nodes,1); % lower diagonal
D = zeros(nodes,1); % diagonal
U = zeros(nodes,1); % upper diagonal 

r = dt/(eta*dz^2*s_old^2); % constant 
ds_dt = (s-s_old)/dt; % change in s

% boundary condition at x = 0, Equation (2.34)
D(1) = k(1)/q(1)^2;
U(1) = -k(2)/q(2)^2;

% internal nodes, Equation (2.32)
for i = 2:length(z)-1
    C = z(i)*dt/(2*s_old*dz);
    D(i) = -1 - 2*r*k(i)*1/q(i)^2;
    L(i) = r*k(i-1)*1/q(i-1)^2 - C*ds_dt;
    U(i) = r*k(i+1)*1/q(i+1)^2 + C*ds_dt;
end

% boundary condition at x = L(t), Equation (2.35)
i = nodes;
D(i) = 1/(2*s_old*dz)*(-1/q(i)^2*(k(i)*(1/q(i)-a(i)) - k(i-1)*(1/q(i-1)-a(i-1))) ...
    - 1/q(i)*k(i)/q(i)^2) - k(i)/q(i)^2 ...
    +(k(i)/8*(-1/q(i)^2*(1/q(i)*(a(i)-a(i-1)) - 1/q(i-1)*(a(i)-a(i-1)) ) + 1/q(i)*(a(i)-a(i-1))*-1/q(i)^2) ...
    + 1/(4)*(a(i)-a(i-1))*(k(i)-k(i-1))*-2/q(i)^3 ...
    + a(i)/8*(-1/q(i)^2*(1/q(i)*(k(i)-k(i-1)) - 1/q(i-1)*(k(i)-k(i-1)) ) + 1/q(i)*(k(i)-k(i-1))*-1/q(i)^2) )*1/(s_old^2*dz^2);
L(i) = 1/(2*q(i)*s_old*dz)*k(i-1)/q(i-1)^2 ...
    + (k(i)/(8*q(i))*-(a(i)-a(i-1))*-1/q(i-1)^2 + a(i)/(8*q(i))*-(k(i)-k(i-1))*-1/q(i-1)^2)*1/(s_old^2*dz^2);


end