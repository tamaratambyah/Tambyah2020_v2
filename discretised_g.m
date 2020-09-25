function g = discretised_g(q,k,a,G,G_old,dt,s,s_old,G_ode,diff,H,multi,is_C1)
global z dz eta nodes
% FUNCTION DISCRETISED_G
% returns the discretised version of Equation (2.33)


g = zeros(nodes,1); % initialise output
if multi == 0 && is_C1 == 0
    return;
end

ds_dt = (s-s_old)/dt;

% Neumann boundary condition at x = 0, Equation (2.31)
j = 1;
g(j) = G(j) - G(j+1);

% internal nodes, Equation (2.33)
for j = 2:nodes-1
    
    %%% UPWINDING
    v = 1/(eta*q(j)*s^2*2*dz)*(k(j+1)*(1/q(j+1)-a(j+1)) - k(j-1)*(1/q(j-1)-a(j-1))) ...
        -z(j)/s*ds_dt;
    
    if v > 0
        advect = G(j)-G(j-1);
    elseif v < 0
        advect = G(j+1)-G(j);
    else
        advect = 0;
    end
    
    g(j) = (G(j) - G_old(j))/dt ...
        + v*advect/dz ...
        + G(j)/(eta*s^2*dz^2*q(j))*( k(j-1)*(1/q(j-1)-a(j-1)) - 2*k(j)*(1/q(j)-a(j)) + k(j+1)*(1/q(j+1)-a(j+1)) ...
        -1/(4*q(j))*(q(j+1)-q(j-1))*(k(j+1)*(1/q(j+1)-a(j+1))-k(j-1)*(1/q(j-1)-a(j-1)))) ...
        -diff/(s^2*dz^2)*(G(j+1)-2*G(j)+G(j-1)) ...
        -G_ode(q(j),k(j),a(j),G(j),H(j));
        
end

% Neumann boundary condition at x = L(t), Equation (2.31)
j = nodes;
v = 1/(s^2*eta*q(j)*dz)*(k(j)*(1/q(j)-a(j))-k(j-1)*(1/q(j-1)-a(j-1)))  -z(j)/s*ds_dt;
g(j) = (G(j)-G_old(j))/dt ...
    + v*(G(j)-G(j-1))/dz ...
    + G(j)/(eta*s^2*dz^2)*(1/q(j)-1/q(j-1))*(k(j)*(1/q(j)-a(j))-k(j-1)*(1/q(j-1)-a(j-1))) ...
    -diff/(s^2*dz^2)*(2*G(j-1)-2*G(j)) ...
    - G_ode(q(j),k(j),a(j),G(j),H(j));



end