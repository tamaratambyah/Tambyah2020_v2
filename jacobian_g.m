function [LL,DD,UU] = jacobian_g(q,k,a,G,dt,s,s_old,...
    derv_G_ode,derv_force_G,diff,H,multi,is_C1)
global z dz eta nodes

%%% INTIALISE OUTPUT
LL = zeros(nodes,1); % lower diagonal
DD = zeros(nodes,1); % diagonal
UU = zeros(nodes,1); % upper diagonal

if multi == 0 && is_C1 == 0
    DD = ones(nodes,1); % diagonal
    return;
end

ds_dt = (s-s_old)/dt;

% Neumann boundary condition at x = 0, Equation (2.31)
j = 1;
DD(j) = 1;
UU(j) = -1;


% internal nodes, Equation (2.33)
for j = 2:nodes-1
    
    %%% NON-UPWINDING COMPONENT
    D = 1/dt ...
        + 1/(eta*s^2*dz^2*q(j))*( k(j-1)*(1/q(j-1)-a(j-1)) - 2*k(j)*(1/q(j)-a(j)) + k(j+1)*(1/q(j+1)-a(j+1)) ...
        -1/(4*q(j))*(q(j+1)-q(j-1))*(k(j+1)*(1/q(j+1)-a(j+1))-k(j-1)*(1/q(j-1)-a(j-1)))) ...
        + G(j)/(eta*s^2*dz^2*q(j))*( -2*derv_force_G(q(j),k(j),a(j),G(j),H(j))) ...
        + 2*diff/(s^2*dz^2) ...
        - derv_G_ode(q(j),k(j),a(j),G(j),H(j)) ;
    U = + G(j)/(eta*s^2*dz^2*q(j))*( derv_force_G(q(j+1),k(j+1),a(j+1),G(j+1),H(j+1)) -1/(4*q(j))*(q(j+1)-q(j-1))*derv_force_G(q(j+1),k(j+1),a(j+1),G(j+1),H(j+1))) ...
        -diff/(s^2*dz^2) ;
    L = G(j)/(eta*s^2*dz^2*q(j))*(derv_force_G(q(j-1),k(j-1),a(j-1),G(j-1),H(j-1)) -1/(4*q(j))*(q(j+1)-q(j-1))*-derv_force_G(q(j-1),k(j-1),a(j-1),G(j-1),H(j-1))) ...
        -diff/(s^2*dz^2);
    
    %%% UPWINDING
    v = 1/(eta*q(j)*s^2*2*dz)*(k(j+1)*(1/q(j+1)-a(j+1)) - k(j-1)*(1/q(j-1)-a(j-1))) ...
        -z(j)/s*ds_dt;
    
    if v > 0
        D_ad = 1/dz*(1/(eta*q(j)*s^2*2*dz)*(k(j+1)*(1/q(j+1)-a(j+1)) - k(j-1)*(1/q(j-1)-a(j-1))) ...
            -z(j)/s*ds_dt );
        L_ad =  -1/dz*(1/(eta*q(j)*s^2*2*dz)*(k(j+1)*(1/q(j+1)-a(j+1)) - k(j-1)*(1/q(j-1)-a(j-1))) ...
            -z(j)/s*ds_dt   ) ...
            + (G(j)-G(j-1))/dz*1/(eta*q(j)*s^2*2*dz)*-derv_force_G(q(j-1),k(j-1),a(j-1),G(j-1),H(j-1));
        U_ad = (G(j)-G(j-1))/dz*1/(eta*q(j)*s^2*2*dz)*derv_force_G(q(j+1),k(j+1),a(j+1),G(j+1),H(j+1));
        
    elseif v < 0
        D_ad = -1/dz*(1/(eta*q(j)*s^2*2*dz)*(k(j+1)*(1/q(j+1)-a(j+1)) - k(j-1)*(1/q(j-1)-a(j-1))) ...
            -z(j)/s*ds_dt );
        L_ad = (G(j+1)-G(j))/dz*1/(eta*q(j)*s^2*2*dz)*-derv_force_G(q(j-1),k(j-1),a(j-1),G(j-1),H(j-1));
        U_ad = 1/dz*(1/(eta*q(j)*s^2*2*dz)*(k(j+1)*(1/q(j+1)-a(j+1)) - k(j-1)*(1/q(j-1)-a(j-1))) ...
            -z(j)/s*ds_dt  ) ...
            + (G(j+1)-G(j))*(1/(eta*q(j)*s^2*2*dz)*derv_force_G(q(j+1),k(j+1),a(j+1),G(j+1),H(j+1)));
        
    else
        D_ad = 0;
        U_ad = 0;
        L_ad = 0;
    end
    
    DD(j) = D + D_ad;
    LL(j) = L + L_ad;
    UU(j) = U + U_ad;
    
end


% Neumann boundary condition at x = L(t), Equation (2.31)
j = nodes;
v = 1/(s^2*eta*q(j)*dz)*(k(j)*(1/q(j)-a(j))-k(j-1)*(1/q(j-1)-a(j-1)))  -z(j)/s*ds_dt;
DD(j) = 1/dt ...
    + v*1/dz ...
    + (G(j)-G(j-1))/dz*( 1/(s^2*eta*q(j)*dz)*derv_force_G(q(j),k(j),a(j),G(j),H(j))) ...
    + 1/(s^2*dz^2*eta)*(1/q(end)-1/q(end-1))*(k(j)*(1/q(j)-a(j))-k(j-1)*(1/q(j-1)-a(j-1))) ...
    + G(j)/(s^2*dz^2*eta)*(1/q(end)-1/q(end-1))*derv_force_G(q(j),k(j),a(j),G(j),H(j)) ...
    - derv_G_ode(q(j),k(j),a(j),G(j),H(j)) ...
    - diff/(s^2*dz^2)*-2;
LL(j) = v*-1/dz ...
    + (G(j)-G(j-1))/dz*( 1/(s^2*eta*q(j)*dz)*-derv_force_G(q(j-1),k(j-1),a(j-1),G(j-1),H(j-1))) ...
    + G(j)/(s^2*dz^2*eta)*(1/q(end)-1/q(end-1))*-derv_force_G(q(j-1),k(j-1),a(j-1),G(j-1),H(j-1)) ...
    - diff/(s^2*dz^2)*2;



end