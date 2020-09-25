function run_continuum_simulation(L,t_end,trecord,multi)
% FUNCTION RUN_CONTINUUM_SIMULATION
% returns the solution to Equations (2.32-2.37)
fprintf('\nRunning continuum simulation\n')

global N eta % cellular parameters
global diff diff_2 % diffusivities
global n1 n2 n3 n4 % source term parameters for Turing patterns
global b beta gamma n G_T phi G_h p l0 % source term parameters for gtpase
global dz dt max_iters err_tol z nodes tol % newton solver parameters

%% CHEMICAL SOURCE TERMS

% from Equation (3.4,3.5)
if multi == 1 % turing patterns
    C1_ode = @(q,k,a,G,H) + n1 - n2*G + n3*G^2*H;
    derv_C1_ode = @(q,k,a,G,H) -n2 + n3*2*G*H;
    C2_ode = @(q,k,a,H,G) n4 - n3*G^2*H;
    derv_C2_ode = @(q,k,a,H,G) -n3*G^2;
    derv_force = @(q,k,a,G,H) 0;
    a_func = @(G) 1 + 0.*G;
    k_func = @(G) 1 + 0.*G;
    C2_IC = 0.5.*ones(length(z),1); % uniform initial condition for C2
else % gtpase from Equation (3.1)
    C1_ode = @(q,k,a,G,H) (b + beta*(1./q - a) + gamma*G^n/(1+G^n))*(G_T-G) - G;
    derv_C1_ode = @(q,k,a,G,H) (gamma*n*G^(n-1)/(1+G^n)^2 + beta*phi*p*G^(p-1)*G_h^p/(G_h^p+G^p)^2  )*(G_T - G)...
        -(b + beta*(1./q - a) + gamma*G^n/(1+G^n)) -1;
    derv_force = @(q,k,a,G,H) 0.05*1/q - 0.05*a - k*-phi.*p.*G_h.^p.*G.^(p-1)./(G_h.^p + G.^p).^2;
    C2_ode = 0;
    derv_C2_ode = 0;
    a_func = @(G) l0 - phi.*G.^p./(G_h.^p + G.^p);
    k_func = @(G) 1 + 0.05.*G;
    C2_IC = zeros(length(z),1); % uniform initial condition for C2
end

%% INITIAL CONDITIONS

q_IC = N/L+0*z; % uniform density initial condition
C1_IC = ones(length(z),1); % uniform chemical initial condition for C1
s_IC = L; % initial interface positions
a_IC = a_func(C1_IC);
k_IC = k_func(C1_IC);

%% INITIALISE STORAGE OF SOLUTIONS
q_all = zeros(length(q_IC), length(trecord));
k_all = zeros(length(k_IC), length(trecord));
a_all = zeros(length(k_IC), length(trecord));
C1_all = zeros(length(k_IC), length(trecord));
C2_all = zeros(length(k_IC), length(trecord));
s_all = zeros(1,length(trecord));

%% NEWTON-RAPHSON METHOD

tic
t = 0; % start time
record = 1; % counter to record solution

% pass variable
q_update = q_IC;
s_update = s_IC;
C1_update = C1_IC;
C2_update = C2_IC;
k_update = k_IC;
a_update = a_IC;

while t < t_end % time loop
    t = t + dt;
    q_old = q_update;
    s_old = s_update;
    C1_old = C1_update;
    C2_old = C2_update;
    
    %%% ITERATIONS
    for w = 1:max_iters
        q = q_update;
        s = s_update;
        C1 = C1_update;
        C2 = C2_update;
        k = k_update;
        a = a_update;
        
        %%% SOLVE FOR Q
        [L,D,U] = jacobian(q,k,nodes,dt,a,s_old,s);
        rhs = -discretised_func(q,q_old,k,a,nodes,dt,s_old,s);
        dq = tridi(L,D,U,rhs); % solve solution using thomas algorthim
        q_update = q + dq;
        res_q = norm(dq,inf); % update residual
        
        %%% SOLVE FOR C1
        [L,D,U] = jacobian_g(q_update,k,a,C1,dt,s,s_old,derv_C1_ode,derv_force,diff,C2,multi,1);
        rhs = -discretised_g(q_update,k,a,C1,C1_old,dt,s,s_old,C1_ode,diff,C2,multi,1);
        dC = tridi(L,D,U,rhs);
        C1_update = C1 + dC;
        res_C1 = norm(dC,inf);    % update residual
        
        %%% SOLVE FOR C2
        [L,D,U] = jacobian_g(q_update,k,a,C2,dt,s,s_old,derv_C2_ode,derv_force,diff_2,C1_update,multi,0);
        rhs = -discretised_g(q_update,k,a,C2,C2_old,dt,s,s_old,C2_ode,diff_2,C1_update,multi,0);
        dC = tridi(L,D,U,rhs);
        C2_update = C2 + dC;
        res_C2 = norm(dC,inf);    % update residual
        
        %%% COMPUTE K and A
        a_update = a_func(C1_update);
        k_update = k_func(C1_update);
        
        
        %%% SOLVE FOR TISSUE LENGTH
        s_update = s_old + dt*1/(eta*q_update(end)*s_old*dz)*(k_update(end)*(1/q_update(end)-a_update(end)) ...
            - k_update(end-1)*(1/q_update(end-1)-a_update(end-1)));
        
        % check if error tolerances met
        if res_q <= err_tol && res_C1 <= err_tol && res_C2 <= err_tol
            iters = w;
            break;
        end
        
    end
    
    % record solution
    if abs(trecord(record) - t) <= tol
        q_all(:,record) = q_update;
        k_all(:,record) = k_update;
        a_all(:,record) = a_update;
        C1_all(:,record) = C1_update;
        C2_all(:,record) = C2_update;
        s_all(:,record) = s_update;
        record = record + 1;
    end
    
    if mod(t,1) < tol
        fprintf('\nt = %f; iters = %f\n', t,iters);
    end
    
end
toc

if multi == 1
    save('continuum_solution_turing')
else
    save('continuum_solution_gtpase')
end
end

