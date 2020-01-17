function [S_sol] = socp(P,S_sol,Ad,Bd,Cd,zd,Delta,i)

% Exceptions.
%#ok<*EQEFF>
%#ok<*VUNUS>
%#ok<*IDISVAR>
%#ok<*NODEF>

% Begin CVX.
cvx_begin
cvx_solver(P.cvx_solver)
cvx_precision(P.cvx_precision)
cvx_quiet(P.cvx_quiet)

% CVX variables.
variable s(1,P.K) nonnegative
variable x(P.n_x,P.K)
variable u(P.n_u,P.K)
variable nu(P.n_x,P.K-1)
variable eta(P.n_obs,P.K) nonnegative

% Cost function.
minimize(sum(s)+P.SCvx.lambda*norm(nu(:),1)+P.SCvx.lambda*sum(eta(:)))

% Trust region.
dx = x-S_sol(i).x;
du = u-S_sol(i).u;
norm([dx(:); du(:)],1) <= Delta;

% Time loop.
for k=1:P.K
    % Local variables.
    s_k = s(k);
    x_k = x(:,k);
    u_k = u(:,k);
    eta_k = eta(:,k);
    
    x0_k = S_sol(i).x(:,k);
    u0_k = S_sol(i).u(:,k);
    
    dx_k = dx(:,k);
    du_k = du(:,k);
    
    if (k < P.K)
        nu_k = nu(:,k);
        
        x_kp1 = x(:,k+1);
        u_kp1 = u(:,k+1);
        
        x0_kp1 = S_sol(i).x(:,k+1);
        u0_kp1 = S_sol(i).u(:,k+1);
        
        dx_kp1 = dx(:,k+1);
        du_kp1 = du(:,k+1);
    end
    
    r_k = x_k(P.i_r);
    v_k = x_k(P.i_v);
    T_k = u_k(P.i_T);
    
    % Initial conditions.
    if (k == 1)
        r_k == P.r_i;
        v_k == P.v_i;
        T_k == P.T_i;
    end
    
    % Final Conditions.
    if (k == P.K)
        r_k == P.r_f;
        v_k == P.v_f;
        T_k == P.T_f;
    end
    
    % Dynamics.
    if (k < P.K)
        x_kp1 == Ad{k}*x_k+Bd{k}*u_k+Cd{k}*u_kp1+zd{k}+nu_k;
    end
    
    % Non-convex keepout zones.
    [f_obs,A_obs,B_obs] = expr_obs(P,x0_k,u0_k);
    for j=1:P.n_obs
        f_obs{j}+A_obs{j}*dx_k+B_obs{j}*du_k <= eta_k(j);
    end
    
    % Lossless convexification.
    norm(T_k,2) <= s_k;
    
    % Thrust constraints.
    P.T_min <= s_k;
    s_k <= P.T_max;
    cos(P.theta_max)*s_k <= T_k(1);
    
    % Altitude constraint.
    r_k(1) == 0;
end

% End CVX.
cvx_end

% Store new solution.
S_sol(i+1).t = linspace(0,P.t_f,P.K);
S_sol(i+1).s = s;
S_sol(i+1).x = full(x);
S_sol(i+1).u = full(u);
S_sol(i+1).nu = full(nu);
S_sol(i+1).L = cvx_optval;
S_sol(i+1).J = compute_J(P,S_sol(i+1));
S_sol(i+1).r = NaN;
S_sol(i+1).Delta = norm(Delta,2);
S_sol(i+1).accepted = NaN;

% Print result summary.
fprintf('|')
fprintf(' i = %03d |',i);
% fprintf(' L = %07.4e |',S_sol(i+1).L)
% fprintf(' J = %07.4e |',S_sol(i+1).J)
fprintf(' nu = %07.2e |',norm(S_sol(i+1).nu(:),1))
% fprintf(' Del = %07.2e |',S_sol(i+1).Delta)
fprintf(' nIters = %03d |',cvx_slvitr)
fprintf(' %s |',cvx_status)
