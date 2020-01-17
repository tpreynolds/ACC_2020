function [P,S_sol,S_den,Ad,Bd,Cd,zd,varargout] = initialize(P)

S_sol = struct;
S_sol.t = linspace(0,P.t_f,P.K);
S_sol.s = zeros(1,P.K);
S_sol.x = zeros(P.n_x,P.K);
S_sol.u = zeros(P.n_u,P.K);
S_sol.nu = zeros(P.n_x,P.K-1);
S_sol.eta = zeros(P.n_obs,P.K);
S_sol.L = NaN;
S_sol.J = NaN;
S_sol.r = NaN;
S_sol.Delta = NaN;
S_sol.accepted = NaN;

for i=1:3
  S_sol.x(P.i_r(i),:) = linspace(P.r_i(i),P.r_f(i),P.K);
  S_sol.x(P.i_v(i),:) = linspace(P.v_i(i),P.v_f(i),P.K);
end
S_sol.u(P.i_T(1),:) = P.m*norm(P.g);

[S_den,Ad,Bd,Cd,zd] = propagate(P,S_sol,'step');
S_sol.J = compute_J(P,S_sol);

varargout{1} = P.SCvx.Delta_i;
varargout{2} = 1;
