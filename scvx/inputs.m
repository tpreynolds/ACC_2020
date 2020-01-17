function [P] = inputs()

P = struct;

P.DEG2RAD = pi/180;
P.RAD2DEG = 180/pi;

P.cvx_solver = 'sdpt3';
P.cvx_precision = 'low';
P.cvx_quiet = true;

P.i_r = (1:3);
P.i_v = (1:3)+P.i_r(end);
P.i_T = (1:3);
P.n_x = P.i_v(end);
P.n_u = P.i_T(end);

P.K = 31;
% P.K = 43;
P.K_sub = 15;

P.SCvx.max_iter = 50;
P.SCvx.Delta_i = 1;
P.SCvx.lambda = 1e5;
P.SCvx.alpha = 2.0;
P.SCvx.beta = 3.2;
% P.SCvx.dL_tol = 1e-3;
P.SCvx.dJ_tol = 1e-3;
P.SCvx.rho_0 = 0.00;
P.SCvx.rho_1 = 0.2;
% P.SCvx.rho_2 = 0.9;
P.SCvx.rho_2 = 0.7;

P.t_f = 3;
% P.t_f = 4.2;
P.g = [-9.81 0 0]';
P.theta_max = 45*P.DEG2RAD;
P.T_min = 1;
P.T_max = 4;
P.m = 0.3;
P.C_d = 0.5;

P.R_obs = [1 1];
P.r_obs = [3 0.45; 7 -0.45]';
% P.R_obs = [1 1 1];
% P.r_obs = [3 0.45; 7 -0.45; 10 0.3]';
P.n_obs = length(P.R_obs);

P.r_i = [0 0 0]';
P.v_i = [0 0.5 0]';
P.T_i = -P.m*P.g;

P.r_f = [0 10 0]';
% P.r_f = [0 14 0]';
P.v_f = [0 0.5 0]';
P.T_f = -P.m*P.g;
