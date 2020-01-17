function [S_den,varargout] = propagate(P,S_sol,mode)

if isequal(mode,'step')
  I_x = 1:P.n_x;
  I_Ad = (1:P.n_x^2)+I_x(end);
  I_Bd = (1:P.n_x*P.n_u)+I_Ad(end);
  I_Cd = (1:P.n_x*P.n_u)+I_Bd(end);
  I_zd = (1:P.n_x)+I_Cd(end);
  n_X = P.n_x+P.n_x^2+2*P.n_x*P.n_u+P.n_x;

  Ad = cell(P.K-1,1);
  Bd = cell(P.K-1,1);
  Cd = cell(P.K-1,1);
  zd = cell(P.K-1,1);

  S_den = struct;
  S_den.t = cell(P.K-1,1);
  S_den.x = cell(P.K-1,1);
  S_den.u = cell(P.K-1,1);

  for k=1:P.K-1
    x_i = S_sol.x(:,k);
    Ad_i = eye(P.n_x,P.n_x);
    Bd_bar_i = zeros(P.n_x,P.n_u);
    Cd_bar_i = zeros(P.n_x,P.n_u);
    zd_bar_i = zeros(P.n_x,1);

    t_k = S_sol.t(k);
    t_kp1 = S_sol.t(k+1);
    t = linspace(t_k,t_kp1,P.K_sub)';

    X_i = zeros(n_X,1);
    X_i(I_x) = x_i;
    X_i(I_Ad) = reshape(Ad_i,P.n_x^2,1);
    X_i(I_Bd) = reshape(Bd_bar_i,P.n_x*P.n_u,1);
    X_i(I_Cd) = reshape(Cd_bar_i,P.n_x*P.n_u,1);
    X_i(I_zd) = zd_bar_i;

    [~,X] = ode45(@(t,X)derivs_step(P,S_sol,t_k,t_kp1,X,t),t,X_i);

    S_den.t{k} = t;
    S_den.x{k} = X(:,I_x)';
    S_den.u{k} = interp1(S_sol.t,S_sol.u',t,'linear')';

    Ad_k = reshape(X(end,I_Ad),P.n_x,P.n_x);
    Bd_bar_k = reshape(X(end,I_Bd),P.n_x,P.n_u);
    Cd_bar_k = reshape(X(end,I_Cd),P.n_x,P.n_u);
    zd_bar_k = reshape(X(end,I_zd),P.n_x,1);

    Ad{k} = Ad_k;
    Bd{k} = Ad_k*Bd_bar_k;
    Cd{k} = Ad_k*Cd_bar_k;
    zd{k} = Ad_k*zd_bar_k;
  end

  varargout = cell(4,1);
  varargout{1} = Ad;
  varargout{2} = Bd;
  varargout{3} = Cd;
  varargout{4} = zd;
elseif isequal(mode,'full')
  for k=1:P.K-1
    if (k == 1)
      x_i = S_sol.x(:,k);
    else
      x_i = S_den.x{k-1}(:,end);
    end
    
    t_k = S_sol.t(k);
    t_kp1 = S_sol.t(k+1);
    t = linspace(t_k,t_kp1,P.K_sub)';
    
    [~,x] = ode45(@(t,x)derivs_full(P,S_sol,x,t),t,x_i);
    
    S_den.t{k} = t;
    S_den.x{k} = x';
    S_den.u{k} = interp1(S_sol.t,S_sol.u',t,'linear')';
	end
end

function [x_dot] = derivs_full(P,S_sol,x,t)
u_t = interp1(S_sol.t,S_sol.u',t,'linear')';
[f_t,~,~] = expr_f(P,x,u_t);
x_dot = f_t;

function [X_dot] = derivs_step(P,S_sol,t_k,t_kp1,X,t)
I_x = 1:P.n_x;
I_Ad = (1:P.n_x^2)+I_x(end);
I_Bd = (1:P.n_x*P.n_u)+I_Ad(end);
I_Cd = (1:P.n_x*P.n_u)+I_Bd(end);
I_zd = (1:P.n_x)+I_Cd(end);
n_X = P.n_x+P.n_x^2+2*P.n_x*P.n_u+P.n_x;

x_t = X(I_x);
u_t = interp1(S_sol.t,S_sol.u',t,'linear')';

[f_t,A_t,B_t,z_t] = expr_f(P,x_t,u_t);
Ad_t = reshape(X(I_Ad),P.n_x,P.n_x);

x_dot = f_t;
Ad_dot = A_t*Ad_t;
Bd_bar_dot = (t_kp1-t)/(t_kp1-t_k)*(Ad_t\B_t);
Cd_bar_dot = (t-t_k)/(t_kp1-t_k)*(Ad_t\B_t);
zd_bar_dot = Ad_t\z_t;

X_dot = zeros(n_X,1);
X_dot(I_x) = x_dot;
X_dot(I_Ad) = reshape(Ad_dot,P.n_x^2,1);
X_dot(I_Bd) = reshape(Bd_bar_dot,P.n_x*P.n_u,1);
X_dot(I_Cd) = reshape(Cd_bar_dot,P.n_x*P.n_u,1);
X_dot(I_zd) = zd_bar_dot;
