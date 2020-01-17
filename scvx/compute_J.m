function [J] = compute_J(P,S_sol)

s = S_sol.s;
x = S_sol.x;
u = S_sol.u;

nu = zeros(P.n_x,P.K-1);
[~,Ad,Bd,Cd,zd] = propagate(P,S_sol,'step');

for k=1:P.K-1
	nu(:,k) = Ad{k}*x(:,k)+Bd{k}*u(:,k)+Cd{k}*u(:,k+1)+zd{k}-x(:,k+1);
end

eta = zeros(P.n_obs,P.K);
for k=1:P.K
  [f_obs,~,~] = expr_obs(P,x(:,k),u(:,k));
  for j=1:P.n_obs
    eta(j,k) = max(0,f_obs{j});
  end
end

J = sum(s)+P.SCvx.lambda*norm(nu(:),1)+P.SCvx.lambda*sum(eta(:));
