function [f,A,B] = expr_obs(P,x,~)

H_r = zeros(2,P.n_x);
H_r(:,2:3) = eye(2);

f = cell(P.n_obs,1);
A = cell(P.n_obs,1);
B = cell(P.n_obs,1);

for j=1:P.n_obs
  f{j} = P.R_obs(j)-norm(H_r*x-P.r_obs(:,j),2);
  A{j} = -(H_r*x-P.r_obs(:,j))'/norm(H_r*x-P.r_obs(:,j),2)*H_r;
  B{j} = zeros(1,P.n_u);
end
