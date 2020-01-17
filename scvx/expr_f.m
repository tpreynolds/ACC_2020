function [f,A,B,z] = expr_f(P,x,u)

n_x = P.n_x;
n_u = P.n_u;
i_r = P.i_r;
i_v = P.i_v;

I_3x3 = eye(3);
v = x(i_v);
a = (1/P.m)*u+P.g-P.C_d*norm(v,2)*v;

f = zeros(n_x,1);
f(i_r) = v;
f(i_v) = a;

A = zeros(n_x,n_x);
A(i_r,i_v) = I_3x3;
A(i_v,i_v) = zeros(3,3)-P.C_d*((v*v')/norm(v,2)+norm(v,2)*I_3x3);

B = zeros(n_x,n_u);
B(i_v,:) = (1/P.m)*I_3x3;

z = zeros(n_x,1);
z(:) = f-A*x-B*u;
