clear all %#ok<CLALL>
clc

P = inputs();

[P,S_sol,S_den] = SCvx(P);

set(0,'defaulttextinterpreter','latex')
plot_all(P,S_sol,S_den)

save('SCvx.mat','P','S_den','S_sol')
