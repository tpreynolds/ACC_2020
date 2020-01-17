function [] = plot_all(P,S_sol,~)

fs = 20;
ms = 10;
t = S_sol(end).t;
rU = S_sol(end).x(P.i_r(1),:);
rE = S_sol(end).x(P.i_r(2),:);
rN = S_sol(end).x(P.i_r(3),:);
vU = S_sol(end).x(P.i_v(1),:);
vE = S_sol(end).x(P.i_v(2),:);
vN = S_sol(end).x(P.i_v(3),:);
TU = S_sol(end).u(P.i_T(1),:);
TE = S_sol(end).u(P.i_T(2),:);
TN = S_sol(end).u(P.i_T(3),:);

figure(1),clf(1)

subplot(3,1,1),hold on,grid on
plot(t,rU,'color','r','marker','.','markersize',ms)
plot(t,rE,'color','g','marker','.','markersize',ms)
plot(t,rN,'color','b','marker','.','markersize',ms)
set(gca,'fontsize',fs)
xlabel('Time [s]','fontsize',fs)
ylabel('Position [m]','fontsize',fs)
legend({'Up','East','North'},'fontsize',fs)

subplot(3,1,2),hold on,grid on
plot(t,vU,'color','r','marker','.','markersize',ms)
plot(t,vE,'color','g','marker','.','markersize',ms)
plot(t,vN,'color','b','marker','.','markersize',ms)
set(gca,'fontsize',fs)
xlabel('Time [s]','fontsize',fs)
ylabel('Velocity [m/s]','fontsize',fs)
legend({'Up','East','North'},'fontsize',fs)

subplot(3,1,3),hold on,grid on
plot(t,TU,'color','r','marker','.','markersize',ms)
plot(t,TE,'color','g','marker','.','markersize',ms)
plot(t,TN,'color','b','marker','.','markersize',ms)
set(gca,'fontsize',fs)
xlabel('Time [s]','fontsize',fs)
ylabel('Thrust [N]','fontsize',fs)
legend({'Up','East','North'},'fontsize',fs)

figure(2),clf(2)

T_mag = sum(TU.^2+TE.^2+TN.^2,1).^0.5;
T_tilt = acosd(TU./T_mag);

subplot(3,1,1),hold on,grid on
plot3(rE,rN,rU,'color','r','marker','.','markersize',ms)
for k=1:P.K
  fT = 0.5;
  plot3(rE(k)+fT*[0 TE(k)], ...
        rN(k)+fT*[0 TN(k)], ...
        rU(k)+fT*[0 TU(k)],'color','b','marker','none')
end
for i=1:P.n_obs
  thetas = linspace(0,360,100);
  cir_E = P.R_obs(i)*cosd(thetas)+P.r_obs(1,i);
  cir_N = P.R_obs(i)*sind(thetas)+P.r_obs(2,i);
  plot3(cir_E,cir_N,zeros(size(cir_E)),'color','k')
end
set(gca,'fontsize',fs)
xlabel('East [m]','fontsize',fs)
ylabel('North [m]','fontsize',fs)
zlabel('Up [m]','fontsize',fs)
axis equal

subplot(3,1,2),hold on,grid on
plot(t,T_tilt,'color','r','marker','.','markersize',ms)
plot(t,P.RAD2DEG*P.theta_max*ones(size(t)),'color','k')
set(gca,'fontsize',fs)
xlabel('Time [s]','fontsize',fs)
ylabel('Thrust Tilt [deg]','fontsize',fs)

subplot(3,1,3),hold on,grid on
plot(t,T_mag,'color','r','marker','.','markersize',ms)
plot(t,P.T_min*ones(size(t)),'color','k')
plot(t,P.T_max*ones(size(t)),'color','k')
set(gca,'fontsize',fs)
xlabel('Time [s]','fontsize',fs)
ylabel('Thrust Magnitude [N]','fontsize',fs)

X_star = [S_sol(end).x(:); S_sol(end).u(:)];
delta = zeros(length(S_sol),1);

for i=1:length(delta)
  X_i = [S_sol(i).x(:); S_sol(end).u(:)];
  delta(i) = norm(X_i-X_star,2);
end

figure(3),hold on,grid on
set(gca,'yscale','log')
plot(0:length(S_sol)-2,delta(1:end-1),'color','b','marker','.','markersize',10)
% xlabel('SCvx Iteration')
set(gca,'fontsize',fs)
xlabel('Iteration','fontsize',fs)
ylabel('$\|X^k-X^*\|_2$','fontsize',fs)
legend({'SCvx'},'fontsize',fs,'location','best')

saveas(1,'SCvx_Pos-Vel-Thrust.png')
saveas(2,'SCvx_Traj-Tilt-Thrust.png')
saveas(3,'SCvx_Convergence.png')
