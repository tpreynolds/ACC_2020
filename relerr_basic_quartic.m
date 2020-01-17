clear variables;
set(0,'defaulttextinterpreter','latex')

%% setup

S = struct;
C = struct;
% SCvx parameters
S.rho0 = 0.1;
S.rho1 = 0.9; 
S.rho2 = 1.0; 
S.alpha = 1.5; % shrink rate
S.beta = 2.0; % growth rate
S.r_l = 0.005; % minimum trust region size
S.w = 4e2;%2e2; % penalty weight
S.type = '2';
S.ineq = true;
S.plot = true;
S.quiet = false;
S.iter_max  = 120;
S.tol       = 1e-3; % tolerance for dJ_k (convergence)
S.feas_tol  = 1e-3;

S.X_all = zeros(2,S.iter_max);

% constraints & cost
C.c = [1;1]; % cost function
C.a1 = [1; -2]; % a1 and a2 control the linear inequality
C.a2 = [-2; 2]; % a2 is on the left, a1 on the right
m = (C.a2(2)-C.a1(2))/(C.a2(1)-C.a1(1)); % slope

C.f_cost    = @(x,y)(C.c(1)*x + C.c(2)*y); % cost function
C.d_iq    = [-m;1]; % ineq. constraint is d' * x >= b;
C.b_iq    = C.a2(2) - m*C.a2(1); % y-intercept
C.f_iq    = @(x)(C.b_iq - dot(C.d_iq,x));
C.c_eq    = [1;2;-1.2;-2;0];
C.f_eq    = @(x)(C.c_eq(1)*x.^4 ...
                    + C.c_eq(2).*x.^3 ...
                    + C.c_eq(3).*x.^2 ...
                    + C.c_eq(4).*x ...
                    + C.c_eq(5)); % eq. constraint RHS

% initial guess
S.x = [1.5;1.5]; % [0.5,1.5]
S.J = C.c'*S.x + S.w*penalty(S.x);
S.r = NaN(S.iter_max,1);
S.dJ = NaN(S.iter_max,1);
S.dL = NaN(S.iter_max,1);
S.dx = NaN(S.iter_max,1);
S.rho = NaN(S.iter_max,1);
S.reject = NaN(S.iter_max,1);
% initial trust region size
S.r(1) = 0.1;
S.reject(1) = false; % obviously accept the first iter

% bound states
S.x_min = [-2;-2];
S.x_max = [2;2];

%% solve
S = relerr_function(S,C);
S.cvrg_iter = 65;
%% plot
set(0,'defaultAxesTickLabelInterpreter','latex','defaultLegendInterpreter','latex')
xs_1 = [-0.7372; 0.3163];
xs_2 = [0.5288; -1.0192];
    
% final plot
ic_col = [0,0.6,1];%./255;
fc_col = [1,0,1];%./255;
for k = 1:3
    N_col(k,:) = linspace(ic_col(k),fc_col(k),S.cvrg_iter);
end
X1 = -3:0.01:2;
X2 = -3:0.01:2;
[X,Y] = meshgrid(X1);

colR = [76,0,153]./255;
colL = [0,130,0]./255;%[0.8510,0.3725,0.0078];

% Iterates + Cost
figure(8), clf, hold on, box on
h1=plot(X1,C.f_eq(X1),'k','LineWidth',1);
if (S.ineq)
    h2=area([C.a1(1) C.a2(1)],[C.a1(2) C.a2(2)],-3,'FaceColor','k','FaceAlpha',0.3);
end
for k = 1:S.cvrg_iter
    plot(S.X_all(1,k),S.X_all(2,k),'o','MarkerSize',5,...
        'MarkerFaceColor',N_col(:,k),'MarkerEdgeColor',N_col(:,k),'HandleVisibility','off');
end
h3=contour(X1,X2,C.f_cost(X,Y),-4:0.2:4,'k:');
h4=plot(xs_1(1),xs_1(2),'b*','MarkerSize',10);
h5=plot(xs_2(1),xs_2(2),'r*','MarkerSize',10);
set(gca,'FontSize',12)
legend('equality constraint','inequality constraint','cost lines',...
    'local min','global min','Location','northwest','FontSize',12)
colormap(N_col')
hC = colorbar('TickLabelInterpreter','latex','FontSize',14);
hC.Label.String = 'Iteration';
hC.Label.Interpreter = 'latex';
hC.Label.FontSize = 14;
caxis([0 min(S.iter_max,S.cvrg_iter)]);
set(gca,'Xlim',[-2,2],'Ylim',[-2,2]);
xlabel('$z_1$','FontSize',16)
ylabel('$z_2$','FontSize',16)
% title('SCvx Method')
% subplot(3,1,3), hold on, grid on, box on
% x_data = 1:S.cvrg_iter;
% plot(x_data,S.cost(x_data),'b','LineWidth',1,'HandleVisibility','off')
% for k = 1:S.cvrg_iter
%     if (S.feas(k))
%         plot(x_data(k),S.cost(k),'ko',...
%             'MarkerSize',5,'MarkerFaceColor','k','HandleVisibility','off')
%     else
%         plot(x_data(k),S.cost(k),'ro',...
%             'MarkerSize',5,'MarkerFaceColor','r','HandleVisibility','off')
%     end
% end
% ylabel('$c^{\top}z$','FontSize',18)
% xlabel('Iteration Number','FontSize',16)


gcol = 0.1*[0,1,0]+0.9*[1,1,1];
ycol = 0.1*[1,1,0]+0.9*[1,1,1];
rcol = 0.1*[1,0,0]+0.9*[1,1,1];

figure(9)%, clf
subplot(2,1,2), cla, hold on, grid on, box on
F = fill([0 S.iter_max S.iter_max 0],[S.rho2 S.rho2 2 2],...
        rcol,'LineStyle','none','FaceAlpha',0.8);
fill([0 S.iter_max S.iter_max 0],[S.rho1 S.rho1 S.rho2 S.rho2],...
        ycol,'LineStyle','none','FaceAlpha',0.8); 
fill([0 S.iter_max S.iter_max 0],[-1 -1 S.rho1 S.rho1],...
        gcol,'LineStyle','none','FaceAlpha',0.8); 
plot(1:S.iter_max,S.rho,'k-o','MarkerSize',4,'MarkerFaceColor','k','LineWidth',1);
plot([1 S.iter_max],[S.rho0 S.rho0],'r--','LineWidth',1)
plot([1 S.iter_max],[S.rho1 S.rho1],'r--','LineWidth',1)
plot([1 S.iter_max],[S.rho2 S.rho2],'r--','LineWidth',1)
set(gca,'Ylim',[-0.25 1.25],'FontSize',14)
set(gca,'Xlim',[0 min(S.iter_max,S.cvrg_iter)])
xlabel('Iteration','FontSize',14)
ylabel('$\rho$','FontSize',16)
title('Relative Error','FontSize',15)

% Iterates + TR
% figure(9)
% subplot(3,1,[1 2]), hold on, box on
% plot(X1,C.f_eq(X1),'k','LineWidth',1)
% if (S.ineq)
%     area([C.a1(1) C.a2(1)],[C.a1(2) C.a2(2)],-3,'FaceColor','k','FaceAlpha',0.3)
% end
% for k = 1:S.cvrg_iter
%     plot(S.X_all(1,k),S.X_all(2,k),'o','MarkerSize',5,...
%         'MarkerFaceColor',N_col(:,k),'MarkerEdgeColor',N_col(:,k))
% end
% contour(X1,X2,C.f_cost(X,Y),-4:0.2:4,'k:')
% % plot(xs_1(1),xs_1(2),'b*','MarkerSize',10)
% % plot(xs_2(1),xs_2(2),'r*','MarkerSize',10)
% set(gca,'Xlim',[-2,2],'Ylim',[-2,2]);
% xlabel('$z_1$','FontSize',18)
% ylabel('$z_2$','FontSize',18)
% title('SCvx Method','FontSize',16)
% %
% subplot(3,1,3), hold on, grid on, box on
% x_data = 1:S.cvrg_iter+1;
% % yyaxis left
% plot(x_data ,S.r(x_data),'o','Color',colL,...
%     'MarkerSize',5,'MarkerFaceColor',colL,'HandleVisibility','off')
% ylabel('$\eta$','FontSize',18)
% y_ub = max(max(S.r),max(S.dx));
% for k = 1:S.cvrg_iter+1
%    if (~S.reject(k))
%        plot(x_data(k),1.25*y_ub,'kx','MarkerSize',9,'HandleVisibility','off')
%    else
%        plot(x_data(k),1.25*y_ub,'rx','MarkerSize',9,'HandleVisibility','off')
%    end
% end
% set(gca,'Ycolor',colL,'Ylim',[0 1.5*y_ub])
% % yyaxis right
% % plot(x_data ,S.dx(x_data),'*','Color',colR,...
% %         'MarkerSize',5,'MarkerFaceColor',colR)
% % ylabel('$\|z-\bar{z}\|$','FontSize',18)
% % set(gca,'Ycolor',colR,'Ylim',[0 1.5*y_ub])
% xlabel('Iteration Number','FontSize',16)

% figure, hold on, grid on, box on
% plot(x_data,S.cost(x_data),'b','LineWidth',1,'HandleVisibility','off')
% for k = 1:S.cvrg_iter
%     if (S.feas(k))
%         plot(x_data(k),S.cost(k),'ko',...
%             'MarkerSize',5,'MarkerFaceColor','k','HandleVisibility','off')
%     else
%         plot(x_data(k),S.cost(k),'ro',...
%             'MarkerSize',5,'MarkerFaceColor','r','HandleVisibility','off')
%     end
% end
% ylabel('Original Cost','FontSize',16)
% xlabel('Iteration Number','FontSize',16)

% penalty function
function val = penalty(x)
    val = norm(x,1);
end