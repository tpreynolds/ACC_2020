clear variables;
set(0,'defaulttextinterpreter','latex')

%% setup

S = struct;
C = struct;

% parameters
S.rho0 = 0.0;
S.rho1 = 0.1; 
S.rho2 = 0.9; 
S.alpha = 1.5; % shrink rate
S.beta  = 2.0; % growth rate
S.r_l   = 0.001; % minimum trust region size
S.w     = 1e3; % penalty weight
S.w_eta = 1e2;
S.w_eta_min = 1e2;

% options
S.type = 'ptr';
S.ineq = true;
S.plot = true;
S.quiet = false;
S.iter_max  = 50;
S.tol       = 1e-5;  % tolerance for dJ (convergence)
S.feas_tol  = 1e-3;  % tolerance to be considered feasible wrt constraints

% define constraints and cost
C.c = [0;1]; % cost function
C.a1 = [1; -2]; % a1 and a2 control the linear inequality
C.a2 = [-2; 2]; % a2 is on the left, a1 on the right
C.m = (C.a2(2)-C.a1(2))/(C.a2(1)-C.a1(1)); % slope of the line
C.f_cost = @(x,y)(C.c(1)*x+C.c(2)*y); % cost function
C.d_iq  = [ -C.m; 1 ]; % ineq. constraint is d' * x >= b;
C.b_iq  = C.a2(2) - C.m*C.a2(1); % y-intercept
C.f_iq  = @(x)(C.b_iq - dot(C.d_iq,x)); % ineq constraint LHS: f_iq(x) <= 0
C.c_eq  = [1;2;-1.2;-2;0];
C.f_eq  = @(x)( C.c_eq(1)*x.^4 ...
                + C.c_eq(2).*x.^3 ...
                + C.c_eq(3).*x.^2 ...
                + C.c_eq(4).*x ...
                + C.c_eq(5)); % eq constraint RHS: f_eq(x) == 0

% initial guess
S.x     = [-1.51;1.01];
S.J     = C.c'*S.x + S.w*penalty(S.x);
S.r     = NaN(S.iter_max,1);
S.dJ    = NaN(S.iter_max,1);
S.dL    = NaN(S.iter_max,1);
S.dx    = NaN(S.iter_max,1);
S.eta   = NaN(S.iter_max,1);
S.X_all = NaN(2,S.iter_max);
S.feas  = false;

% r_k = 0.2; eta = r_k;
% J_k = c'*x_k + w*penalty(x_k);
% bound states
S.x_min = [-2;-2];
S.x_max = [2;2];

%% solve

S = ptr_function(S,C);

%% plot
xs_1 = [-0.7372; 0.3163];
xs_2 = [0.6429; -1.0795];

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
colL = [0,152,0]./255;%[0.8510,0.3725,0.0078];

% Iterates + Cost
figure(8)
subplot(3,1,[1 2]), hold on, box on
plot(X1,C.f_eq(X1),'k','LineWidth',1)
if (S.ineq)
    area([C.a1(1) C.a2(1)],[C.a1(2) C.a2(2)],-3,'FaceColor','k','FaceAlpha',0.3)
end
for k = 1:S.cvrg_iter
    plot(S.X_all(1,k),S.X_all(2,k),'o','MarkerSize',5,...
        'MarkerFaceColor',N_col(:,k),'MarkerEdgeColor',N_col(:,k))
end
contour(X1,X2,C.f_cost(X,Y),-4:0.2:4,'k:')
plot(xs_1(1),xs_1(2),'b*','MarkerSize',10)
plot(xs_2(1),xs_2(2),'r*','MarkerSize',10)
set(gca,'Xlim',[-2,2],'Ylim',[-2,2]);
xlabel('$z_1$','FontSize',18)
ylabel('$z_2$','FontSize',18)
title('PTR Method','FontSize',16)
subplot(3,1,3), hold on, grid on, box on
x_data = 1:S.cvrg_iter;
plot(x_data,S.cost(x_data),'b','LineWidth',1,'HandleVisibility','off')
for k = 1:S.cvrg_iter
    if (S.feas(k))
        plot(x_data(k),S.cost(k),'ko',...
            'MarkerSize',5,'MarkerFaceColor','k','HandleVisibility','off')
    else
        plot(x_data(k),S.cost(k),'ro',...
            'MarkerSize',5,'MarkerFaceColor','r','HandleVisibility','off')
    end
end
ylabel('$c^{\top}z$','FontSize',18)
xlabel('Iteration Number','FontSize',16)

% iterates + TR
figure(9)
subplot(3,1,[1 2]), hold on, box on
plot(X1,C.f_eq(X1),'k','LineWidth',1)
if (S.ineq)
    area([C.a1(1) C.a2(1)],[C.a1(2) C.a2(2)],-3,'FaceColor','k','FaceAlpha',0.3)
end
for k = 1:S.cvrg_iter
    plot(S.X_all(1,k),S.X_all(2,k),'o','MarkerSize',5,...
        'MarkerFaceColor',N_col(:,k),'MarkerEdgeColor',N_col(:,k))
end
contour(X1,X2,C.f_cost(X,Y),-4:0.2:4,'k:')
set(gca,'Xlim',[-2,2],'Ylim',[-2,2]);
xlabel('$x_1$','FontSize',14)
ylabel('$x_2$','FontSize',14)
title('PTR Method','FontSize',16)
%
subplot(3,1,3), hold on, grid on, box on
x_data = 1:S.cvrg_iter;
yyaxis left
plot(x_data ,S.eta(x_data),'o','Color',colL,...
    'MarkerSize',5,'MarkerFaceColor',colL)
ylabel('$\eta$','FontSize',18)
y_ub = max(max(S.eta),max(S.dx));
for k = 1:S.cvrg_iter
   if (abs(S.dx(k)-S.eta(k))<S.feas_tol)
       plot(x_data(k),1.25*y_ub,'kx','MarkerSize',5)
   else
       plot(x_data(k),1.25*y_ub,'rx','MarkerSize',5)
   end
end
set(gca,'Ycolor',colL,'Ylim',[0 1.5*y_ub])
yyaxis right
plot(x_data ,S.dx(x_data),'*','Color',colR,...
        'MarkerSize',5,'MarkerFaceColor',colR)
ylabel('$\|x-\bar{x}\|$','FontSize',18)
set(gca,'Ycolor',colR,'Ylim',[0 1.5*y_ub])
xlabel('Iteration Number','FontSize',16)
% xticks(x_data)

%% other functions

% penalty function
function val = penalty(x)
    val = norm(x,1);
end
