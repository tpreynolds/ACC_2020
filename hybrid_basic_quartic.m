clear variables;
set(0,'defaulttextinterpreter','latex',...
      'defaultAxesTickLabelInterpreter','latex',...
      'defaultLegendInterpreter','latex','defaultAxesFontSize',12)

%% setup

S = struct;
C = struct;

% SCP parameters
S.tr_shrink = 1.5;      % shrink rate
S.tr_grow   = 2.0;      % growth rate
S.r_l       = 0.005;    % minimum trust region size
S.w         = 4e2;      % penalty weight
S.type      = '2';
S.ineq      = true;
S.plot      = true;
S.scp_algo  = 'rel_err';
if strcmp(S.scp_algo,'rel_dec')
    S.rho0      = 0.0;
    S.rho1      = 0.1;
    S.rho2      = 0.9;
    S.rho_ub    = S.rho2 - 0.001;
    S.rho_lb    = S.rho1 + 0.001;
else
    S.rho0      = 0.1;
    S.rho1      = 0.9;
    S.rho2      = 1.0;
    S.rho_ub    = S.rho1 - 0.001;
    S.rho_lb    = S.rho0 + 0.001;
end
S.prev_event    = false;
% FSQP parameters
S.alpha = 0.1;
S.beta  = 0.5;
S.delta = 2.0;
S.gamma = 1.0;
S.eta   = 0.1;
S.kappa = 2.1;
S.tau   = 2.5;
S.M     = 10.0;
S.wsqp  = 5e0;          % FSQP penalty weight
S.cvrg_iter = NaN;
S.iter_max  = 60;       % global iterations
S.tol       = 1e-2;     % tolerance (convergence)
S.feas_tol  = 2e-3;     % tolerance (feasibility)

% constraints & cost
C.c         = [1;1];                        % cost function
C.f_cost    = @(x,y)(C.c(1)*x + C.c(2)*y);  % cost function
C.lb        = [-2;-2];                      % state lower bound
C.ub        = [2;2];                        % state upper bound
C.a1        = [1; -2];          % a1 and a2 control the linear inequality
C.a2        = [-2; 2];          % a2 is on the left, a1 on the right
C.m         = (C.a2(2)-C.a1(2))/(C.a2(1)-C.a1(1)); % slope
C.d_iq      = [-C.m;1];                     % ineq. constraint is b-d'*x<=0;
C.b_iq      = C.a2(2) - C.m*C.a2(1);        % y-intercept
C.f_iq      = @(x)(C.b_iq - dot(C.d_iq,x)); % ineq. constraint RHS
C.ineq      = S.ineq;
C.c_eq      = [1;2;-1.2;-2;0];              % eq. constraint coefficients
C.f_eq      = @(x)(C.c_eq(1)*x.^4 + C.c_eq(2).*x.^3 + C.c_eq(3).*x.^2 + ...
                    C.c_eq(4).*x + C.c_eq(5));      % eq. constraint RHS
C.rts       = roots([4*3*C.c_eq(1) 3*2*C.c_eq(2) 2*1*C.c_eq(3)]);
C.xs_1      = [-0.7372; 0.3163];   % local minimum
C.xs_2      = [0.5288; -1.0192];   % global minimum

% initial guess and data matrices
S.x     = [1.5;1.5];
S.dx    = NaN(S.iter_max,1);
S.X_all = [S.x,NaN(2,S.iter_max-1)];
S.r     = NaN(S.iter_max,1);
S.rho   = NaN(S.iter_max,1);
S.dJ    = NaN(S.iter_max,1);
S.dL    = NaN(S.iter_max,1);
S.rej   = NaN(S.iter_max,1);
% initial trust region size
S.r(1)   = 0.1;
S.J      = C.c'*S.x + S.w*penalty(S.x);
S.rej(1) = false; % obviously accept the first iter

%% Solve
S = hybrid_function(S,C);

%% Plot

% final plot
ic_col = [0,0.6,1];%./255;
fc_col = [1,0,1];%./255;
lim     = min(S.iter_max,S.cvrg_iter);
N_col   = zeros(3,lim);
for k = 1:3
    N_col(k,:) = linspace(ic_col(k),fc_col(k),lim);
end
X1 = -3:0.01:2;
X2 = -3:0.01:2;
[X,Y] = meshgrid(X1);

% Iterates + Cost
figure(8), clf, hold on, box on
h1=plot(X1,C.f_eq(X1),'k','LineWidth',1);
if (S.ineq)
    h2=area([C.a1(1) C.a2(1)],[C.a1(2) C.a2(2)],-3,'FaceColor','k','FaceAlpha',0.3);
end
for k = 1:lim
    plot(S.X_all(1,k),S.X_all(2,k),'o','MarkerSize',5,...
        'MarkerFaceColor',N_col(:,k),'MarkerEdgeColor',N_col(:,k),'HandleVisibility','off');
end
h3=contour(X1,X2,C.f_cost(X,Y),-4:0.2:4,'k:');
h4=plot(C.xs_1(1),C.xs_1(2),'b*','MarkerSize',10);
h5=plot(C.xs_2(1),C.xs_2(2),'r*','MarkerSize',10);
set(gca,'FontSize',14)
legend('equality constraint','inequality constraint','cost lines',...
    'local min','global min','Location','northwest','FontSize',14)
colormap(N_col')
hC = colorbar('TickLabelInterpreter','latex','FontSize',15);
hC.Label.String = 'Iteration';
hC.Label.Interpreter = 'latex';
hC.Label.FontSize = 16;
% hC.Ticks = [hC.Ticks,lim];
caxis([0 lim]);
set(gca,'Xlim',[-2,2],'Ylim',[-2,2]);
xlabel('$z_1$','FontSize',18)
ylabel('$z_2$','FontSize',18)
