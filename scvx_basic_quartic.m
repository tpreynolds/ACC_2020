clear variables;
set(0,'defaulttextinterpreter','latex')

%% setup

S = struct;
C = struct;
% SCvx parameters
S.rho0 = 0.0;
S.rho1 = 0.1; 
S.rho2 = 0.9; 
S.alpha = 1.5; % shrink rate
S.beta = 2.0; % growth rate
S.r_l = 0.001; % minimum trust region size
S.w = 1e3; % penalty weight
S.type = 'inf';
S.ineq = true;
S.plot = true;
S.quiet = false;
S.iter_max = 50;
S.tol = 1e-5; % tolerance for dJ_k (convergence)

S.X_all = zeros(2,S.iter_max);

% constraints & cost
C.c = [0;1]; % cost function
C.a1 = [1; -2]; % a1 and a2 control the linear inequality
C.a2 = [-2; 2]; % a2 is on the left, a1 on the right
m = (C.a2(2)-C.a1(2))/(C.a2(1)-C.a1(1)); % slope

C.cost    = @(x,y)(C.c(1)*x + C.c(2)*y); % cost function
C.d_iq    = [-m;1]; % ineq. constraint is d' * x >= b;
C.b_iq    = C.a2(2) - m*C.a2(1); % y-intercept
C.c_eq    = [1;2;-1.2;-2;0];
C.f_eq    = @(x)(C.c_eq(1)*x.^4 ...
                    + C.c_eq(2).*x.^3 ...
                    + C.c_eq(3).*x.^2 ...
                    + C.c_eq(4).*x ...
                    + C.c_eq(5)); % eq. constraint RHS

% initial guess
S.x = [1.5;1.5];
S.J = C.c'*S.x + S.w*penalty(S.x);
S.r = NaN(S.iter_max,1);
S.dJ = NaN(S.iter_max,1);
S.dL = NaN(S.iter_max,1);
S.dx = NaN(S.iter_max,1);
S.reject = NaN(S.iter_max,1);
% initial trust region size
S.r(1) = 0.1;
S.reject(1) = false; % obviously accept the first iter

% bound states
S.x_min = [-2;-2];
S.x_max = [2;2];

%% solve
S = scvx_function(S,C);

%% plot

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
contour(X1,X2,C.cost(X,Y),-4:0.2:4,'k:')
set(gca,'Xlim',[-2,2],'Ylim',[-2,2]);
xlabel('$x_1$','FontSize',14)
ylabel('$x_2$','FontSize',14)
title('SCvx Method','FontSize',16)
%
subplot(3,1,3), hold on, grid on, box on
x_data = 1:S.cvrg_iter;
yyaxis left
plot(x_data ,S.r(x_data),'o','Color',colL,...
    'MarkerSize',5,'MarkerFaceColor',colL)
ylabel('$\rho$','FontSize',18)
for k = 1:S.cvrg_iter
   if (S.reject(k))
       plot(x_data(k),3,'rx','MarkerSize',5)
   else
       plot(x_data(k),3,'kx','MarkerSize',5)
   end
end
set(gca,'Ycolor',colL,'Ylim',[0 4])
yyaxis right
plot(x_data ,S.dx(x_data),'*','Color',colR,...
        'MarkerSize',5,'MarkerFaceColor',colR)
ylabel('$\|x-\bar{x}\|$','FontSize',18)
set(gca,'Ycolor',colR,'Ylim',[0 4])
xlabel('Iteration Number','FontSize',16)
% xticks(x_data)

% penalty function
function val = penalty(x)
    val = norm(x,1);
end