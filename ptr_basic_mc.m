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
S.type = 'tpr';
S.ineq = true;
S.plot = false;
S.quiet = true;
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
S.x     = [-1.5;1.5];
S.J     = C.c'*S.x + S.w*penalty(S.x);
S.r     = NaN(S.iter_max,1);
S.dJ    = NaN(S.iter_max,1);
S.dL    = NaN(S.iter_max,1);
S.dx    = NaN(S.iter_max,1);
S.eta   = NaN(S.iter_max,1);
S.X_all = NaN(2,S.iter_max);
S.feas  = false;

% bound states
S.x_min = [-2;-2];
S.x_max = [2;2];

% initial guess matrix
X_ic = S.x_min(1):0.25:S.x_max(1);
Y_ic = S.x_min(2):0.25:S.x_max(2);

%% solve

MC = struct;
MC.N_x  = numel(X_ic);
MC.N_y  = numel(Y_ic);
MC.N    = MC.N_x * MC.N_y;
MC.ic   = cell(MC.N_x,MC.N_y);
MC.xf   = cell(MC.N_x,MC.N_y);
MC.cvrg_iter = cell(MC.N_x,MC.N_y);
MC.tr   = cell(MC.N_x,MC.N_y);
MC.dx   = cell(MC.N_x,MC.N_y);
MC.feas = cell(MC.N_x,MC.N_y);

fprintf('Running Monte Carlo...\n')
for i = 1:MC.N_x
    fprintf('X_ic(%d/%d) and {',i,MC.N_x)
    for j = 1:MC.N_y
        if (j<MC.N_y)
            fprintf('Y_ic(%d), ',j)
        else
            fprintf('Y_ic(%d)}\n',j)
        end
        S.x     = [X_ic(i);Y_ic(j)];
        S.J     = C.c'*S.x + S.w*penalty(S.x);
        S.r     = NaN(S.iter_max,1);
        S.dJ    = NaN(S.iter_max,1);
        S.dL    = NaN(S.iter_max,1);
        S.dx    = NaN(S.iter_max,1);
        % initial trust region size
        S.r(1)  = 0.2;
        % solve
        S   = ptr_function(S,C);
        % get data
        MC.ic{i,j}          = S.X_all(:,1);
        MC.xf{i,j}          = S.x;
        MC.cvrg_iter{i,j}   = S.cvrg_iter;
        MC.tr{i,j}          = S.r;
        MC.dx{i,j}          = S.dx;
        MC.feas{i,j}        = S.feas;
    end
end

%% post process
xs_1 = [-0.7372; 0.3163];
id_xs_1 = 1;
xs_2 = [0.6429; -1.0795];
id_xs_2 = 2;
id_infeas = -1;
id_other = 3;

% check to see if we got to either local optima
for i = 1:MC.N_x
    for j = 1:MC.N_y
        x_final = MC.xf{i,j};
        if (norm(xs_1-x_final)<1e-2)
            MC.result{i,j} = id_xs_1;
        elseif (norm(xs_2-x_final)<1e-2)
            MC.result{i,j} = id_xs_2;
        else
            if (MC.feas{i,j})
                MC.result{i,j} = id_other;
            else
                MC.result{i,j} = id_infeas;
            end
        end
    end
end

% save workspace to data folder
save(strcat('data/mc_ptr_',S.type));

%% plot

% final plot
% ic_col = [0,0.6,1];%./255;
% fc_col = [1,0,1];%./255;
% for k = 1:3
%     N_col(k,:) = linspace(ic_col(k),fc_col(k),S.cvrg_iter);
% end
X1 = -3:0.01:2;
X2 = -3:0.01:2;
[X,Y] = meshgrid(X1);

figure(1), hold on, box on
plot(X1,C.f_eq(X1),'k','LineWidth',1)
if (S.ineq)
    area([C.a1(1) C.a2(1)],[C.a1(2) C.a2(2)],-3,'FaceColor','k','FaceAlpha',0.3)
end
% plot(xs_1(1),xs_1(2),'b*','MarkerSize',5,'MarkerFaceColor','b')
% plot(xs_2(1),xs_2(2),'r*','MarkerSize',5,'MarkerFaceColor','r')
for i = 1:MC.N_x
    for j = 1:MC.N_y
        if (MC.result{i,j}==id_xs_1)
            plot(MC.ic{i,j}(1),MC.ic{i,j}(2),'bo',...
                'MarkerSize',5,'MarkerFaceColor','b','HandleVisibility','off')
        elseif (MC.result{i,j}==id_xs_2)
            plot(MC.ic{i,j}(1),MC.ic{i,j}(2),'ro',...
                'MarkerSize',5,'MarkerFaceColor','r','HandleVisibility','off')
        elseif (MC.result{i,j}==id_infeas)
            plot(MC.ic{i,j}(1),MC.ic{i,j}(2),'kx','MarkerSize',6,'HandleVisibility','off')
        else
            plot(MC.ic{i,j}(1),MC.ic{i,j}(2),'ko','MarkerSize',5,'HandleVisibility','off')
        end
    end
end
contour(X1,X2,C.f_cost(X,Y),-4:0.2:4,'k:')
set(gca,'Xlim',[-2,2],'Ylim',[-2,2]);
xlabel('$x_1$','FontSize',16)
ylabel('$x_2$','FontSize',16)
title('PTR Method','FontSize',16)

% penalty function
function val = penalty(x)
    val = norm(x,1);
end
