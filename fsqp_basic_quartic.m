clear variables;
set(0,'defaulttextinterpreter','latex',...
      'defaultAxesTickLabelInterpreter','latex',...
      'defaultLegendInterpreter','latex','defaultAxesFontSize',12)

%% setup

S = struct;
C = struct;

% FSQP parameters
S.alpha = 0.1;
S.beta  = 0.5;
S.delta = 2.0;
S.gamma = 1.0;
S.eta   = 0.1;
S.kappa = 2.1;
S.tau   = 2.5;
S.M     = 10.0;
S.w     = 1.5e0; % penalty weight
S.ineq  = true;
S.plot  = true;
S.quiet = false;
S.cvrg_iter = NaN;
S.iter_max  = 20;
S.tol       = 1e-3; % tolerance (convergence)
S.feas_tol  = 1e-3;

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
rts = roots([4*3*C.c_eq(1) 3*2*C.c_eq(2) 2*1*C.c_eq(3)]); % inflection pts

% initial guess and data matrices
S.x     = [0.78211;-0.9675];%[0.9471;-0.4676];
S.dx    = NaN(S.iter_max,1);
S.X_all = [S.x,NaN(2,S.iter_max-1)];

% Initialize figure
C.xs_1 = [-0.7372; 0.3163];   % local minimum
C.xs_2 = [0.5288; -1.0192];   % global minimum
if (S.plot)
    X1 = -3:0.01:2;
    X2 = -3:0.01:2;
    A = 0:0.01:2*pi;
    [X,Y] = meshgrid(X1);
    figure(1), clf(1)
    figure(1), hold on
    plot(X1,C.f_eq(X1),'k','LineWidth',1)
    if (S.ineq)
        area([C.a1(1) C.a2(1)],[C.a1(2) C.a2(2)],-3,'FaceColor','k','FaceAlpha',0.3)
    end
    contour(X1,X2,C.f_cost(X,Y),-4:0.2:4,'k:')
    plot(S.x(1),S.x(2),'ro','MarkerSize',5,'MarkerFaceColor','r')
    plot(C.xs_1(1),C.xs_1(2),'b*','MarkerSize',10);
    plot(C.xs_2(1),C.xs_2(2),'r*','MarkerSize',10);
    set(gca,'Xlim',[C.lb(1),C.ub(1)],'Ylim',[C.lb(2),C.ub(2)]);
    xlabel('$x_1$','FontSize',14)
    ylabel('$x_2$','FontSize',14)
    title('FSQP Method','FontSize',16)
end


%% Solve 

% print progress
fprintf('Iter: %02d  |  z = (%04.2f,%04.2f)  |  Dz = %02.2e\n',...
            1,S.x(1),S.x(2),0)

[~,~,Hx] = get_f(S.x,C);

% form the initial approximate Lagrangian
nu = 0.2;   % multiplier from the equality constraint
HL = -S.w * Hx + nu * Hx;

for iter = 2:S.iter_max
    
    [fx,dfdx,Hx]    = get_f(S.x,C);
    [gx,dgdx]       = get_g(S.x,C);
    [phi_c,dphi_c]  = get_phi(S.x,C,S.w);
    
    % solve for d0
    cvx_begin quiet
        cvx_solver('ecos')
        variable d0(2,1)
        dual variable dual_lam

        minimize( 0.5*d0'*HL*d0 + dphi_c'*d0 )

        fx + dfdx'*d0 <= 0.0;
        dual_lam : gx + dgdx'*d0 <= 0.0;
    cvx_end

    % solve for d1 and gamma
    cvx_begin quiet
        cvx_solver('ecos')
        variable d1(2,1)
        variable gamm

        minimize( 0.5*S.eta*dot(d0-d1,d0-d1)+gamm )

        dphi_c' * d1 <= gamm;
        fx + dfdx' * d1 <= gamm;
        gx + dgdx' * d1 <= gamm;
    cvx_end

    % compute rho_k and d_k
    temp    = norm(d0)^(S.kappa);
    temp2   = max(0.5,norm(d1)^(S.tau));
    rho     = temp/(temp+temp2);
    d       = (1-rho)*d0 + rho*d1;
    
    % compute d_tilde
    [fd,~,~]    = get_f(S.x+d,C);
    [gd,~]      = get_g(S.x+d,C);
    temp        = -min(0.01*norm(d),norm(d)^(S.tau));
    
    cvx_begin quiet
        cvx_solver('ecos')
        variable dd(2,1)

        minimize( 0.5*(d+dd)'*HL*(d+dd) + dphi_c'*(d+dd) )

        fd + dfdx' * dd <= -temp;
        gd + dgdx' * dd <= -temp;
    cvx_end

    if (norm(dd) > norm(d))
        dd = zeros(2,1);
    end
    
    % arc search
    t     = 1.0;
    temp1 = 1.0;
    temp2 = 1.0;
    % loop
    for k_ = 1:10
        z = S.x + t*d + t^2*dd;
        [phi_c_z,~] = get_phi(z,C,S.w);
        [f_x_,~,~]  = get_f(z,C);
        temp1 = (phi_c_z - phi_c - S.alpha * t * dot(dphi_c,d) <= 0.0 );
        temp2 = (f_x_ <= 0.0);
        
        if (temp1(1) && temp2(1))
            break;
        else
            t = S.beta * t;
        end
    end
    
    % solve for mu
    cvx_begin quiet
    cvx_solver('ecos')
    variable mub
    minimize( (C.c-dual_lam(1)*C.d_iq-mub*dfdx)'*(C.c-dual_lam(1)*C.d_iq-mub*dfdx) )
    cvx_end
    
    % update penalty parameter
    if ( (S.w+mub >= gamm) || (S.w*norm(HL*d0)>S.M) )
        S.w = S.w;
    else
        S.w = max(gamm-mub,S.delta*S.w);
    end
    
    % update the Hessian approximation
    [~,~,Hz] = get_f(z,C);
    
    % update the Hessian
    if (z(1)<=rts(1) || z(1)>=rts(2))
        HL = -S.w * Hz + nu * Hz;
    else
        HL = zeros(2,2);
    end

    % print progress
    fprintf('Iter: %02d  |  z = (%04.2f,%04.2f)  |  Dz = %02.2e\n',...
        iter,z(1),z(2),norm(z-S.x))
    
    % save the iterate
    S.dx(iter)      = norm(S.x - z);
    S.x             = z;
    S.X_all(:,iter) = z;
    
    if (S.plot)
        % plot the point
        figure(1)
        plot(S.x(1),S.x(2),'ro','MarkerSize',5,'MarkerFaceColor','r')
        pause(0.1)
    end
    
    % exit criteria
    if (S.dx(iter) < S.tol)
        S.cvrg_iter = iter;
        fprintf('Converged.\n')
        break;
    end
    
end
    

%% PLOT

% colors
ic_col = [0,0.6,1];%./255;
fc_col = [1,0,1];%./255;
for k = 1:3
    N_col(k,:) = linspace(ic_col(k),fc_col(k),S.cvrg_iter);
end
X1 = -3:0.01:2;
X2 = -3:0.01:2;
[X,Y] = meshgrid(X1);

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
h4=plot(C.xs_1(1),C.xs_1(2),'b*','MarkerSize',10);
h5=plot(C.xs_2(1),C.xs_2(2),'r*','MarkerSize',10);
set(gca,'FontSize',14)
legend('equality constraint','inequality constraint','cost lines',...
    'local min','global min','Location','northwest','FontSize',12)
colormap(N_col')
hC = colorbar('TickLabelInterpreter','latex','FontSize',14);
hC.Label.String = 'Iteration';
hC.Label.Interpreter = 'latex';
hC.Label.FontSize = 14;
caxis([0 min(S.cvrg_iter,S.iter_max)]);
set(gca,'Xlim',[-2,2],'Ylim',[-2,2]);
xlabel('$z_1$','FontSize',16)
ylabel('$z_2$','FontSize',16)


%% HELPER FUNCTIONS

% penalty function
function val = penalty(x)
    val = norm(x,1);
end

function [fz,dfdz,Hfz] = get_f(z,C)
    fz = z(2) - C.f_eq(z(1));
    dfdz = [ -4*C.c_eq(1)*z(1)^3-3*C.c_eq(2)*z(1)^2-2*C.c_eq(3)*z(1)-1*C.c_eq(4);
                1 ];
    Hfz = [ -4*3*C.c_eq(1)*z(1)^2-3*2*C.c_eq(2)*z(1)-2*1*C.c_eq(3)    0;
                                0                               0 ];
end

function [gz,dgdz] = get_g(z,C)

gz      = zeros(5,1);
dgdz    = zeros(2,5);

gz(1)   = C.b_iq - dot(C.d_iq,z);
gz(2:3) = z - C.ub;
gz(4:5) = C.lb - z;

dgdz(:,1)       = -C.d_iq;
dgdz(:,2:3)     =  eye(2);
dgdz(:,4:5)     = -eye(2);

end

function [phi_c,dphi_c] = get_phi(z,C,w)

[fz,dfdz,~] = get_f(z,C);

phi_c   = dot(C.c,z) - w * fz;
dphi_c  = C.c - w * dfdz;

end