function S = hybrid_function(S,C)

% check initial feasibility
f_iq_x = C.f_iq(S.x);
f_eq_x = C.f_eq(S.x);
defect = penalty(f_eq_x) + max(f_iq_x,0);
if (defect < S.feas_tol)
    S.feas(1) = true;
else
    S.feas(1) = false;
end

% get trust region type (may contain the string 'dual' somewhere)
if (strfind(S.type,'inf'))
    type = 'inf';
elseif (strfind(S.type,'1'))
    type = '1';
elseif (strfind(S.type,'2'))
    type = '2';
else
    type = '2sqr';
end

% Initialize figure
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
    title('Hybrid Method','FontSize',15)
end

fprintf('SCP...using %s\n',S.scp_algo)
fprintf('Iter: %03d | x = (%0.2f,%0.2f) | s = %s | dJ = %s | dL = %s | rho = %s | step: %s, %s \n',...
    1,S.x(1),S.x(2),"    -    ","   -    ","   -    "," -  ","accepted","keep");

% main loop
algo = 'scp';
for iter = 1:S.iter_max
   
    S.X_all(:,iter) = S.x;
    if (S.plot)
        % plot the point
        figure(1)
        plot(S.x(1),S.x(2),'ro','MarkerSize',5,'MarkerFaceColor','r')
        switch type
            case 'inf'
                rectangle('Position',[S.x(1)-S.r(iter) S.x(2)-S.r(iter) 2*S.r(iter) 2*S.r(iter)],...
                    'LineWidth',1,'EdgeColor','b')
            case '1'
                plot([S.x(1),S.x(1)+S.r(iter),S.x(1),S.x(1)-S.r(iter),S.x(1)],...
                    [S.x(2)+S.r(iter),S.x(2),S.x(2)-S.r(iter),S.x(2),S.x(2)+S.r(iter)],...
                    'b','LineWidth',1)
            otherwise
                plot(S.x(1)+S.r(iter).*cos(A),S.x(2)+S.r(iter).*sin(A),'b','LineWidth',1)
        end
        pause(0.1)
    end
    
    % step
    switch algo
        case 'scp'
            [S,crawl] = scp_step(S,C,type,iter);
        case 'fsqp'
            [S,crawl] = fsqp_step(S,C,iter);
%             fprintf('Back to SCP...\n');
    end
    
    % crawling?
    if (crawl)
        fprintf('Crawling detected, taking FSQP step...\n');
        algo    = 'fsqp';
%     else
%         algo = 'scp';
    end
    
    % converged
    if ~(isnan(S.cvrg_iter))
       break; 
    end
    
end

end

%% Extra functions

function S = update_scales(S,C)
    Px = eye(numel(S.x));
    qx = zeros(size(S.x));
    
    for k = 1:numel(S.x)
       Px(k,k)  = C.ub(k) - C.lb(k);
       qx(k)    = 0.5*( C.ub(k)+C.lb(k) );
    end
    
    S.Px = Px;
    S.qx = qx;
end

function y = saturate(u,low,high)
    n = numel(u);
    y = u;
    for i = 1:n
        if (abs(y(i))>high)
            y(i) = high;
        end
        if (abs(y(i))<low)
            y(i) = low;
        end
    end
end

function [fz,dfdz,Hfz] = get_f(z,C)
    fz   = z(2) - C.f_eq(z(1));
    dfdz = [ -4*C.c_eq(1)*z(1)^3-3*C.c_eq(2)*z(1)^2-2*C.c_eq(3)*z(1)-1*C.c_eq(4);
                1 ];
    Hfz  = [ -4*3*C.c_eq(1)*z(1)^2-3*2*C.c_eq(2)*z(1)-2*1*C.c_eq(3)    0;
                                0                                      0 ];
end

function [gz,dgdz] = get_g(z,C)
    switch C.ineq
        case true
            gz          = zeros(5,1);
            dgdz        = zeros(2,5);
            gz(5)       = C.b_iq - dot(C.d_iq,z);
            dgdz(:,5)   = -C.d_iq;
        case false
            gz      = zeros(4,1);
            dgdz    = zeros(2,4);
    end

    gz(1:2) = z - C.ub;
    gz(3:4) = C.lb - z;

    dgdz(:,1:2)     =  eye(2);
    dgdz(:,3:4)     = -eye(2);
end

function [phi_c,dphi_c] = get_phi(z,C,w)
    [fz,dfdz,~] = get_f(z,C);
    phi_c       = dot(C.c,z) - w * fz;
    dphi_c      = C.c - w * dfdz;
end

function [S,crawl] = scp_step(S,C,type,iter)

    crawl = false;

    % solve the sub problem
    [fx,dfdx,~] = get_f(S.x,C);
    [gx,dgdx]   = get_g(S.x,C);
    S = update_scales(S,C);
    cvx_begin quiet
        cvx_solver('ecos')
        variable x(2,1)
        variable s(1)
        dual variable dual_mu;
    
        minimize( C.c'*(S.Px*x+S.qx) + S.w*penalty(s) )
    
        subject to
    
        dual_mu : fx + dfdx'*((S.Px*x+S.qx)-S.x) + s(1) == 0.0;
        gx + dgdx'*((S.Px*x+S.qx)-S.x) <= 0.0;
                
        switch type
            case '2'
                norm((S.Px*x+S.qx)-S.x,2) <= S.r(iter);
            case 'inf'
                norm((S.Px*x+S.qx)-S.x,inf) <= S.r(iter);
            case '1'
                norm((S.Px*x+S.qx)-S.x,1) <= S.r(iter);
            otherwise
                ((S.Px*x+S.qx)-S.x)'*((S.Px*x+S.qx)-S.x) <= S.r(iter)^2;
        end      
    cvx_end
    
    switch type
        case '2'
            S.dx(iter) = norm((S.Px*x+S.qx)-S.x,2);
        case 'inf'
            S.dx(iter) = norm((S.Px*x+S.qx)-S.x,inf);
        case '1'
            S.dx(iter) = norm((S.Px*x+S.qx)-S.x,1);
        otherwise
            S.dx(iter) = sqrt(((S.Px*x+S.qx)-S.x)'*((S.Px*x+S.qx)-S.x));
    end
    
    % rescale solution
    x       = (S.Px*x+S.qx);
    S.mu    = dual_mu;
       
    switch S.scp_algo
        
        case 'rel_dec'
            % compute ratio
            [fz,~,~] = get_f(x,C);
            if (S.ineq)
                S.dJ(iter) = S.J - (C.c'*x + S.w*(penalty(fz)+max(C.b_iq-dot(C.d_iq,x),0)));
            else
                S.dJ(iter) = S.J - (C.c'*x + S.w*(penalty(fz)));
            end
            S.dL(iter) = saturate(S.J - cvx_optval,1e-8,inf);
            rho = S.dJ(iter)/S.dL(iter);
            S.J = C.c'*x + S.w*penalty(fz);
            
            % update trust region
            if (rho<S.rho0)
                S.r(iter+1)     = S.r(iter) / S.tr_shrink;
                S.rej(iter+1)   = true;
                verdict         = 'rejected';
                tr              = 'shrink';
            else
                verdict         = 'accepted';
                S.rej(iter+1)   = false;
                if (rho<S.rho1)
                    S.x         = x;
                    S.r(iter+1) = S.r(iter) / S.tr_shrink;
                    tr          = 'shrink';
                elseif (rho>=S.rho1 && rho<S.rho2)
                    S.x         = x;
                    S.r(iter+1) = S.r(iter);
                    tr          = 'keep';
                else
                    S.x         = x;
                    S.r(iter+1) = S.tr_grow * S.r(iter);
                    tr          = 'grow';
                end
            end
            
        case 'rel_err'
            % compute ratio
            [fz,~,~] = get_f(x,C);
            if (S.ineq)
                dJ_temp     = penalty(fz)+max(C.b_iq-dot(C.d_iq,x),0);
                S.dJ(iter)  = abs(1+C.c'*x+S.w*(dJ_temp)-cvx_optval-1);
            else
                S.dJ(iter)  = abs(1+C.c'*x+S.w*(penalty(fz))-cvx_optval-1);
            end
            S.dL(iter) = saturate(abs(1+cvx_optval),1e-8,inf);
            rho = S.dJ(iter)/S.dL(iter);
                       
            % update trust region
            if (rho>S.rho2)
                S.r(iter+1)         = S.r(iter) / S.tr_shrink;
                S.rej(iter+1)       = true;
                verdict             = 'rejected';
                tr                  = 'shrink';
            else
                verdict             = 'accepted';
                S.rej(iter+1)       = false;
                if (rho>=S.rho1 && rho<=S.rho2)
                    S.x         = x;
                    S.r(iter+1) = S.r(iter) / S.tr_shrink;
                    tr          = 'shrink';
                elseif (rho>=S.rho0 && rho<S.rho1)
                    S.x         = x;
                    S.r(iter+1) = S.r(iter);
                    tr          = 'keep';
                else
                    S.x         = x;
                    S.r(iter+1) = S.tr_grow * S.r(iter);
                    tr          = 'grow';
                end
            end           
    end % switch case scp_algo
    % saturate trust region
    if (S.r(iter+1)<S.r_l)
        S.r(iter+1) = S.r_l;
        tr = 'limitedL';
    end
    if (S.r(iter+1)>10)
        S.r(iter+1) = 10;
        tr = 'limitedU';
    end
    S.rho(iter+1) = saturate(rho,0,1);
    
    % print progress
    fprintf('Iter: %03d | x = (%0.2f,%0.2f) | s = %2.2e | dJ = %2.2e | dL = %2.2e | rho = %0.2f | step: %s, %s \n',...
            iter+1,S.x(1),S.x(2),penalty(s),S.dJ(iter),S.dL(iter),rho,verdict,tr);
   
    % check feasibility
    [gz,~] = get_g(x,C);
    defect = penalty(fz) + sum(max(gz,0));
    if (defect < S.feas_tol)
        S.feas(iter+1) = true;
    else
        S.feas(iter+1) = false;
    end
    S.cost(iter+1) = C.c'*x;
        
    % check convergence
    if ( (abs(S.dJ(iter)) < S.tol) && (S.feas(iter+1)) )
        fprintf('Converged!\n')
        S.cvrg_iter = iter;
    end
    
    % check crawling
    if (iter > 3)
        temp1 = strcmp(tr,'shrink');
        temp2 = (~(S.rej(iter)||S.rej(iter-1)||S.rej(iter-2)));
        temp31 = (S.rho(iter)<S.rho_ub&&S.rho(iter)>S.rho_lb);
        temp32 = (S.rho(iter-1)<S.rho_ub&&S.rho(iter-1)>S.rho_lb);
        temp33 = (S.rho(iter-2)<S.rho_ub&&S.rho(iter-2)>S.rho_lb);
        if ( temp1 && temp2 && temp31 && temp32 && temp33 )
            if (S.prev_event && S.r(iter+1)<S.r(1))
                crawl = true;               % this has happened before
            else
                S.prev_event = true;
            end
        end
    end
    
end % scvx_step

function [S,crawl] = fsqp_step(S,C,iter)

    [fx,dfdx,Hx]    = get_f(S.x,C);
    [gx,dgdx]       = get_g(S.x,C);
    [phi_c,dphi_c]  = get_phi(S.x,C,S.wsqp);
    
    % form the initial approximate Lagrangian
    if (S.x(1)<=C.rts(1) || S.x(1)>=C.rts(2))
        HL = -S.wsqp * Hx + S.mu * Hx;
    else
        HL = zeros(2,2);
    end
    
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
    [fd,dfdd,~]    = get_f(S.x+d,C);
    [gd,dgdd]      = get_g(S.x+d,C);
    temp        = -min(0.01*norm(d),norm(d)^(S.tau));
    
    cvx_begin quiet
        cvx_solver('ecos')
        variable dd(2,1)

        minimize( 0.5*(d+dd)'*HL*(d+dd) + dphi_c'*(d+dd) )

        fd + dfdd' * dd <= -temp;
        gd + dgdd' * dd <= -temp;
    cvx_end

    if (norm(dd) > norm(d))
        dd = zeros(2,1);
    end
    
    % arc search
    t = 1.0;
    for k_ = 1:10
        z = S.x + t*d + t^2*dd;
        [phi_c_z,~] = get_phi(z,C,S.wsqp);
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
        minimize( (C.c-dual_lam(5)*C.d_iq-mub*dfdx)'*(C.c-dual_lam(5)*C.d_iq-mub*dfdx) )
    cvx_end
    
    % update penalty parameter
    if ( (S.wsqp+mub >= gamm) || (S.wsqp*norm(HL*d0)>S.M) )
        S.wsqp = S.wsqp;
    else
        S.wsqp = max(gamm-mub,S.delta*S.wsqp);
    end
    
    % check feasibility
    [g_x_,~]    = get_g(z,C);
    defect      = penalty(f_x_) + sum(max(g_x_,0));
    if (defect < S.feas_tol)
        S.feas(iter+1) = true;
    else
        S.feas(iter+1) = false;
    end
        
    % print progress
    fprintf('Iter: %03d | z = (%04.2f,%04.2f) | Dz = %02.2e | feas = %d\n',...
                iter+1,z(1),z(2),norm(z-S.x),S.feas(iter+1))
   
    % save the iterate
    S.dx(iter+1)        = norm(S.x - z);
    S.x                 = z;
    S.X_all(:,iter+1)   = z;
    S.cost(iter+1)      = C.c'*z; 
    S.r(iter+1)         = S.r(iter);
    S.J                 = C.c'*z + S.w*defect;
    S.rej(iter+1)       = false;
    
    if (S.plot)
        % plot the point
        figure(1)
        plot(S.x(1),S.x(2),'ro','MarkerSize',5,'MarkerFaceColor','r')
        pause(0.1)
    end
    
    % exit criteria
    if (S.dx(iter+1) < S.tol)
        S.cvrg_iter = iter;
        fprintf('Converged.\n')
    end
    % output crawl
    crawl = false;

end % fsqp_function
