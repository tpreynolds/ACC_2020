function [S] = relerr_function(S,C)

% check initial feasibility
f_iq_x = C.f_iq(S.x);
f_eq_x = C.f_eq(S.x);
defect = penalty(f_eq_x) + max(f_iq_x,0);
if (defect < S.feas_tol)
    S.feas(1) = true;
else
    S.feas(1) = false;
end
S.cost(1) = C.c'*S.x;

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
    set(gca,'Xlim',[-2,2],'Ylim',[-2,2]);
    xlabel('$x_1$','FontSize',14)
    ylabel('$x_2$','FontSize',14)
    title('SCvx Method','FontSize',16)
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

% loop through
if (~S.quiet)
    fprintf('SCvx...\n')
    % print progress
    fprintf('Iter: 1 | x = (%0.2f,%0.2f) | s = %s | dJ = %s | dL = %s | rho = %s | step: %s, %s \n',...
        S.x(1),S.x(2),"    -    ","   -    ","   -    "," -  ","accepted","keep");
end
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
    
    % solve the sub problem
    [f_k, dfdx_k] = f_vals(S.x,C.f_eq,C.c_eq);
    S = update_scales(S);
    cvx_begin quiet
        cvx_solver('ecos')
        variable x(2,1)
        variable s(1,1)
        dual variable dual_mu;
    
        minimize( C.c'*(S.Px*x+S.qx) + S.w*penalty(s) )
    
        subject to
    
        dual_mu : f_k + dfdx_k'*((S.Px*x+S.qx)-S.x) + s(1) == 0;
        if (S.ineq)
            C.b_iq - dot(C.d_iq,(S.Px*x+S.qx)) <= 0;
        end
        
        S.x_min <= (S.Px*x+S.qx) <= S.x_max;
        
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
    x = (S.Px*x+S.qx);
       
    % compute ratio
    [f_x,~] = f_vals(x,C.f_eq,C.c_eq);
    if (S.ineq)
        temp = penalty(f_x)+max(C.b_iq-dot(C.d_iq,x),0);
        if (temp<1e-5)
            dJ_temp = 0.0;
        else
            dJ_temp = temp;
        end
        S.dJ(iter) = abs(1+C.c'*x + S.w*(dJ_temp) - cvx_optval-1);
    else
        S.dJ(iter) = abs(1+C.c'*x + S.w*(penalty(f_x)) - cvx_optval-1);
    end
    S.dL(iter) = saturate(abs(1+cvx_optval),1e-8,inf);
    rho = S.dJ(iter)/S.dL(iter);
    
    % update the linearized & non-linear costs
%     L_k = cvx_optval;
    
%     S.J = C.c'*x + S.w*penalty(f_x);
    
    % update trust region
    if (rho>S.rho2)
        S.r(iter+1) = S.r(iter) / S.alpha;
        S.reject(iter+1) = true;
        verdict = 'rejected';
        tr = 'shrink';
    else
        verdict = 'accepted';
        S.reject(iter+1) = false;
        if (rho>=S.rho1 && rho<=S.rho2)
            S.x = x;
            S.r(iter+1) = S.r(iter) / S.alpha;
            tr = 'shrink';
        elseif (rho>=S.rho0 && rho<S.rho1)
            S.x = x;
            S.r(iter+1) = S.r(iter);
            tr = 'keep';
        else
            S.x = x;
            S.r(iter+1) = S.beta * S.r(iter);
            tr = 'grow';
        end
    end
    %
    if (S.r(iter+1)<S.r_l) 
        S.r(iter+1) = S.r_l;
        tr = 'limitedL';
    end
    if (S.r(iter+1)>10)
        S.r(iter+1) = 10;
        tr = 'limitedU';
    end
    S.rho(iter) = saturate(rho,0,1);
           
    % update the penalty weight using dual variables if called for
    if (strfind(S.type,'dual'))
        S.w = 15 * abs(dual_mu);
    end
    
    % print progress
    if (~S.quiet)
        fprintf('Iter: %d | x = (%0.2f,%0.2f) | s = %2.2e | dJ = %2.2e | dL = %2.2e | rho = %0.2f | step: %s, %s \n',...
            iter+1,S.x(1),S.x(2),penalty(s),S.dJ(iter),S.dL(iter),rho,verdict,tr);
    end
   
    % check feasibility
    f_iq_x = C.f_iq(x);
    defect = penalty(f_x) + max(f_iq_x,0);
    if (defect < S.feas_tol)
        S.feas(iter+1) = true;
    else
        S.feas(iter+1) = false;
    end
    S.cost(iter+1) = C.c'*x;
        
    % check convergence
    if ( (norm(x-[0.52878234; -1.01920895]) < S.tol) && (S.feas(iter+1)) )
        if (~S.quiet)
            fprintf('Converged!\n')
        end
        S.cvrg_iter = iter;
        break
    else
        S.cvrg_iter = iter;
    end
    
    % break out after 5 iters of no progress
    if (iter>5)
        if (abs(S.dx(iter))<S.tol)
            if ( (abs(S.dx(iter-1))<S.tol) ...
                    && (abs(S.dx(iter-2))<S.tol) ...
                    && (abs(S.dx(iter-3))<S.tol) ...
                    && (abs(S.dx(iter-4))<S.tol) ...
                    && (abs(S.dx(iter-5))<S.tol) )
                if (~S.quiet)
                    fprintf('Insufficient Progress.\n')
                end
                break
            end
        end
    end

end % scvx loop

end

% penalty function
function val = penalty(x)
    val = norm(x,1);
end

function [f,dfdx] = f_vals(x,f_eq,c_eq)
    f = x(2) - f_eq(x(1));
    dfdx = [ -4*c_eq(1).*x(1).^3 + ...
             -3*c_eq(2).*x(1).^2 + ...
             -2*c_eq(3).*x(1) + ...
             -1*c_eq(4); 1 ]; 
end

function S = update_scales(S)
    Px = eye(numel(S.x));
    qx = zeros(size(S.x));
    
    for k = 1:numel(S.x)
       Px(k,k) = S.x_max(k) - S.x_min(k);
       qx(k) = 0.5*(S.x_max(k)+S.x_min(k));
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

