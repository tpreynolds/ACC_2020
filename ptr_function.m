function [S] = ptr_function(S,C)

if (S.plot)
    % initialize figure
    X1 = -3:0.01:2;
    X2 = -3:0.01:2;
    A   = 0:0.01:2*pi;
    [X,Y] = meshgrid(X1);
    figure(2), clf(2)
    figure(2), hold on
    plot(X1,C.f_eq(X1),'k','LineWidth',1)
    if (S.ineq)
        area([C.a1(1) C.a2(1)],[C.a1(2) C.a2(2)],-3,'FaceColor','k','FaceAlpha',0.3)
    end
    contour(X1,X2,C.f_cost(X,Y),-4:0.2:4,'k:')
    set(gca,'Xlim',[-2,2],'Ylim',[-2,2]);
    xlabel('$x_1$','FontSize',14)
    ylabel('$x_2$','FontSize',14)
    title('PTR Method','FontSize',16)
end

% loop through
if (~S.quiet)
    fprintf('PTR...\n')
    fprintf('Iter: 1 | x = (%0.2f,%0.2f) | s = %s | dJ = %s | eta: %s \n',...
            S.x(1),S.x(2),"    -    ","   -    ","   -    ");
end

for iter = 1:S.iter_max
    
    S.X_all(:,iter) = S.x;
    
    % plot the point
    if (S.plot)
        figure(2)
        plot(S.x(1),S.x(2),'ro','MarkerSize',5,'MarkerFaceColor','r')
        if (iter>1)
            plot(S.x(1)+S.eta(iter-1).*cos(A),S.x(2)+S.eta(iter-1).*sin(A),...
                    'b','LineWidth',1)
        end
        pause(0.1)
    end
    
    % solve the sub problem
    [f_k, dfdx_k] = f_vals(S.x,C.f_eq,C.c_eq);
    S = update_scales(S);
    
    cvx_begin quiet
        cvx_solver('ecos')
        variable x(2,1)
        variable s(2,1)
        variable eta nonnegative
        dual variable dual_mu;
    
        minimize( C.c'*(S.Px*x+S.qx) + S.w*penalty(s) + S.w_eta*penalty(eta) )
    
        subject to
    
        dual_mu : f_k + dfdx_k'*((S.Px*x+S.qx)-S.x) + s(1) == 0;
        if (S.ineq)
            dot(C.d_iq,(S.Px*x+S.qx)) >= C.b_iq;
        end
        S.x_min <= (S.Px*x+S.qx) <= S.x_max;
        (x-(S.Px\(S.x-S.qx)))'*(x-(S.Px\(S.x-S.qx))) <= eta;  
        
    cvx_end
    
    % compute (scaled) change in state
    S.dx(iter) = (x-(S.Px\(S.x-S.qx)))'*(x-(S.Px\(S.x-S.qx)));
    
    % rescale solution
    x = (S.Px*x+S.qx);
       
    % compute ratio
    [f_x, ~] = f_vals(x,C.f_eq,C.c_eq);
    f_iq_x = C.f_iq(x);
    S.dJ(iter) = S.J - (C.c' * x + S.w*(penalty(f_x) + penalty(max(f_iq_x,0))));
    S.dL(iter) = saturate(S.J - cvx_optval,1e-8,inf);
    rho = S.dJ(iter)/S.dL(iter);
    
    % update the non-linear cost
    S.J = C.c' * x + S.w*(penalty(f_x) + penalty(max(f_iq_x,0)));
    
    % update the solution
    S.x = x;
    S.eta(iter) = eta;
    
    % check feasibility
    defect = penalty(f_x) + max(f_iq_x,0);
    if (defect < S.feas_tol)
        S.feas = true;
    else
        S.feas = false;
    end
           
    % update the penalty weight(s) if called for
    if (~isempty(strfind(S.type,'tpr')))
        S.w_eta = 1.0/defect;
        if (S.w_eta<S.w_eta_min)
            S.w_eta = S.w_eta_min;
        end
    end
    if (strncmp(S.type,'dual',4))
        S.w = 15*abs(dual_mu);
    end

    % print progress
    if (~S.quiet)
        fprintf('Iter: %d | x = (%0.2f,%0.2f) | s = %2.2e | dJ = %2.2e | eta: %2.2e | feas: %u \n',...
            iter+1,S.x(1),S.x(2),penalty(s),S.dJ(iter),eta,S.feas);
    end
   
    % check convergence
    if ( (abs(S.dJ(iter)) < S.tol)&&(S.feas) )
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
        if (abs(S.dJ(iter))<S.tol)
            if ( (abs(S.dJ(iter-1))<S.tol) ...
                    && (abs(S.dJ(iter-2))<S.tol) ...
                    && (abs(S.dJ(iter-3))<S.tol) ...
                    && (abs(S.dJ(iter-4))<S.tol) ...
                    && (abs(S.dJ(iter-5))<S.tol) )
                if (~S.quiet)
                    fprintf('Insufficient Progress.\n')
                end
                break
            end
        end
    end      

end % main loop

end % ptr_function

% penalty function
function val = penalty(x)
    val = norm(x,1);
end

% main f function
function [f,dfdx] = f_vals(x,f_eq,c_eq)
    f = x(2) - f_eq(x(1));
    dfdx = [ -4*c_eq(1).*x(1).^3 + ...
             -3*c_eq(2).*x(1).^2 + ...
             -2*c_eq(3).*x(1) + ...
             -1*c_eq(4); 1 ]; 
end

% scaling function
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

% saturation to avoid degenerate cases
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
