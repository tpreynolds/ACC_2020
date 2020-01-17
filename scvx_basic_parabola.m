clear variables; close all;
set(0,'defaulttextinterpreter','latex')

% SCvx parameters
rho0 = 0.0;
rho1 = 0.1; 
rho2 = 0.9; 
alpha = 1.5; % shrink rate
beta = 2.0; % growth rate
r_l = 0.001; % minimum trust region size
w = 1e3; % penalty weight
type = 'inf';
ineq = true;
c = [1;1]; % cost function
a1 = [1; -0.1]; % a1 and a2 control the linear inequality
a2 = [-1; 0.6];
m = (a2(2)-a1(2))/(a2(1)-a1(1)); % slope
b = a2(2) - m*a2(1); % y-intercept
d = [-m;1]; % constraint is d' * x >= b;

iter_max = 50;
tol = 1e-5; % tolerance for dJ_k (convergence)
small = 1e-8;
x_true = [-0.5,0.25]; % true solution for no inequality constraint

% initial guess
x_k = [0.2;0.8];
r_k = 0.1; eta = r_k;
J_k = c'*x_k + w*penalty(x_k);

% initialize figure
X1 = -1:0.01:1;
X2 = -1:0.01:1;
A   = 0:0.01:2*pi;
constraint = @(x)(x.^2); % constraint function
cost = @(x,y)(c(1)*x+c(2)*y); % cost function
[X,Y] = meshgrid(X1);
figure(1), hold on
plot(X1,constraint(X1),'k','LineWidth',1)
if (ineq)
    area([a1(1) a2(1)],[a1(2) a2(2)],-1,'FaceColor','k','FaceAlpha',0.3)
end
contour(X1,X2,cost(X,Y),-2:0.2:2,'k:')
set(gca,'Xlim',[-1,1],'Ylim',[-0.6,1]);
xlabel('$x_1$','FontSize',14)
ylabel('$x_2$','FontSize',14)
title('SCvx Method','FontSize',16)

% loop through
fprintf('SCvx...\n')
% print progress
fprintf('Iter: 1 | x = (%0.2f,%0.2f) | s = %s | dJ = %s | dL = %s | rho = %s | step: %s, %s \n',...
        x_k(1),x_k(2),"    -    ","   -    ","   -    "," -  ","accepted","keep");
for iter = 1:iter_max
    
    % plot the point
    figure(1)
    plot(x_k(1),x_k(2),'ro','MarkerSize',5,'MarkerFaceColor','r')
    switch type
        case 'inf'
            rectangle('Position',[x_k(1)-r_k x_k(2)-r_k 2*r_k 2*r_k],...
                'LineWidth',1,'EdgeColor','b')
        case '1'
            plot([x_k(1),x_k(1)+r_k,x_k(1),x_k(1)-r_k,x_k(1)],...
                 [x_k(2)+r_k,x_k(2),x_k(2)-r_k,x_k(2),x_k(2)+r_k],...
                 'b','LineWidth',1)
        otherwise
            plot(x_k(1)+r_k.*cos(A),x_k(2)+r_k.*sin(A),'b','LineWidth',1)
    end
    pause(0.1)
    
    % solve the sub problem
    [f_k, dfdx_k] = f_vals(x_k);
    cvx_begin quiet
        cvx_solver('ecos')
        variable x(2,1)
        variable s
        dual variable mu;
    
        minimize( c'*x + w*penalty(s) )
    
        subject to
    
        mu : f_k + dfdx_k'*(x-x_k) + s == 0;
        if (ineq)
            dot(d,x) >= b;
        end
        switch type
            case '2'
                norm(x-x_k,2) <= r_k;
            case 'inf'
                norm(x-x_k,inf) <= r_k;
            case '1'
                norm(x-x_k,1) <= r_k;
            otherwise
                (x-x_k)'*(x-x_k) <= r_k^2;
        end      
    cvx_end
       
    % compute ratio
    dJk = J_k - (c'*x + w*penalty(x(2)-x(1)^2));
    dLk = saturate(J_k - cvx_optval,small,inf);
    rho_k = dJk/dLk;
    
    % update the linearized & non-linear costs
    L_k = cvx_optval;
    J_k = c'*x + w*penalty(x(2)-x(1)^2);
    
    % update trust region
    if (rho_k<rho0)
        r_k = r_k / alpha;
        verdict = 'rejected';
        tr = 'shrink';
    else
        verdict = 'accepted';
        if (rho_k<rho1)
            x_k = x;
            r_k = r_k / alpha;
            tr = 'shrink';
        elseif (rho_k>=rho1 && rho_k<rho2)
            x_k = x;
            r_k = r_k;
            tr = 'keep';
        else
            x_k = x;
            r_k = beta * r_k;
            tr = 'grow';
        end
    end
    if (r_k<r_l) 
        r_k = r_l;
        tr = 'limited';
    end
    if (r_k>10)
        r_k = 10;
        tr = 'limited';
    end
        
    % update the penalty weight using dual variables
%     if (mu>0)
%         w = 15*abs(mu);
%     end
    
    % print progress
    fprintf('Iter: %d | x = (%0.2f,%0.2f) | s = %2.2e | dJ = %2.2e | dL = %2.2e | rho = %0.2f | step: %s, %s \n',...
            iter+1,x_k(1),x_k(2),s(1),dJk,dLk,rho_k,verdict,tr);
   
    % check convergence
    if ((abs(dJk) < tol)&&(w*penalty(x(2)-x(1)^2))<1e-3)
        fprintf('Converged!\n')
        break
    end

end % scvx loop

% penalty function
function val = penalty(x)
    val = norm(x,1);
end

function [f,dfdx] = f_vals(x)
    f = x(2) - x(1).^2;
    dfdx = [ -2.*x(1); 1 ]; 
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
