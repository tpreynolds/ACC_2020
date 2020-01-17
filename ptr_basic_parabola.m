clear variables; close all;
set(0,'defaulttextinterpreter','latex')

% SCvx parameters
rho0 = 0.0;
rho1 = 0.1; 
rho2 = 0.9; 
alpha = 2.0; % shrink rate
beta = 2.0; % growth rate
r_l = 0.001; % minimum trust region size
w = 1e3; % vc weight
w_eta = 1e0; % tr weight
c = [1;1]; % cost function
iter_max = 50;
tol = 1e-6; % tolerance for dJ_k (convergence)
x_true = [-0.5,0.25]; % true solution

% initial guess
x_k = [0.4;0.4];
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
contour(X1,X2,cost(X,Y),-2:0.2:2,'k:')
set(gca,'Xlim',[-1,1],'Ylim',[-0.6,1]);
xlabel('$x_1$','FontSize',14)
ylabel('$x_2$','FontSize',14)
title('FTR Method','FontSize',16)

% loop through
fprintf('SCvx...\n')
for iter = 1:iter_max
    
    % plot the point
    figure(1)
    plot(x_k(1),x_k(2),'ro','MarkerSize',5,'MarkerFaceColor','r')
    plot(x_k(1)+eta.*cos(A),x_k(2)+eta.*sin(A),'b','LineWidth',1)
    pause(0.1)
    
    % solve the sub problem
    cvx_begin quiet
%         cvx_solver('ecos')
        variable x(2,1)
        variable s
        variable eta nonnegative
        dual variable mu;
    
        minimize( c'*x + w*penalty(s) + w_eta*penalty(eta) )
    
        subject to
    
        mu : x_k(1)^2 - 2*x_k(1)*x(1) + x(2) + s == 0;
        (x-x_k)'*(x-x_k) <= eta;
    
    cvx_end
       
    % compute ratio
    dJk = J_k - (c'*x + w*penalty(x(2)-x(1)^2));
    dLk = J_k - cvx_optval;
    rho_k = dJk/dLk;
    
    % update the linearized & non-linear costs
    L_k = cvx_optval;
    J_k = c'*x + w*penalty(x(2)-x(1)^2);
    
    % update the solution
    x_k = x;
        
    % update the penalty weight
%     if (mu>0)
%         w = 10*mu;
%     end
    
    % print progress
    fprintf('Iter: %d | x = (%0.2f,%0.2f) | s = %2.2e | dJk = %2.2e | dLk = %2.2e | rhok = %0.2f | eta = %2.2e \n',...
            iter,x_k(1),x_k(2),s,dJk,dLk,rho_k,eta);
   
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
