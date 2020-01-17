function [P,S_sol,S_den] = SCvx(P)

% SCvx step 0.
fprintf('SCvx: Initializing\n')
P.alg = 'SCvx';
[P,S_sol,S_den,Ad,Bd,Cd,zd,Delta,i] = initialize(P);

% SCvx algorithm.
fprintf('SCvx: Starting\n')
while (i < P.SCvx.max_iter+1)
  % SCvx step 1.
  S_sol = socp(P,S_sol,Ad,Bd,Cd,zd,Delta,i);
  i = i+1;
  
  % SCvx step 2.
  dJ_i = S_sol(i-1).J - S_sol(i).J;
  dL_i = S_sol(i-1).J - S_sol(i).L;
  
  if (abs(dJ_i) < P.SCvx.dJ_tol) % Exit SCvx.
    fprintf(' dJ = %+07.2e |',dJ_i)
    fprintf(' dL = %+07.2e |',dL_i)
    fprintf(' r = %+07.2e |',1)
    fprintf(' SCvx exiting    |\n')
    S_sol(i).r = 1;
    S_sol(i).accepted = true;
    S_den(i) = propagate(P,S_sol(i),'full');
    break;
  else % Keep going.
    S_sol(i).r = dJ_i/abs(dL_i);
    fprintf(' dJ = %+07.2e |',dJ_i)
    fprintf(' dL = %+07.2e |',dL_i)
    fprintf(' r = %+07.2e |',S_sol(i).r)
    
    if (S_sol(i).r < P.SCvx.rho_0) % Reject step.
      fprintf(' SCvx rejecting  |\n')
      
      % Shrink trust region radius.
      Delta = Delta/P.SCvx.alpha;
      
      % Transfer sparse solution.
      S_sol(i).t = S_sol(i-1).t;
      S_sol(i).x = S_sol(i-1).x;
      S_sol(i).u = S_sol(i-1).u;
      S_sol(i).J = S_sol(i-1).J;
      S_sol(i).accepted = false;
      
      % Transfer dense solution.
      S_den(i) = S_den(i-1);
    else
      if (S_sol(i).r < P.SCvx.rho_1) % Shrink trust region radius.
        fprintf(' SCvx shrinking  |\n')
        Delta = Delta/P.SCvx.alpha;
      elseif (S_sol(i).r < P.SCvx.rho_2) % Keep trust region radius.
        fprintf(' SCvx keeping    |\n')
      else % Increase trust region radius.
        fprintf(' SCvx increasing |\n')
        Delta = P.SCvx.beta*Delta;
      end
      
      % Propagate new solution.
      [S_den(i),Ad,Bd,Cd,zd] = propagate(P,S_sol(i),'step');
      S_sol(i).accepted = true;
    end
  end
end

fprintf('SCvx: Done\n')