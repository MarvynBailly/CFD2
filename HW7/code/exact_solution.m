function [rho, u, p, e] = exact_solution(x, gamma, M_start, p_start, rho_start, shock, x_shock)
    len = length(x);
    A = 1.398 + 0.347 * tanh(0.8 * (x - 4.0));
    
    % 1. Initial total pressure and density
    p0 = p_start * (1 + (gamma-1)/2 * M_start^2)^(gamma/(gamma-1));
    rho0 = rho_start * (1 + (gamma-1)/2 * M_start^2)^(1/(gamma-1));
    
    % 2. Compute A*
    Astar = A(1) * M_start * (2/(gamma+1) * (1 + (gamma-1)/2 * M_start^2))^(-(gamma+1)/(2*(gamma-1)));
    
    % 3. Set initial values
    M = zeros(1, len);
    p = zeros(1, len);
    rho = zeros(1, len);
    u = zeros(1, len);
    
    M(1) = M_start;
    p(1) = p_start;
    rho(1) = rho_start;
    u(1) = M_start * sqrt(gamma * p_start / rho_start);
    
    [~, shock_index] = min(abs(x - x_shock));


    % 4. March downstream
    for j = 2:len
        % Solve Area-Mach relation
        A_ratio = A(j) / Astar;
        func = @(Mval) (1/Mval) * (2/(gamma+1)*(1+(gamma-1)/2*Mval^2))^((gamma+1)/(2*(gamma-1))) - A_ratio;
    
        % Use previous Mach as guess
        M_guess = M(j-1);
        M(j) = fzero(func, M_guess);
    
        % Update primitive variables (isentropic relations)
        p(j) = p0 / (1 + (gamma-1)/2 * M(j)^2)^(gamma/(gamma-1));
        rho(j) = rho0 / (1 + (gamma-1)/2 * M(j)^2)^(1/(gamma-1));
        u(j) = M(j) * sqrt(gamma * p(j) / rho(j));
    
        % If shock at this location
        if (shock == 1) && (j == shock_index)
            disp('Shock detected!');
            % App_starty normal shock jump relations
            M1 = M(j);
            M2 = sqrt( (2 + (gamma-1)*M1^2) / (2 * gamma * M1^2 - (gamma-1)) );
            M(j) = M2;
    
            % Update post-shock stagnation quantities
            Astar = A(j) * M2 * (2/(gamma+1)*(1+(gamma-1)/2*M2^2))^(-(gamma+1)/(2*(gamma-1)));
    
            p(j) = p(j) * (1 + 2*gamma/(gamma+1)*(M1^2 - 1));
            rho(j)= rho(j) * ((gamma+1)*M1^2) / (2 + (gamma-1)*M1^2);
            u(j) = u(j) / ((gamma+1)*M1^2 / (2 + (gamma-1)*M1^2));
            p0 = p(j) * (1 + (gamma-1)/2 * M2^2)^(gamma/(gamma-1));
            rho0 = rho(j) * (1 + (gamma-1)/2 * M2^2)^(1/(gamma-1));
        end
    end
    
    % 5. Compute energy
    e = p/(gamma-1) + 0.5 * rho .* u.^2;
    
    end

% function [rho, u, p, e, M, A] = exact_solution(x, gamma, M0, p0, rho0, x_shock)
%     % Compute exact solution with or without a shock
    
%     Nx = length(x);
    
%     % Area function
%     A = calcarea(x);
    
%     % Upstream A*
%     A_star = A(1) * M0 * (2/(gamma+1)*(1+(gamma-1)/2*M0^2))^(-(gamma+1)/(2*(gamma-1)));
    
%     % Allocate
%     M = zeros(1, Nx);
%     rho = zeros(1, Nx);
%     u = zeros(1, Nx);
%     p = zeros(1, Nx);
    
%     % Determine if a shock should be app_startied
%     app_starty_shock = (x_shock > min(x)) && (x_shock < max(x));
    
%     if app_starty_shock
%         % Find shock index
%         [~, j_shock] = min(abs(x - x_shock));
%     else
%         j_shock = Nx + 1; % No shock
%     end
    
%     % Set initial Mach number
%     M(1) = M0;
    
%     % March along x
%     for j = 2:Nx
%         if j < j_shock
%             % Before shock: supersonic branch
%             A_ratio = A(j) / A_star;
%             M(j) = solve_mach_number(A_ratio, gamma, true);
%         elseif j == j_shock && app_starty_shock
%             % At shock: app_starty shock jump
%             M1 = M(j-1);
%             M2_sq = (1 + 0.5*(gamma-1)*M1^2) / (gamma*M1^2 - 0.5*(gamma-1));
%             M(j) = sqrt(M2_sq);
    
%             % Update downstream A*
%             A_star = A(j) * M(j) * ( (2/(gamma+1)*(1+(gamma-1)/2*M(j)^2))^(-(gamma+1)/(2*(gamma-1))) );
%         else
%             % After shock or full domain (no shock): subsonic branch
%             A_ratio = A(j) / A_star;
%             M(j) = solve_mach_number(A_ratio, gamma, false);
%         end
%     end
    
%     % Now compute primitive variables
%     for j = 1:Nx
%         if j < j_shock
%             p_p0 = (1 + 0.5*(gamma-1)*M(j)^2)^(-gamma/(gamma-1));
%             rho_rho0 = (1 + 0.5*(gamma-1)*M(j)^2)^(-1/(gamma-1));
    
%             p(j) = p0 * p_p0;
%             rho(j) = rho0 * rho_rho0;
%         elseif j == j_shock && app_starty_shock
%             % App_starty shock jump in pressure and density
%             M1 = M(j-1);
%             % p_ratio = (1 + 2*gamma/(gamma+1)*(M1^2 - 1));
%             % rho_ratio = ((gamma+1)*M1^2) / (2 + (gamma-1)*M1^2);

%             p_ratio = (1 + 2*gamma/(gamma+1)*(M1^2 - 1));
%             rho_ratio = ((gamma+1)*M1^2) / (2 + (gamma-1)*M1^2);
    
%             p(j) = p(j-1) * p_ratio;
%             rho(j) = rho(j-1) * rho_ratio;
%         else
%             % After shock or no shock
%             p_p2 = (1 + 0.5*(gamma-1)*M(j)^2)^(-gamma/(gamma-1));
%             rho_rho2 = (1 + 0.5*(gamma-1)*M(j)^2)^(-1/(gamma-1));
    
%             p(j) = p(j-1) * p_p2;
%             rho(j) = rho(j-1) * rho_rho2;
%         end
    
%         a = sqrt(gamma * p(j) / rho(j));
%         u(j) = M(j) * a;
%     end
    
%     % Compute energy
%     e = p/(gamma-1) + 0.5 * rho .* u.^2;
    
%     end




% function M = solve_mach_number(A_ratio, gamma, supersonic)
%     % Solve for M given area ratio A/A*
    
%     fun = @(Mval) -(1/Mval) * (2/(gamma+1)*(1 + (gamma-1)/2*Mval^2))^((gamma+1)/(2*(gamma-1))) + A_ratio;
    
%     if supersonic
%         interval = [1.0001, 5];
%     else
%         interval = [0.01, 0.999];
%     end


%     M = fzero(fun, interval);
    
% end
        

% function [rho, u, p, x] = exact_solution(rho_sp, u_sp, e_sp, method, gamma, fsmach, p0, rho0, xsh)
%     jmax = 1000;
%     dx = 10./(jmax-1);
%     x = 0:dx:10;
%     area = calcarea(x);
%     cfl = 0.9;
%     res_list = [];

%     if method == 1
%         %%%%%% generate initial conditions %%%%%%
%         % no shock
%         xsh = -1;
%         % use predictor-corrector method 
%         march_type = 1;
%         % use space marching
%         [rho_sp,u_sp,p_sp,e_sp,amach_sp] = spacemarch(gamma,fsmach,p0,rho0,xsh,x,area,march_type);
%       elseif method == 2
%         %%%%%% generate initial conditions %%%%%%
%         % shock
%         xsh = 4;
%         % use predictor-corrector method 
%         march_type = 1;
%         % use space marching
%         [rho_sp,u_sp,p_sp,e_sp,amach_sp] = spacemarch(gamma,fsmach,p0,rho0,xsh,x,area,march_type);
%       elseif method == 3
%         %%%%%% generate initial conditions %%%%%%
%         % shock
%         xsh = 4;
%         % use predictor-corrector method 
%         march_type = 1;
%         unsteady = 1;
%         % use space marching
%         [rho_sp,u_sp,p_sp,e_sp,amach_sp] = spacemarch(gamma,fsmach,p0,rho0,xsh,x,area,march_type);
%       end
%       p_end = p_sp(end);
%       Qh = convert_to_conservative(rho_sp, u_sp, e_sp, area);
%       max_iter = 4000;

%   %%%%%% Exp_starticit Time Marching %%%%%%
%   for i = 1:20
%     % compute primitive variables
%     Q = convert_to_primitive(Qh, gamma, area);

%     % compute speed of sound
%     c = sqrt(gamma .* Q(3,:) ./ Q(1,:));

%     % compute local timesteps
%     dt = compute_timestep(cfl, dx, Q(2,:), c);  

%     % compute fluxes
%     [Fhp, Fhm] = steger_warming_flux(Q, area, gamma);

%     % compute residuals 
%     res = compute_residual(Qh, Q, Fhp, Fhm, area, dx, dt, x);
%     res_list = [res_list, norm(res, 'fro')/sqrt(jmax)];
%     % app_starty boundary conditions at the exit using compatibility conditions
%     p_exit = p_end;
%     res = boundary_condition(Q, c, dt, dx, area, gamma, res, p_exit);
    
    
%     % update conservative variables
%     Qh = Qh + res;
    
%     % check if L2 of residual is 5 orders below the initial residual
%     % if i > 1 && norm(res, 'fro') < 1e-5 * norm(res_list(1), 'fro')
%     %     break;
%     % end
%   end
%   Q = convert_to_primitive(Qh, gamma, area);
%     rho = Q(1, :);
%     u = Q(2, :);
%     p = Q(3, :);
% end


function Qh = convert_to_conservative(rho, u, e, area)
    N = size(rho, 2);
    Qh = zeros(3, N);
    % convert to conservative variables
    Qh(1, :) = rho .* area;
    Qh(2, :) = rho .* u .* area;
    Qh(3, :) = e .* area;
  end
  
  function Q = convert_to_primitive(Qh, gamma, area)
    N = size(Qh, 2);
    Q = zeros(3, N);
    % convert to primitive variables
    Q(1, :) = Qh(1, :) ./ area;
    Q(2, :) = Qh(2, :) ./ Qh(1, :);
    Q(3, :) = (gamma - 1) .* (Qh(3,:) ./ area - 0.5 .* Q(1,:) .* Q(2,:).^2);
  end
  
  function dt = compute_timestep(cfl, dx, u, c)
    % compute local timesteps
    dt_local = cfl * dx ./ (abs(u) + abs(c));
    dt = min(dt_local, [], "all");
  end
  
  function residual = compute_residual(Qh, Q, Fhp, Fhm, area, dx, dt, x)
    % compute residuals
    N = size(Qh, 2);
    residual = zeros(3, N);
    for j = 2:size(Qh, 2)-1
      residual(:, j) = -dt / dx * (Fhp(:, j) - Fhp(:, j-1) + Fhm(:, j+1) - Fhm(:, j));
      residual(2, j) = residual(2,j)  + dt / dx * (Q(3,j) * abs((0.5*(area(j+1)+area(j)) - 0.5*(area(j)+area(j-1)))));
      % residual(2, j) = residual(2,j)  + dt / dx * (Q(3,j) * (calcarea(j*dx + 0.5 *dx) - calcarea(j*dx - 0.5 *dx)));
    end  
  end