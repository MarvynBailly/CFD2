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
    % e = p/(gamma-1) + 0.5 * rho .* u.^2;
    
end