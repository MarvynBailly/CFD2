function res = boundary_condition(Q, c, dt, dx, area, gamma, res, p_end)
 % apply boundary conditions at the exit using compatibility conditions
 % compute R1, R2, and R3
    u_jmax = Q(2, end);
    u_jmaxm1 = Q(2, end-1);
    rho_jmax = Q(1, end);
    rho_jmaxm1 = Q(1, end-1);
    p_jmax = Q(3, end);
    p_jmaxm1 = Q(3, end-1);

    R1 = - u_jmax * (dt/dx) * ((rho_jmax - rho_jmaxm1) - (1/c(end)^2) * (p_jmax - p_jmaxm1));
    R2 = - (dt/dx) * ((u_jmax + c(end)) * (p_jmax - p_jmaxm1) + (rho_jmax * c(end)) * (u_jmax - u_jmaxm1) ...
        - (rho_jmax * u_jmax * c(end)^2 / area(end)) * (area(end) - area(end-1)));
    R3 = - (dt/dx) * ((u_jmax - c(end)) * (p_jmax - p_jmaxm1) - (rho_jmax * c(end)) * (u_jmax - u_jmaxm1) ...
        - (rho_jmax * u_jmax * c(end)^2 / area(end)) * (area(end) - area(end-1)));

    % R1 = - u_jmax * (dt/dx) * ((Q(1, end) - Q(1, end-1)) - (1/c(end)^2) *(Q(3, end) - Q(3, end-1)));
    %    - (1/c(end)^2) * (gamma * Q(3, end) * Q(1, end) / area(end)) * - (Q(1, end) * Q(2, end) / area(end));

    % R2 = - (dt/dx) * ((u_jmax + c(end)) * (Q(3, end) - Q(3, end-1) + (Q(1,end) * c(end)) * (Q(2, end) - Q(2, end-1)))...
    %  - Q(1,end) * Q(2,end) * c(end)^2 / area(end) * (area(end) - area(end-1)));
    
    % R3 = - (dt/dx) * ((u_jmax - c(end)) * (Q(3, end) - Q(3, end-1) - (Q(1,end) * c(end)) * (Q(2, end) - Q(2, end-1)))...
    %  - Q(1,end) * Q(2,end) * c(end)^2 / area(end) * (area(end) - area(end-1)));

    if (u_jmaxm1 >= c(end-1))
        disp('supersonic')
        dp = 0.5 * (R2 + R3);
    else
        disp('subsonic')
        % change this to sin once converged ?
        dp = p_end - p_jmax;
    end

    drho = R1 + dp / c(end)^2;
    du = (R2 - dp) / (rho_jmax * c(end));

    % apply to conserved variables
   res(1, end) =  area(end) * drho;
   res(2, end) =  area(end) * ((rho_jmax + drho) * (u_jmax + du) - rho_jmax * u_jmax);
   res(3, end) =  area(end) * ( (p_jmax + dp)/(gamma -1) + 0.5 * (rho_jmax + drho) * (u_jmax + du)^2 - p_jmax/(gamma -1) - 0.5 * rho_jmax * u_jmax^2);  
end