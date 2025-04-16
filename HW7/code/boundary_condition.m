function res = boundary_condition(Q, c, dt, dx, area, gamma, res, p_end)
 % apply boundary conditions at the exit using compatibility conditions
 % compute R1, R2, and R3
    u_jmax = Q(2, end);

    R1 = - u_jmax * (dt/dx) * ((Q(1, end) - Q(1, end-1)) - (1/c(end)^2) *(Q(3, end) - Q(3, end-1)));
    %    - (1/c(end)^2) * (gamma * Q(3, end) * Q(1, end) / area(end)) * - (Q(1, end) * Q(2, end) / area(end));

    R2 = - (dt/dx) * ((u_jmax + c(end)) * (Q(3, end) - Q(3, end-1) + (Q(1,end) * c(end)) * (Q(2, end) - Q(2, end-1)))...
     - Q(1,end) * Q(2,end) * c(end)^2 / area(end) * (area(end) - area(end-1)));
    
    R3 = - (dt/dx) * ((u_jmax - c(end)) * (Q(3, end) - Q(3, end-1) - (Q(1,end) * c(end)) * (Q(2, end) - Q(2, end-1)))...
     - Q(1,end) * Q(2,end) * c(end)^2 / area(end) * (area(end) - area(end-1)));

    if (Q(2, end-1) >= c(end-1))
        % disp('supersonic')
        dp = 0.5 * (R2 + R3);
    else
        % change this to sin once converged ?
        dp = p_end - Q(3, end);
    end

    drho = R1 + dp / c(end)^2;
    du = (R2 - dp) / (Q(1, end) * c(end));

    % apply to conserved variables
   res(1, end) =  area(end) * drho;
   res(2, end) =  area(end) * ((Q(1, end) + drho) * (Q(2, end) + du) - Q(1, end) * Q(2, end));
   res(3, end) =  area(end) * ( (Q(3, end) + dp)/(gamma -1) + 0.5 * (Q(1, end) + drho) * (Q(2, end) + du)^2 - Q(3, end)/(gamma -1) - 0.5 * Q(1, end) * Q(2, end)^2);  
end