function res = boundary_condition(Q, c, dt, dx, area, gamma, res, p_end)
 % apply boundary conditions at the exit using compatibility conditions
 % seems like the only difference must be in here???
 % compute R1, R2, and R3
    u_jmax = Q(2, end);
    u_jmaxm1 = Q(2, end-1);
    rho_jmax = Q(1, end);
    rho_jmaxm1 = Q(1, end-1);
    p_jmax = Q(3, end);
    p_jmaxm1 = Q(3, end-1);

    % R1 = - u_jmax * (dt/dx) * ((rho_jmax - rho_jmaxm1) - (1/c(end)^2) * (p_jmax - p_jmaxm1));
    % R2 = - (dt/dx) * ((u_jmax + c(end)) * (p_jmax - p_jmaxm1) + (rho_jmax * c(end)) * (u_jmax - u_jmaxm1) ...
    %     - (rho_jmax * u_jmax * c(end)^2 / area(end)) * (area(end) - area(end-1)));
    % R3 = - (dt/dx) * ((u_jmax - c(end)) * (p_jmax - p_jmaxm1) - (rho_jmax * c(end)) * (u_jmax - u_jmaxm1) ...
    %     - (rho_jmax * u_jmax * c(end)^2 / area(end)) * (area(end) - area(end-1)));

    u = Q(2, :);
    rho = Q(1, :);
    p = Q(3, :);
    dtdx = dt/dx;
    jmax = length(area);
    j=jmax;
    speed = sqrt(gamma*p(j)/rho(j));
    speedm = sqrt(gamma*p(j-1)/rho(j-1));
    alpha1 = 0.5*(u(j-1)+u(j))*dtdx;
    alpha2 = 0.5*(u(j-1)+u(j)+speedm+speed)*dtdx;
    alpha3 = 0.5*(u(j-1)+u(j)-speedm-speed)*dtdx;
    sp2inv = 0.5*(1./speed^2+1./speedm^2);
    rhospeed = 0.5*(rho(j)*speed+rho(j-1)*speedm);
    acomp = 0.5*(rho(j)*speed*u(j)*speed/area(j) + ...
                rho(j-1)*speedm*u(j-1)*speedm/area(j-1));
    r1 = -alpha1*((rho(j)-rho(j-1))-sp2inv*(p(j)-p(j-1)));
    r2 = -alpha2*((p(j)-p(j-1))+rhospeed*(u(j)-u(j-1))) ...
                 -acomp*(area(j)-area(j-1))*dtdx;
    r3 = -alpha3*((p(j)-p(j-1))-rhospeed*(u(j)-u(j-1))) ...
                 -acomp*(area(j)-area(j-1))*dtdx;


    if (u(j-1) >= speedm)   % supersonic exit
        dp = 0.5*(r2+r3);
    else
        dp = p_end - p_jmax;
    end
    drho = r1+dp*sp2inv;
    du = (r2-dp)/(rhospeed);
    res(1,j) = area(j)*drho;
    res(2,j) = area(j)*( (rho(j)+drho)*(u(j)+du) - rho(j)*u(j) );
    res(3,j) = area(j)*( ( (p(j)+dp)/(gamma-1)+0.5*(rho(j)+drho)*(u(j)+du)^2 ) ...
                        -( p(j)/(gamma-1)+0.5*rho(j)*u(j)^2 ) );




%     if (u_jmaxm1 >= c(end-1))
%         dp = 0.5 * (R2 + R3);
%     else
%         % dp = p_end - p_jmax;
%         dp = 0;
%     end

%     drho = R1 + dp / c(end)^2;
%     du = (R2 - dp) / (rho_jmax * c(end));

%     % apply to conserved variables
%    res(1, end) =  area(end) * drho;
%    res(2, end) =  area(end) * ((rho_jmax + drho) * (u_jmax + du) - rho_jmax * u_jmax);
%    res(3, end) =  area(end) * ( (p_jmax + dp)/(gamma -1) + 0.5 * (rho_jmax + drho) * (u_jmax + du)^2 - p_jmax/(gamma -1) - 0.5 * rho_jmax * u_jmax^2);  
end