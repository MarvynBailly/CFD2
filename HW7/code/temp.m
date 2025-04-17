% temp.m
clear
close all
clc

method = 3;
save_plots = 0;

fsmach = 1.265;   % Mach number at the entrance 
rho0 = 0.5;  % density at the entrance 
p0 = 0.379;   % pressure at the entrance 
gamma = 1.4;    % ratio of specific heats


cfl = 0.9; 
max_iter = 4000;
residual_history = zeros(max_iter, 1);

unsteady = 0;
converged = 0;

jmaxes = [200];

for jmax = jmaxes
  dx = 10./(jmax-1);
  x = 0:dx:10;
  area = calcarea(x);
  t = 0;
  res_list = [];
  % error = []; Need exact solution?

  % no shock
  if method == 1
    %%%%%% generate initial conditions %%%%%%
    % no shock
    xsh = -1;
    % use predictor-corrector method 
    march_type = 1;
    % use space marching
    [rho_sp,u_sp,p_sp,e_sp,amach_sp] = spacemarch(gamma,fsmach,p0,rho0,xsh,x,area,march_type);
  elseif method == 2
    %%%%%% generate initial conditions %%%%%%
    % shock
    xsh = 4;
    % use predictor-corrector method 
    march_type = 1;
    % use space marching
    [rho_sp,u_sp,p_sp,e_sp,amach_sp] = spacemarch(gamma,fsmach,p0,rho0,xsh,x,area,march_type);
  elseif method == 3
    %%%%%% generate initial conditions %%%%%%
    % shock
    xsh = 4;
    % use predictor-corrector method 
    march_type = 1;
    unsteady = 1;
    % use space marching
    [rho_sp,u_sp,p_sp,e_sp,amach_sp] = spacemarch(gamma,fsmach,p0,rho0,xsh,x,area,march_type);
  end


  p_end = p_sp(end);  
  % Qh = convert_to_conservative(rho0*rho_sp./rho_sp, u_sp(1)*u_sp./u_sp, e_sp(1)*e_sp./e_sp, area);
  Qh = convert_to_conservative(rho_sp, u_sp, e_sp, area);
  
  %%%%%% Explicit Time Marching %%%%%%
  % use Steger-Warming method to compute fluxes
  for i = 1:max_iter
    % compute primitive variables
    Q = convert_to_primitive(Qh, gamma, area);

    % compute speed of sound
    c = sqrt(gamma .* Q(3,:) ./ Q(1,:));

    % compute local timesteps
    dt = compute_timestep(cfl, dx, Q(2,:), c);  

    % compute fluxes
    [Fhp, Fhm] = steger_warming_flux(Q, area, gamma);

    % compute residuals 
    res = compute_residual(Qh, Q, Fhp, Fhm, area, dx, dt, x);

    % apply boundary conditions at the exit using compatibility conditions
    if unsteady == 1 && converged == 1
      p0 = p_end;
      amp = 0.3;
      T = 500 * dt;
      p_exit = p0 * (1 + amp * sin(2 * pi * t / T));
    else 
      p_exit = p_end;
    end

    res = boundary_condition(Q, c, dt, dx, area, gamma, res, p_exit);

    
    % update conservative variables
    Qh = Qh + res;
    t = t + dt;
    res_list = [res_list, norm(res, 'fro')];
    
    % check if L2 of residual is 5 orders below the initial residual
    if i > 1 && norm(res, 'fro') < 1e-3 * norm(res_list(1), 'fro')
      converged = 1;
      max_iter = max_iter + 2000;
      if unsteady == 1
        disp('beginning unsteady simulation')
      end
    end
    figure(2)
    plot(x, Q(3,:), 'r-', 'LineWidth', 2); hold on
    plot(x, p_sp); hold off
    title(['Pressure, t=',num2str(t)])

    % pause(0.1)
  end

  % plot things
  Q = convert_to_primitive(Qh, gamma, area);

  figure(1)
  plot(x, Q(1,:), 'r-', 'LineWidth', 2); hold on
  plot(x, rho_sp); hold off
  title('Density')
  grid on

  figure(2)
  plot(x, Q(2,:), 'r-', 'LineWidth', 2); hold on
  plot(x, u_sp); hold off
  title('Velocity')
  grid on

  figure(3)
  plot(x, Q(3,:), 'r-', 'LineWidth', 2); hold on
  plot(x, p_sp); hold off
  title('Pressure')
  grid on

  figure(4)
  plot(res_list, 'r-', 'LineWidth', 2)
  grid on

end




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


