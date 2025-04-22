% temp.m
clear
close all
clc

method = 1;
animation = 0;
save_plots = 0;

fsmach = 1.265;   % Mach number at the entrance 
rho0 = 0.5;  % density at the entrance 
p0 = 0.379;   % pressure at the entrance 
gamma = 1.4;    % ratio of specific heats

iteration_chop = []; % For 2.1, 3.1, 4.1, 4.4
early_stop = 1; % For 2.3, 2.4, 3.3, 3.4
final_error_list = []; % For 2.4, 3.4

cfl = 0.9; 
max_iter = 4000;
residual_history = zeros(max_iter, 1);

p_history = [];
t_history = [];

unsteady = 0;
converged = 0;

jmaxes = [61];

for jmax = jmaxes
  dx = 10./(jmax-1);
  x = 0:dx:10;
  area = calcarea(x);
  t = 0;
  res_list = [];
  err_list = [];

  % no shock
  if method == 1
    %%%%%% generate initial conditions %%%%%%
    % no shock
    xsh = -1;
    % use predictor-corrector method 
    march_type = 1;
    % use space marching
    [rho_sp,u_sp,p_sp,e_sp,amach_sp] = spacemarch(gamma,fsmach,p0,rho0,xsh,x,area,march_type);
    plot_p_vs_t = 0;
  elseif method == 2
    %%%%%% generate initial conditions %%%%%%
    % shock
    xsh = 4;
    % use predictor-corrector method 
    march_type = 1;
    % use space marching
    [rho_sp,u_sp,p_sp,e_sp,amach_sp] = spacemarch(gamma,fsmach,p0,rho0,xsh,x,area,march_type);
    plot_p_vs_t = 0;
  elseif method == 3
    %%%%%% generate initial conditions %%%%%%
    % shock
    xsh = 4;
    % use predictor-corrector method 
    march_type = 1;
    unsteady = 1;
    % use space marching
    [rho_sp,u_sp,p_sp,e_sp,amach_sp] = spacemarch(gamma,fsmach,p0,rho0,xsh,x,area,march_type);
    amplitude = 0.01; % for graph 4.X
    plot_p_vs_t = 1; % for graph 4.2 and 4.3
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
    if method == 3
        dt = 0.017695209471134; % Steady state CFL
    else
        dt = compute_timestep(cfl, dx, Q(2,:), c);
    end

    % compute fluxes
    [Fhp, Fhm] = steger_warming_flux(Q, area, gamma);

    % compute residuals 
    res = compute_residual(Qh, Q, Fhp, Fhm, area, dx, dt, x);

    % apply boundary conditions at the exit using compatibility conditions
    if unsteady == 1 %&& converged == 1
      p0 = p_end;
      p_exit = p0 * (1 + amplitude * sin(i*2*pi/500));
    else 
      p_exit = p_end;
    end

    res = boundary_condition(Q, c, dt, dx, area, gamma, res, p_exit);
    err = [rho_sp;u_sp;p_sp] - Q;
    
    % update conservative variables
    Qh = Qh + res;
    t = t + dt;
    res_list = [res_list, norm(res, 'fro')/sqrt(jmax)];
    err_list = [err_list, norm(err, 'fro')/sqrt(jmax)];
    

    if animation
        figure(2)
        plot(x, Q(3,:), 'r-', 'LineWidth', 2); hold on
        plot(x, p_sp); hold off
        title(['Pressure, t=',num2str(t)])
    end

    if ismember(i, iteration_chop) % for graph 4.1 and 4.4
        figure(6)
        plot(x, Q(3,:), 'LineWidth', 1); hold on
        grid on
    end

    if plot_p_vs_t % for graph 4.2 and 4.3
        p_history = [p_history, Q(3,x==4.5)];
        t_history = [t_history, t-dt];
    end

    % check if L2 of residual is 5 orders below the initial residual
    if early_stop && res_list(end) < 1e-5 * res_list(1)
      disp(['residual is 5 orders below the initial residual by iteration ', num2str(i)]);
      figure(8)
      plot(x, Q(3,:), 'LineWidth', 1); hold on
      %converged = 1;
      %max_iter = max_iter + 2000;
      %if unsteady == 1
      %  disp('beginning unsteady simulation')
      %end
      break
    end
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
  semilogy(res_list, '-', 'LineWidth', 2); hold on
  semilogy(err_list, '-', 'LineWidth', 2);
  legend('residual', 'error')
  grid on

  figure(6)
  legends = num2str(iteration_chop');
  if method==1 || method==2
      plot(x, p_sp, 'k-', 'LineWidth', 1); hold on
      legends = [legends; "Exact"];
  end
  legend(legends, Location="best")
  xlabel('x')
  ylabel('p')

  if plot_p_vs_t % for graph 4.2 and 4.3
     figure(7)
     plot(t_history, p_history, 'r-', 'LineWidth', 2)
     xlabel('t')
     ylabel('p')
     grid on
  end

  final_error_list = [final_error_list, err_list(end)];
end

if early_stop
  figure(8)
  legends = num2str(jmaxes');
  legends = [legends; "Exact"];
  plot(x, p_sp, 'k-', 'LineWidth', 1); hold on
  legend(legends, Location="best")
  xlabel('x')
  ylabel('p')
  grid on
end

figure(9)
loglog(jmaxes, final_error_list, 'o-', 'LineWidth', 2);
ylabel('Final error')
xlabel('JMAX')
grid on


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


