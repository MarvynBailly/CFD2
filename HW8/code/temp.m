% temp.m
clear
close all
clc

method = 2;
animation = 0;
save_plots = 0;
implicit = 0;

fsmach = 1.265;   % Mach number at the entrance 
rho0 = 0.5;  % density at the entrance 
p0 = 0.379;   % pressure at the entrance 
gamma = 1.4;    % ratio of specific heats

% iteration_chop = []; % For 2.1, 3.1, 4.1, 4.4
early_stop = 1; % For 2.3, 2.4, 3.3, 3.4
% final_error_list = []; % For 2.4, 3.4

plot_iter = 50;
unsteady_iters = 500;
break_iter = 0;
unsteady_start_iter = 0;

cfl = 0.9; 
max_iter = 4000;
residual_history = zeros(max_iter, 1);

p_history = [];
rho_history = [];
iter_history = [];
t_history = [];

unsteady = 0;
converged = 0;

jmaxes = [181];

plot_p_vs_t = 0;

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
    xsh = 4;
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
    plot_p_vs_t = 1;
    amplitude = 0.03; % for graph 4.X
    unsteady = 1;
  elseif method == 3
    %%%%%% generate initial conditions %%%%%%
    % shock
    xsh = 4;
    % use predictor-corrector method 
    march_type = 1;
    unsteady = 1;
    % use space marching
    [rho_sp,u_sp,p_sp,e_sp,amach_sp] = spacemarch(gamma,fsmach,p0,rho0,xsh,x,area,march_type);
    amplitude = 0.03; % for graph 4.X
    plot_p_vs_t = 1; % for graph 4.2 and 4.3
  end


  % [rho_ex, u_ex, p_ex, e_ex] = exact_solution(x, gamma, fsmach, p0, rho0, 1, xsh);
  load('exact_solution.mat');


  p_end = p_sp(end);  
  % Qh = convert_to_conservative(rho0*rho_sp./rho_sp, u_sp(1)*u_sp./u_sp, e_sp(1)*e_sp./e_sp, area);
  Qh = convert_to_conservative(rho_sp, u_sp, e_sp, area);
  
  %%%%%% Explicit Time Marching %%%%%%
  % use Steger-Warming method to compute fluxes
  for i = 1:max_iter
    Q = convert_to_primitive(Qh, gamma, area);
    c = sqrt(gamma .* Q(3,:) ./ Q(1,:));

    % compute local timesteps
    if method == 3
        dt = 0.017695209471134; % Steady state CFL
    else
        dt = compute_timestep(cfl, dx, Q(2,:), c);
    end

    % apply boundary conditions at the exit using compatibility conditions
    if unsteady == 1 && converged == 1
      p0 = p_end;
      p_exit = p0 * (1 + amplitude * sin(i*2*pi/500));
    else 
      p_exit = p_end;
    end
    
    res_jmax = boundary_condition(Q, c, dt, dx, area, gamma, p_exit, implicit);

    [Fhp, Fhm] = steger_warming_flux(Q, area, gamma);
    if implicit
        % [Jp, Jm] = Flux_Jacobian(Qh, gamma);
        [Jp, Jm] = Flux_Jacobian_Spectral(Q, gamma);
        res = compute_residual_Im(Qh, Q, Fhp, Fhm, area, dx, dt, x, Jp, Jm, res_jmax);
    else
        res = compute_residual(Qh, Q, Fhp, Fhm, area, dx, dt, x);
    end

    res(:,end) = res_jmax;
    % err = [rho_sp;u_sp;p_sp] - Q;
    err = [rho_exact;u_exact;p_exact] - Q;
    
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

        figure(3)
        plot(x, Q(1,:), 'r-', 'LineWidth', 2); hold on
        plot(x, rho_sp); hold off
        title(['Density, t=',num2str(t)])
    end

    % if ismember(i, iteration_chop) % for graph 4.1 and 4.4
    %     figure(6)
    %     plot(x, Q(3,:), 'LineWidth', 1); hold on
    %     grid on
    % end

    % if plot_p_vs_t % for graph 4.2 and 4.3
    %     p_history = [p_history, Q(3,x==4.5)];
    %     t_history = [t_history, t-dt];
    % end

    if mod(i - unsteady_start_iter, plot_iter) == 0 && converged == 1
      p_history = [p_history, Q(3,:).'];
      rho_history = [rho_history, Q(1,:).'];
      t_history = [t_history, t-dt];
      iter_history = [iter_history, i];
    end

    % check if we should break
    if i == break_iter
      break
    end

    % check if L2 of residual is 5 orders below the initial residual
    if early_stop && res_list(end) < 1e-5 * res_list(1)
      disp(['residual is 5 orders below the initial residual by iteration ', num2str(i)]);
      % figure(8)
      % plot(x, Q(3,:), 'LineWidth', 1); hold on
      converged = 1;
      break_iter = i + unsteady_iters;
      unsteady_start_iter = i + 1;
      if unsteady == 1
        disp('beginning unsteady simulation')
      end
      % break
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
  plot(x, p_exact);
  plot(x, p_sp); hold off
  title('Pressure')
  legend('Numerical', 'Exact', 'Steady')
  grid on

  figure(4)
  semilogy(res_list, '-', 'LineWidth', 2); hold on
  semilogy(err_list, '-', 'LineWidth', 2);
  legend('residual', 'error')
  grid on

  % if plot_p_vs_t % for graph 4.2 and 4.3
  if implicit == 1
    name = 'Implicit';
  else 
    name = 'Explicit';
  end
  figure(7);
  hold on;
  legend_entries = cell(size(t_history));
  for j = 1:length(t_history)
      plot(x(:), p_history(:, j), 'LineWidth', 1.5);
      legend_entries{j} = sprintf('iter = %d', round(iter_history(j)));
  end
  xlabel('x'); ylabel('Pressure');
  title(['Pressure Distribution Snapshots ', name]);
  legend(legend_entries, 'Location', 'best');
  grid on;


  figure(8);
  hold on;
  legend_entries = cell(size(t_history));
  for j = 1:length(t_history)
      plot(x(:), rho_history(:, j), 'LineWidth', 1.5);
      legend_entries{j} = sprintf('iter = %d', round(iter_history(j)));
  end
  xlabel('x'); ylabel('Density');
  title(['Density Distribution Snapshots ', name]);
  legend(legend_entries, 'Location', 'best');
  grid on;




  % final_error_list = [final_error_list, err_list(end)];
end

% if early_stop
%   figure(8)
%   legends = num2str(jmaxes');
%   legends = [legends; "Exact"];
%   plot(x, p_sp, 'k-', 'LineWidth', 1); hold on
%   legend(legends, Location="best")
%   xlabel('x')
%   ylabel('p')
%   grid on
% end

% figure(9)
% loglog(jmaxes, final_error_list, 'o-', 'LineWidth', 2);
% ylabel('Final error')
% xlabel('JMAX')
% grid on


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

function [Jp, Jm] = Flux_Jacobian(Qh, gamma)
    N = size(Qh,2);
    J = zeros(N,3,3);
    
    J(:,1,2) = 1;
    J(:,2,1) = (gamma-3)/2 * (Qh(2,:)./Qh(1,:)).^2;
    J(:,2,2) = (3-gamma) * (Qh(2,:)./Qh(1,:));
    J(:,2,3) = gamma-1;
    J(:,3,1) = (gamma-1)*(Qh(2,:)./Qh(1,:)).^3 - gamma.*Qh(3,:).*Qh(2,:)./Qh(1,:).^2;
    J(:,3,2) = gamma*Qh(3,:)./Qh(1,:) - 3/2*(gamma-1)*(Qh(2,:)./Qh(1,:)).^2;
    J(:,3,3) = gamma*Qh(2,:)./Qh(1,:);

    Jp = 0.5 * (J+abs(J)); 
    Jm = 0.5 * (J-abs(J));
end



function [Jp, Jm] = Flux_Jacobian_(Q, gamma)
    N = size(Q,2);
    J = zeros(N,3,3);
    
    e = Q(3,:)/(gamma-1) + 0.5*Q(1,:).*Q(2,:).^2;
    J(:,1,2) = 1;
    J(:,2,1) = (gamma-3)/2 * Q(2,:).^2;
    J(:,2,2) = (3-gamma) * Q(2,:);
    J(:,2,3) = gamma-1;
    J(:,3,1) = (gamma-1)*Q(2,:).^3 - gamma*e.*Q(2,:)./Q(1,:);
    J(:,3,2) = gamma*e.*Q(1,:) - 3/2*(gamma-1)*Q(2,:).^2;
    J(:,3,3) = gamma*Q(2,:);

    Jp = 0.5 * (J+abs(J));
    Jm = 0.5 * (J-abs(J));
end


function [Jp, Jm] = Flux_Jacobian_Spectral(Q, gamma)
  rho = Q(1,:);
  u   = Q(2,:);
  p   = Q(3,:);
  e   = p./(gamma - 1) + 0.5 * rho .* u.^2;

  m = length(rho);
  J = zeros(m, 3, 3);  % Full Jacobian at each point

  for j = 1:m
      % Fill the 3x3 Jacobian matrix at point j
      J(j,:,:) = [
          0,                         1,                       0;
          0.5*(gamma-3)*u(j)^2,      (3-gamma)*u(j),          gamma - 1;
          (gamma-1)*u(j)^3 - gamma*e(j)*u(j)/rho(j), ...
          gamma*e(j)/rho(j) - 1.5*(gamma-1)*u(j)^2, ...
          gamma*u(j)
      ];
  end

  % Spectral radius σ = |u| + a
  a = sqrt(gamma * p ./ rho);
  sigma = abs(u) + a;

  % Broadcast σ * I at each point
  Jp = 0.5 * (J + reshape(sigma, [m 1 1]) .* eye3(m));
  Jm = 0.5 * (J - reshape(sigma, [m 1 1]) .* eye3(m));
end

function I3 = eye3(N)
  % Returns N stacked 3×3 identity matrices: size N×3×3
  I = eye(3);
  I3 = repmat(reshape(I, 1, 3, 3), N, 1, 1);
end