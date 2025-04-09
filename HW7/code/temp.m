% temp.m
clear
close all
clc

method = 1;
save_plots = 0;

fsmach = 1.265;   % Mach number at the entrance 
rho0 = 0.5;  % density at the entrance 
p0 = 0.379;   % pressure at the entrance 
gamma = 1.4;    % ratio of specific heats


cfl = 0.9; 
max_iter = 2000;
residual_history = zeros(max_iter, 1);

jmaxes = [61];

for jmax = jmaxes
  dx = 10./(jmax-1);
  x = 0:dx:10;
  area = calcarea(x);

  % no shock
  if method == 1
    %%%%%% generate initial conditions %%%%%%
    % no shock
    xsh = -1;
    % use predictor-corrector method 
    march_type = 1;
    % use space marching
    [rho_sp,u_sp,p_sp,e_sp,amach_sp] = spacemarch(gamma,fsmach,p0,rho0,xsh,x,area,march_type);
  end

  % plot initial conditions
  figure(1)
  subplot(2,1,1)
  plot(x, rho_sp, 'r-', 'LineWidth', 2)
  title('Density')
  xlabel('x')
  ylabel('Density')
  grid on
  subplot(2,1,2)
  plot(x, p_sp, 'r-', 'LineWidth', 2)
  title('Pressure')
  xlabel('x')
  ylabel('Pressure')
  grid on

  % form conservative variables
  Qh = zeros(3, jmax);
  Fh = zeros(3, jmax);

  for j = 1:jmax
    Fh(1, j) = rho_sp(j) * u_sp(j) * area(j);
    Fh(2, j) = (rho_sp(j) * u_sp(j)^2 + p_sp(j)) * area(j);
    Fh(3, j) = (e_sp(j) + p_sp(j)) * u_sp(j) * area(j);

    Qh(1, j) = rho_sp(j) * area(j);
    Qh(2, j) = rho_sp(j) * u_sp(j) * area(j);
    Qh(3, j) = e_sp(j)* area(j);
  end


  %%%%%% Explicit Time Marching %%%%%%
  % use Steger-Warming method to compute fluxes
  res = zeros(3, jmax);
  for t = 1:20
    % recover primitive variables
    rho = Qh(1, :) ./ area;
    rhou = Qh(2, :) ./ area;        
    e = Qh(3, :) ./ area;
    
    u = rhou ./ rho;
    u2 = u.^2;
    p = (gamma - 1) .* (e - 0.5 * rho .* u2);

    % compute speed of sound and enthalpy
    c = sqrt(gamma .* p ./ rho);

    % compute local timesteps
    dt_local = cfl * dx ./ (abs(u) + abs(c));
    dt = min(dt_local, [], "all");
    

    % compute fluxes
    [Fhp, Fhm] = steger_warming_flux(Qh, area, gamma);

    % FIX BOUNDARY CONDITIONS
    % compute residuals 
    for j = 2:jmax-1
      res(:, j) = -dt / dx * (Fhp(:, j) - Fhp(:, j-1) + Fhm(:, j+1) - Fhm(:, j));
      res(2, j) = res(2,j) + dt / dx * (p(j) * (0.5*(area(j+1)-area(j)) - 0.5*(area(j)-area(j-1))));
    end
  
    % fix to inlet?
    Qh(:,1) = [rho_sp(1); rho_sp(1)*u_sp(1); e_sp(1)] * area(1);
    % extrap for now
    Qh(:,jmax) = Qh(:,jmax-1);

    % update conservative variables
    Qh = Qh + res;
  end
end


% plot results
% recover primitive variables
rho = Qh(1, :) ./ area;
rhou = Qh(2, :) ./ area;        
e = Qh(3, :) ./ area;

u = rhou ./ rho;
u2 = u.^2;
p = (gamma - 1) .* (e - 0.5 * rho .* u2);
size(p)

% plot pressure
figure(2)
plot(x, p, 'r-', 'LineWidth', 2)
