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

  % form conservative variables
  Qh = zeros(3, jmax);
  Q = zeros(3, jmax);

  Qh = convert_to_conservative(Qh, rho_sp, u_sp, e_sp, area);

  %%%%%% Explicit Time Marching %%%%%%
  % use Steger-Warming method to compute fluxes
  for t = 1:10
    % compute primitive variables
    Q = convert_to_primitive(Q, Qh, gamma, area);

    % compute speed of sound and enthalpy
    c = sqrt(gamma .* Q(3,:) ./ Q(1,:));

    % compute local timesteps
    dt = compute_timestep(cfl, dx, Q(2,:), c);    

    % compute fluxes
    [Fhp, Fhm] = steger_warming_flux(Q, area, gamma);

    % FIX BOUNDARY CONDITIONS
    % compute residuals 
    res = compute_residual(Qh, Q, Fhp, Fhm, area, dx, dt);
  
    % fix to inlet?
    Qh(:,1) = [rho_sp(1); rho_sp(1)*u_sp(1); e_sp(1)] * area(1);
    
    % extrap for now
    Qh(:,jmax) = Qh(:,jmax-1);

    % update conservative variables
    Qh = Qh + res;
  end
end

function Qh = convert_to_conservative(Qh, rho, u, e, area)
  % convert to conservative variables
  Qh(1, :) = rho .* area;
  Qh(2, :) = rho .* u .* area;
  Qh(3, :) = e .* area;
end

function Q = convert_to_primitive(Q, Qh, gamma, area)
  % convert to primitive variables
  Q(1, :) = Qh(1, :) ./ area;
  Q(2, :) = Qh(2, :) ./ Qh(2, :);
  Q(3, :) = (gamma - 1) .* (Qh(3,:) ./ area - 0.5 .* Q(1,:) .* Q(2,:).^2);
end

function dt = compute_timestep(cfl, dx, u, c)
  % compute local timesteps
  dt_local = cfl * dx ./ (abs(u) + abs(c));
  dt = min(dt_local, [], "all");
end

function residual = compute_residual(Qh, Q, Fhp, Fhm, area, dx, dt)
  % compute residuals
  N = size(Qh, 2);
  residual = zeros(3, N);
  for j = 2:size(Qh, 2)-1
    residual(:, j) = -dt / dx * (Fhp(:, j) - Fhp(:, j-1) + Fhm(:, j+1) - Fhm(:, j));
    residual(2, j) = residual(2,j)  + dt / dx * (Q(3,j) * (0.5*(area(j+1)-area(j)) - 0.5*(area(j)-area(j-1))));
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

% plot pressure
figure(2)
plot(x, p, 'r-', 'LineWidth', 2)
