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

    %%%%%% Explicit Time Marching %%%%%%
    % use Steger-Warming method

  end
end