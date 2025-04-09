% temp.m
clear
close all
clc

method = 3;
save_plots = 1;

fsmach = 1.265;   % Mach number at the entrance 
rho0 = 0.5;  % density at the entrance 
p0 = 0.379;   % pressure at the entrance 
gamma = 1.4;    % ratio of specific heats

jmaxes = [101, 51];


figure(1);  clf; hold on;
figure(2);  clf; hold on;
figure(3);  clf; hold on;
figure(4);  clf; hold on;

ph1 = {}; ph2 = {}; ph3 = {}; ph4 = {}; % handles for figures 1, 2, 3
legend_entries = {};


for jmax = jmaxes
  dx = 10./(jmax-1);
  x = 0:dx:10;
  area = calcarea(x);

  if jmax == 101
    style = '-';  % dotted line
  else
      style = '--';   % solid line
  end

  % no shock
  if method == 1
    % first order accurate space marching
    xsh = -1;
    march_type = 0;
    [rho_sp,u_sp,p_sp,e_sp,amach_sp] = spacemarch(gamma,fsmach,p0,rho0,xsh,x,area,march_type);
  end

  % no shock
  if method == 2
    % predict and correcter accurate space marching
    xsh = -1;
    march_type = 1;
    [rho_sp,u_sp,p_sp,e_sp,amach_sp] = spacemarch(gamma,fsmach,p0,rho0,xsh,x,area,march_type);
  end

  if method == 3
    % fixed shock and predictor-corrector
    xsh = 4;
    march_type = 1;
    [rho_sp,u_sp,p_sp,e_sp,amach_sp] = spacemarch(gamma,fsmach,p0,rho0,xsh,x,area,march_type);
  end

  figure(1); 
  h1 = plot(x, p_sp, style, 'LineWidth', 2.0);
  ph1{end+1} = h1;

  % Plot and store handle
  figure(2); 
  h2 = plot(x, u_sp, style, 'LineWidth', 2.0);
  ph2{end+1} = h2;

  % Plot and store handle
  figure(3); 
  h3 = plot(x, rho_sp, style, 'LineWidth', 2.0);
  ph3{end+1} = h3;

  figure(4); 
  h4 = plot(x, amach_sp, style, 'LineWidth', 2.0);
  ph4{end+1} = h4;

  legend_entries{end+1} = sprintf('jmax = %d', jmax);
end


switch method
  case 1
      method_name = 'no_shock-fist_order';
  case 2
      method_name = 'no_shock-p_c';
  case 3
      method_name = 'shock-p_c';
end

figure(1);
set(gca, 'FontSize', 16, 'LineWidth', 2.0, 'FontWeight', 'demi');
% vertically strecth the figure
title('Pressure Distribution Along Nozzle');
xlabel('X'); ylabel('p');
grid on;
legend([ph1{:}], legend_entries, 'Location', 'best');
if save_plots == 1 print('-dpng', fullfile('images', [method_name '-pressure.png'])); end

figure(2);
set(gca,'FontSize',16,'LineWidth',2.0,'FontWeight','demi');
title(['Velocity Distribution Along Nozzle'])
grid on;
xlabel(['X']); ylabel(['u']);
legend([ph2{:}], legend_entries, 'Location', 'best');
if save_plots == 1 print('-dpng', fullfile('images', [method_name '-velocity.png'])); end

figure(3);
set(gca,'FontSize',16,'LineWidth',2.0,'FontWeight','demi');
title(['Density Distribution Along Nozzle'])
grid on;
xlabel(['X']); ylabel(['rho']);
legend([ph3{:}], legend_entries, 'Location', 'best');
if save_plots == 1 print('-dpng', fullfile('images', [method_name '-density.png'])); end

figure(4);
set(gca,'FontSize',16,'LineWidth',2.0,'FontWeight','demi');
title(['Mach Number Distribution Along Nozzle'])
grid on;
xlabel(['X']); ylabel(['mach number']);
legend([ph4{:}], legend_entries, 'Location', 'best');
if save_plots == 1 print('-dpng', fullfile('images', [method_name '-mach.png'])); end