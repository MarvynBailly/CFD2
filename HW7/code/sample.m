% sample.m
% Just a few simple calls to the different functions for 
% finding the steady flow using space marching
% Uses: calcarea.m; 
%       spacemarch.m (and hence march.m and shock.m)
%       findshock.m (and hence spacemarch.m, etc.)
clear;
for ifig=1:2
figure(ifig); clf;
end
% set-up the mesh with x-location and area
jmax=61;
dx = 1./(jmax-1);
x = 0:dx:1;
area = calcarea(x);
% Entrance conditions
amach0 = 2.30;  % Mach number at the entrance 
rho0 = 0.2500;  % density at the entrance 
p0 = 0.07000;   % pressure at the entrance 
gamma = 1.4;    % ratio of specific heats
% Target exit pressure
pexit=1./gamma;
% Let's plot pressure distribution for fixed shock location
figure(1);
hold on;
xsh = 0.6;
[rho_sp,u_sp,p_sp,e_sp,amach_sp]=...
            spacemarch(gamma,amach0,p0,rho0,xsh,x,area,1);
set(gca,'FontSize',14,'LineWidth',2.0,'FontWeight','demi');
plot(x,p_sp,'go');
%end
%set(gca,'FontSize',[16],'LineWidth',[2.0],'FontWeight','demi');
title(['Pressure Distribution Along Nozzle (fixed shock locations)'])
grid on;
xlabel(['X']); ylabel(['p']);
% Let's plot pressure distribution for target pressure at exit using different
% space marching methods (predictor-corrector vs. predictor only)
figure(2);
xsh1=0.65;
xsh2=0.85;
% Using predictor-corrector find shock location and plot
[xsh]= findshock(gamma,amach0,p0,rho0,xsh1,xsh2,x,area,1,pexit);
[rho_sp,u_sp,p_sp,e_sp,amach_sp]=...
            spacemarch(gamma,amach0,p0,rho0,xsh,x,area,1);
plot(x,p_sp,'r-','LineWidth',2);
%set(gca,'FontSize',[16],'LineWidth',[2.0],'FontWeight','demi');
title(['Pressure Distribution Along Nozzle (fixed exit pressure)'])
grid on;
hold on
% Using predictor only find shock location and plot
[xsh]= findshock(gamma,amach0,p0,rho0,xsh1,xsh2,x,area,0,pexit);
[rho_sp,u_sp,p_sp,e_sp,amach_sp]=...
            spacemarch(gamma,amach0,p0,rho0,xsh,x,area,0);
set(gca,'FontSize',14,'LineWidth',2.0,'FontWeight','demi');
plot(x,p_sp,'b--','LineWidth',2);
xlabel(['X']); ylabel(['p']);
% Done
