% quasi1d.m
%______________________________________________________________________
% This Is a Program for ENAE464
% Program by DR JAMES D BAEDER
%
% Designed to study the quasi-1D flow in a nozzle w/ Euler equations
%
% Notes:
%  1. The area distribution along the nozzle is defined in calcarea.m
%  2. Entrance of the nozzle is assumed to be supersonic
%     (Set the entrance Mach number, pressure and density)
%  3. Non-dimensionalization of quantities is assumed
%  4. Shock location is calculated to give the required exit pressure
%  5. Initial guess for steady state is calculated using space marching
%     with a predictor-corrector method
%  6. The unsteady equations are solved using flux-vector splitting (Steger-Warming)
%  7. The method is explicit and first order in time and space
%  8. The source term is calculated using the interface areas
%______________________________________________________________________
%
% Calls: calcarea.m; 
%        spacemarch.m (and hence march.m and shock.m)
%        steger_flux.m (or vanleer_flux.m)
%        flux.m
%        loadq.m
%        loadpr.m
%______________________________________________________________________
%
% Set up the initial data
%
clear;
for ifig=1:4
  figure(ifig); clf;
  set(gcf,'Units','pixels','Position',[(10+ifig)*30+10 (10-ifig)*30+10 750 450]);
end
gamma = 1.4;    % ratio of specific heats
global iter;

% Physical parameters
amach0 = 1.265;  % Mach number at the entrance 
rho0 = 0.50;    % density at the entrance 
p0 = 0.379;      % pressure at the entrance 
%
xlength = 10.;  % length of the nozzle
xsh = 4;      % location if use fixed shock location
%
irefine = 4;          % change this parameter to refine the space mesh (and time)
jmax = 40*irefine+1;  % number of mesh points along the nozzle
% needed for time marching
cfl = 0.90;           % Determines time step size, make less than 1
itsteady=1500*irefine;% number of iterations to perform to get steady initial solution
imeth = 1;           % method (1=Steger-Warming)
%
dpvar = 0.10;        % fraction of pressure variation at nozzle exit
timeper= 0.5;         % how fast is the pressure variation 
iper = 0;             % how many periods of pressure oscillation to simulate
% Calculate how many iterations per period of oscillation
% Note: this uses more iterations per period if mesh is refined in space
%       this will keep the cfl number roughly similar!
itperiod = irefine*400*timeper;
itplots = 16; % how many plots in one period
itmax = itsteady+itperiod*iper;



% [rho,u,p,e,amach] = loadpr(q,area,gamma);

err_list = [];


disp([' Getting ready to run ',num2str(itmax),...
      ' iterations, plotting every ',num2str(itperiod/10),]);
% pre-allocate size of some arrays for speed
dq=zeros(3,jmax); normit=zeros(1,itmax); pplot1=zeros(1,itmax);
%
% Initialize the grid
dx = xlength/(jmax-1); % nozzle of length xlength
x = 0:dx:xlength;
%
% calculate the area at each grid point and grid interfaces
area = calcarea(x);
areaint = calcarea(x+0.5*dx);
%
% Let's pick a location to plot (just downstream of shock)
xplot1 = xsh+0.05*xlength;
iloc1 = floor(xplot1/xlength*(jmax-1)+1);
f1=(xplot1-x(iloc1))/dx;
disp([' Plot pressure time history at location: ',num2str(xplot1)]);
% Initialize the solution by space marching (with corrector)
[rho_sp,u_sp,p_sp,e_sp,amach_sp]=...
            spacemarch(gamma,amach0,p0,rho0,xsh,x,area,1);
rho = rho_sp; u = u_sp; p = p_sp; e = e_sp; amach = amach_sp;
% Note the exit pressure from space marching
pend = p(jmax);
disp([' Calculated initial exit pressure is: ',num2str(pend)]);
% Load primitive variables into conservative variables
[q] = loadq(rho,u,e,area);


% compute exact solution
[rho_ex, u_ex, p_ex, e_ex] = exact_solution(x, gamma, amach0, p0, rho0, 1, xsh);




% Now let's iterate
for iter=1:itmax
% calculate dt/dx for steady flow and then freeze
  if( iter < itsteady )
    dtdx = min(cfl./(abs(u)+u./amach));
  end
  if(imeth==1) % 
% Load up the fluxes according to Steger-Warming
      [fluxp,fluxn] = steger_flux(gamma,area,rho,u,p,e);
% Flux at j+1/2
      fluxjp(:,1:jmax-1)=fluxp(:,1:jmax-1)+fluxn(:,2:jmax);
  end
%
% What is the change? (include the source term in the momentum)
  if(imeth<=1)
    dq(:,2:jmax-1)=-dtdx*(fluxjp(:,2:jmax-1)-fluxjp(:,1:jmax-2));
    dq(2,2:jmax-1)=dq(2,2:jmax-1)+dtdx*p(2:jmax-1).*(areaint(2:jmax-1)-areaint(1:jmax-2));
  end
% Boundary condition at j = 1 (assumed supersonic entrance)
% Boundary condition at j = jmax (subsonic or supersonic exit)
  [dq]=calcbc(gamma,dq,p,rho,u,area,dtdx,dpvar,pend,itperiod,iter,itsteady);
% Calculate L2norm of the residual
  normres(iter) = norm(sqrt(dq(1,:).^2+dq(2,:).^2+dq(3,:).^3),2);
  if(rem(iter,100) == 0)
    disp([' L2norm of residual is ',num2str(normres(iter)),...
          ' at iteration ',num2str(iter)]);
  end
% Calculate L2norm of the change
  normit(iter) = norm(sqrt(dq(1,:).^2+dq(2,:).^2+dq(3,:).^3),2);
  if(rem(iter,100) == 0)
%      pause(1)
    disp([' L2norm of change is ',num2str(normit(iter)),...
          ' at iteration ',num2str(iter)]);
  end



% Update conservative variables
  q = q + dq;

  
  
  
  % Load conservative variables into primitive variables
  [rho,u,p,e,amach] = loadpr(q,area,gamma);


  % compute error of numerical against exact solution
      error = norm(rho - rho_ex, 'fro')/sqrt(jmax);
  err_list = [err_list, error];
  % disp([' L2norm of error is ',num2str(error),...
          % ' at iteration ',num2str(iter)]);



  % What is the pressure at xplot1?
  pplot1(iter) = (1-f1)*p(iloc1)+f1*p(iloc1+1);
% let's store steady pressure distribution and exit pressure
  if(iter==itsteady) 
    pend = p(jmax); 
	p_steady = p;
    amach_steady=amach;
  end
% let's plot pressure and Mach every so many iterations
  if(rem(iter-itsteady,itperiod/itplots)==0)
    figure(1)
    plot(x,p_sp,'m-','LineWidth',2.0)
    hold on;
    plot(x,p,'r-','LineWidth',2.0)
    plot(x,p_ex,'g-','LineWidth',2.0)
    if(iter > itsteady) plot(x,p_steady,'g-','LineWidth',2.0); end
    hold off;
% This is the last oscillatory period so overlay pressure plots
    if(iter > itmax-itperiod) 
      hold on; 
    end
    axis([0.,xlength,0.,1.1]);
    title(' Pressure Distribution Along the Nozzle')
    xlabel('Distance along the Nozzle'); ylabel('Pressure')
    set(gca,'FontSize',16,'LineWidth',2.0,'FontWeight','demi');
    figure(4)
    plot(x,amach_sp,'m-','LineWidth',2.0)
    hold on;
    plot(x,amach,'r-','LineWidth',2.0)
    axis([0.,xlength,0.,2.5]);
    if(iter > itsteady) plot(x,amach_steady,'g-','LineWidth',2.0); end
    hold off;
% This is the last oscillatory period so overlay mach plots
    if(iter > itmax-itperiod) 
      hold on; 
    end
    axis([0.,xlength,0.,2.5]);
    title(' Mach Distribution Along the Nozzle')
    xlabel('Distance along the Nozzle'); ylabel('Mach')
    set(gca,'FontSize',16,'LineWidth',2.0,'FontWeight','demi');
  end
end
% finished the simulation, let's plot final pressure distribution on top!
figure(1)
plot(x,p,'y-','LineWidth',2.0)
hold off;
% let's plot pressure time history at location as a function of time
figure(2)
time = (1:iter-itsteady)/itperiod; % Fraction of oscillation period
plot(time,pplot1(itsteady+1:iter),'ro','LineWidth',2.0)
set(gca,'FontSize',16,'LineWidth',2.0,'FontWeight','demi');
title(' Time History at Point in Nozzle Just Downstream of Shock')
xlabel('Time, fraction of period'); ylabel('Pressure')
% let's plot L2 norm of change as a function of iteration
figure(3)
semilogy(1:iter,normit,'ro','LineWidth',2.0)
hold on;
semilogy(1:iter,normres,'r--','LineWidth',2.0)
semilogy(1:iter,err_list,'g-','LineWidth',2.0)
set(gca,'FontSize',16,'LineWidth',2.0,'FontWeight','demi');
title(' Convergence Time History')
xlabel('Iteration'); ylabel('Residual')
% finished!
