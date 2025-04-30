% calcbc.m
%
% Calculate entrance and exit bc
%
function [dq]=calcbc(gamma,dq,p,rho,u,area,dtdx,dpvar,pend,itperiod,iter,itsteady)
% Determine length of the array
  jmax = length(area);
% Assumed supersonic at entrance
  j=1;
  dq(1,j) = 0;
  dq(2,j) = 0;
  dq(3,j) = 0;
%   dq(1,j) = 0.95*dq(1,j+1);
%   dq(2,j) = 0.95*dq(2,j+1);
%   dq(3,j) = 0.95*dq(3,j+1);
% Allow for subsonic or supersonic at exit
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
% Test here and set dp appropriately
  if (u(j-1) >= speedm)   % supersonic exit
    dp = 0.5*(r2+r3);
  else
    if(iter > itsteady )  % subsonic exit
      timenew = (iter-itsteady)/itperiod; % in terms of oscillation period
      timeold = (iter-itsteady-1)/itperiod; % in terms of oscillation period
	  dp = dpvar*pend*(sin(2*pi*timenew)-sin(2*pi*timeold));
	else
	  dp = 0;
	end
  end
  drho = r1+dp*sp2inv;
  du = (r2-dp)/(rhospeed);
  dq(1,j) = area(j)*drho;
  dq(2,j) = area(j)*( (rho(j)+drho)*(u(j)+du) - rho(j)*u(j) );
  dq(3,j) = area(j)*( ( (p(j)+dp)/(gamma-1)+0.5*(rho(j)+drho)*(u(j)+du)^2 ) ...
                     -( p(j)/(gamma-1)+0.5*rho(j)*u(j)^2 ) );
end
