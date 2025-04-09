% findshock.m
%______________________________________________________________________
function [xsh]= findshock(gamma,fsmach,p0,rho0,xsh1,xsh2,x,area,icorr,ptarget);
%
% Let's find shock location to give target exit pressure
% use around 8 iterations of Newton-Raphson?
%
xshold = xsh1;
[rho,u,p,e,amach]=...
          spacemarch(gamma,fsmach,p0,rho0,xshold,x,area,icorr);
pold = p(length(p));
xshnew = xsh2;
for it=1:8
  [rho,u,p,e,amach]=...
            spacemarch(gamma,fsmach,p0,rho0,xshnew,x,area,icorr);
  pnew = p(length(p));
  if(pnew-pold == 0)
    dxsh = 0.;
  else
    dxsh=-(pnew-ptarget)*(xshnew-xshold)/(pnew-pold);
  end
  xshold=xshnew;
  xshnew=xshold+dxsh;
% for debugging the iterations!!
%  disp([num2str(it),' ',num2str(xshold),' ',num2str(xshnew),...
%                    ' ',num2str(pold),' ',num2str(pnew)]);
  pold=pnew;
end
if(abs(pold-ptarget) > 0.0001*ptarget)
  disp([' ** MISSED TARGET PRESSURE BY ** ',num2str(abs(pold-ptarget)/ptarget*100),' %']);
end
xsh = xshnew;
%end
