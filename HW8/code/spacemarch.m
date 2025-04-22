% spacemarch.m
%
% Space march along the nozzle w/ or w/o shock at xsh
%
function [rho,u,p,e,amach]=spacemarch(gamma,fsmach,p0,rho0,xsh,x,area,i_correct)
% Determine length of the array to space march along
  jmax = length(area);
%
% Initialize upstream condition
%
  rho(1) = rho0;
  p(1) = p0;
  amach(1) = fsmach;
  u(1) = amach(1)*sqrt(gamma*p(1)/rho(1));
%
% Let's perform the space marching and note if also need to apply shock jump
%
  for j=2:jmax
%
% Space marching from j-1 to j
%
%   Is the shock within the cell?
%
    if ( x(j-1) <= xsh & x(j) > xsh )  
      %     Shock is in the cell!
      %     Let's calculate the area at the shock
      areash = calcarea(xsh);
      %     Space march up to the shock, unless shock is exactly at x(j-1)
      if( xsh-x(j-1) > 0 )
          [rho_o,u_o,p_o,amach_o] = ...
        march(areash,gamma,rho(j-1),u(j-1),p(j-1),amach(j-1),area(j-1),i_correct);
      else
        rho_o=rho(j-1); u_o=u(j-1); p_o=p(j-1); amach_o=amach(j-1);
      end
      %   Apply the correct shock jump condition
      [rho_n,u_n,p_n,amach_n] = shock(gamma,rho_o,u_o,p_o,amach_o);
      %   Now spacepace march after the shock
      [rho(j),u(j),p(j),amach(j)] = ...
	    march(area(j),gamma,rho_n,u_n,p_n,amach_n,areash,i_correct);
    else
      %   Shock is not in the cell - just space march!
      [rho(j),u(j),p(j),amach(j)] = ...
	    march(area(j),gamma,rho(j-1),u(j-1),p(j-1),amach(j-1),area(j-1),i_correct);
    end
  end
  % Update the energy along the whole nozzle
  e = 0.5*rho.*u.^2+p./(gamma-1);
end
