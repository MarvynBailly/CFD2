% shock.m
%
% Shock Jump Relations for Euler Equations
%
function [rho_n,u_n,p_n,amach_n] = shock(gamma,rho_o,u_o,p_o,amach_o)
%
% Test to make sure that the Mach number ahead of shock is supersonic!
%
  if(amach_o < 1.0) 
    disp([' ** WARNING  ** Mach number ahead of shock is only '...
          ,num2str(amach_o)]);
    disp([' ** WARNING  ** it should be greater than 1!']);
  end
%
% Apply Rankine-Hugoniot Jump for Euler Equations
%
  rho_n = rho_o + rho_o*(amach_o^2-1)/(0.5*(gamma-1)*amach_o^2+1.);
% conservation of momentum
  u_n = u_o*rho_o/rho_n;
% conservation of energy
  p_n = rho_n * ( 0.5*(gamma-1.)/gamma*(u_o^2-u_n^2)+p_o/rho_o );
% update the mach number downstream by definition
  amach_n = u_n/sqrt(gamma*p_n/rho_n);
end
