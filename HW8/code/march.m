% march.m
%
% Space march from one area to the next using predictor-corrector
%
function [rho_n,u_n,p_n,amach_n] = ...
         march(area_n,gamma,rho_o,u_o,p_o,amach_o,area_o,i_correct)

    area_h = 0.5*(area_o+area_n);
%
% predictor
%
    rho_n   = rho_o + (area_n-area_o)/area_h*rho_o*amach_o^2/(1-amach_o^2);
    u_n     = u_o*rho_o*area_o/(rho_n*area_n);
    p_n     = rho_n*( 0.5*(gamma-1.)/gamma*(u_o^2-u_n^2)+p_o/rho_o );
    amach_n = u_n/sqrt(gamma*p_n/rho_n);
%
% corrector (twice)
%
    if (i_correct == 1)
	  for it=1:2
	    rho_h  = 0.5*(rho_o+rho_n);
	    amach_h= 0.5*(amach_o+amach_n);
        rho_n   = rho_o + (area_n-area_o)/area_h*rho_h*amach_h^2/(1-amach_h^2);
        u_n     = u_o*rho_o*area_o/(rho_n*area_n);
        p_n     = rho_n*( 0.5*(gamma-1.)/gamma*(u_o^2-u_n^2)+p_o/rho_o );
        amach_n = u_n/sqrt(gamma*p_n/rho_n);
	  end
    end
end
