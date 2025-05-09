% steger_flux.m
function [fluxp,fluxn] = steger_flux(gamma,area,rho,u,p,e);
	% disp(size(rho));
	% disp(size(u));
	% disp(size(p));

    eps2 = 0.000001;
    gamm1 = gamma-1;
	speed = sqrt(gamma*p./rho);
	alam1 = u;
	alam2 = u+speed;
	alam3 = u-speed;
	alam1 = 0.5*(alam1+sqrt(alam1.^2+eps2));
	alam2 = 0.5*(alam2+sqrt(alam2.^2+eps2));
	alam3 = 0.5*(alam3+sqrt(alam3.^2+eps2));
	fac = rho*0.5/gamma.*area;
	fluxp(1,:) = fac.*( 2.*gamm1*alam1+alam2+alam3 );
	fluxp(2,:) = fac.*( 2.*gamm1*u.*alam1+alam2.*(u+speed)+alam3.*(u-speed) );
	fluxp(3,:) = fac.*( gamm1*u.^2.*alam1+0.5*alam2.*(u+speed).^2+0.5*alam3.*(u-speed).^2 + ...
	                0.5*(3-gamma)/gamm1*(alam2+alam3).*speed.^2 );
	fluxn(1,:) = rho.*u.*area-fluxp(1,:);
	fluxn(2,:) = (rho.*u.^2+p).*area-fluxp(2,:);
	fluxn(3,:) = (e+p).*u.*area-fluxp(3,:);
end
