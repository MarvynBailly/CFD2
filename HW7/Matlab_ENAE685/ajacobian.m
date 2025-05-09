% ajacobian.m
%
function [lmat,umat,bscalinv] = ajacobian(dtdx,gamma,rho,u,p,e)
    % THIS IS ABAR I THINK 
    m=length(rho);
    lmat=zeros(3,3,m);
    umat=zeros(3,3,m);

    gamm1 = gamma-1;
    gamm3 = gamma-3;
	speed = sqrt(gamma*p./rho);

	sprad = 1.05*(abs(u)+speed);
	bscalinv = 1./(1+dtdx*sprad);

    ajac(1,1,:) = zeros(1,m);
	ajac(1,2,:) = ones(1,m);
	ajac(1,3,:) = zeros(1,m);
	ajac(2,1,:) = 0.5*gamm3.*u.^2;
	ajac(2,2,:) = -gamm3.*u;
	ajac(2,3,:) = gamm1;
	ajac(3,1,:) = gamm1*u.^3-gamma.*e.*u./rho;
	ajac(3,2,:) = gamma*e./rho-1.5*gamm1*u.^2;
	ajac(3,3,:) = gamma*u;

    lmat(:,:,2:m-1) = -0.5*dtdx.*ajac(:,:,1:m-2);
    umat(:,:,2:m-1) = 0.5*dtdx.*ajac(:,:,3:m);
    for jj=2:m-1
      lmat(1,1,jj)=lmat(1,1,jj)-0.5*dtdx*sprad(jj-1);
      lmat(2,2,jj)=lmat(2,2,jj)-0.5*dtdx*sprad(jj-1);
      lmat(3,3,jj)=lmat(3,3,jj)-0.5*dtdx*sprad(jj-1);
      umat(1,1,jj)=umat(1,1,jj)-0.5*dtdx*sprad(jj+1);
      umat(2,2,jj)=umat(2,2,jj)-0.5*dtdx*sprad(jj+1);
      umat(3,3,jj)=umat(3,3,jj)-0.5*dtdx*sprad(jj+1);
    end
end
