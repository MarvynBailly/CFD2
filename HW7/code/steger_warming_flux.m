function [Fp, Fm] = steger_warming_flux(Q, gamma)
% Steger-Warming flux vector splitting
%   Q = [rho, rho*u, rho*E] = [rho, rho*u, rho*(e + 0.5*u^2)]
    N = size(Q, 2);
    Fp = zeros(3, N);
    Fm = zeros(3, N);


end