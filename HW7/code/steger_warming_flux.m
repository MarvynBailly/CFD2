function [Fhp, Fhm] = steger_warming_flux(Qh, A, gamma)
% Steger-Warming flux vector splitting
% assumes 1D Euler equations in conservative form
%   Qh = [rho*A, rho*u*A, rho*e*A] 
% so have to divide by A and play around to get back the primitive variables
%   Q = [rho, u, p]
% return fluxes in conservative form by multiplying by A

    N = size(Q, 2);
    Fp = zeros(3, N);
    Fm = zeros(3, N);

    for j = 1:N
        % convert back to primitive variables
        rho = Q(1, j) / A(j);
        rhou = Q(2, j) / A(j);        
        e = Q(3, j) / A(j);
        
        u = rhou / rho;
        u2 = u^2;
        p = (gamma - 1) * (e - 0.5 * rho * u2);

        % compute speed of sound and enthalpy
        c = sqrt(gamma * p / rho);
        H = (e + p) / rho;

        % compute eigenvalues
        lambda1 = u;
        lambda2 = u + c;
        lambda3 = u - c;

        % split the evals
        lambda1p = 0.5 * (lambda1 + abs(lambda1)); lambda1m = 0.5 * (lambda1 - abs(lambda1));
        lambda2p = 0.5 * (lambda2 + abs(lambda2)); lambda2m = 0.5 * (lambda2 - abs(lambda2));
        lambda3p = 0.5 * (lambda3 + abs(lambda3)); lambda3m = 0.5 * (lambda3 - abs(lambda3));
        
        % compute fluxes directly from provided formula
        coef = 1 / (2 * gamma);
        Fp(1, j) = coef * (2*(gamma - 1) * lambda1p + lambda2p + lambda3p);
        Fp(2, j) = coef * (2*(gamma - 1) *u * lambda1p + lambda2 * lambda2p + lambda3 * lambda3p);
        Fp(3, j) = coef * (2*(gamma - 1) * u2 * lambda1p + (H + u*c) * lambda2p + (H - u*c) * lambda3p);
        
        Fm(1, j) = coef * (2*(gamma - 1) * lambda1p + lambda2p + lambda3p);
        Fm(2, j) = coef * (2*(gamma - 1) *u * lambda1p + lambda2 * lambda2p + lambda3 * lambda3p);
        Fm(3, j) = coef * (2*(gamma - 1) * u2 * lambda1p + (H + u*c) * lambda2p + (H - u*c) * lambda3p);
    end

    % convert back to conservative form
    Fhp = Fp .* A;
    Fhm = Fm .* A;
end