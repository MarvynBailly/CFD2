function [Fhp, Fhm] = steger_warming_flux(Q, A, gamma)
% Steger-Warming flux vector splitting
% return fluxes in conservative form by multiplying by A

    N = size(Q, 2);
    Fp = zeros(3, N);
    Fm = zeros(3, N);

    rho = Q(1, :);
    u = Q(2, :);
    p = Q(3, :);

    e = (p / (gamma - 1)) + 0.5 * rho .* u.^2;
    u2 = u.^2;
    c = sqrt(gamma * p ./ rho);
    H = (e + p) ./ rho;

    lambda1 = u;
    lambda2 = u + c;
    lambda3 = u - c;

    lambda1p = 0.5 * (lambda1 + abs(lambda1)); lambda1m = 0.5 * (lambda1 - abs(lambda1));
    lambda2p = 0.5 * (lambda2 + abs(lambda2)); lambda2m = 0.5 * (lambda2 - abs(lambda2));
    lambda3p = 0.5 * (lambda3 + abs(lambda3)); lambda3m = 0.5 * (lambda3 - abs(lambda3));

    % compute fluxes directly from provided formula
    coef = 1 / (2 * gamma);
    Fp(1, :) = coef * (2*(gamma - 1) * lambda1p + lambda2p + lambda3p);
    Fp(2, :) = coef * (2*(gamma - 1) *u .* lambda1p + lambda2 .* lambda2p + lambda3 .* lambda3p);
    Fp(3, :) = coef * (2*(gamma - 1) * u2 .* lambda1p + (H + u.*c) .* lambda2p + (H - u.*c) .* lambda3p);
    
    Fm(1, :) = coef * (2*(gamma - 1) * lambda1p + lambda2p + lambda3p);
    Fm(2, :) = coef * (2*(gamma - 1) * u .* lambda1p + lambda2 .* lambda2p + lambda3 .* lambda3p);
    Fm(3, :) = coef * (2*(gamma - 1) * u2 .* lambda1p + (H + u.*c) .* lambda2p + (H - u.*c) .* lambda3p);

    % convert back to conservative form
    Fhp = Fp .* A;
    Fhm = Fm .* A;
end