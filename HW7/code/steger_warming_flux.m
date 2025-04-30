function [Fhp, Fhm] = steger_warming_flux(Q, A, gamma)
    % it looks like this is all good
    N = size(Q, 2);
    Fhp = zeros(3, N);
    Fhm = zeros(3, N);

    rho = Q(1, :);
    u   = Q(2, :);
    p   = Q(3, :);

    a = sqrt(gamma * p ./ rho);
    H = (p / (gamma - 1) + 0.5 * rho .* u.^2 + p) ./ rho;  % total enthalpy

    eps2 = 0.000001;

    % Split eigenvalues
    lambda1p = 0.5 * (u     + sqrt(u.^2 + eps2));
    lambda2p = 0.5 * (u + a + sqrt((u + a).^2 + eps2));
    lambda3p = 0.5 * (u - a + sqrt((u - a).^2 + eps2));

    lambda1m = 0.5 * (u     - sqrt(u.^2 + eps2));
    lambda2m = 0.5 * (u + a - sqrt((u + a).^2 + eps2));
    lambda3m = 0.5 * (u - a - sqrt((u - a).^2 + eps2));

    for j = 1:N
        % Positive flux
        Fhp(:,j) = 0.5 * (1/gamma) * A(j) * rho(j) * [ ...
            2 * (gamma - 1) * lambda1p(j) + lambda2p(j) + lambda3p(j); ...
            2 * (gamma - 1) * u(j) * lambda1p(j) + (u(j) + a(j)) * lambda2p(j) + (u(j) - a(j)) * lambda3p(j); ...
            (gamma - 1) * u(j)^2 * lambda1p(j) + ...
            (H(j) + u(j)*a(j)) * lambda2p(j) + (H(j) - u(j)*a(j)) * lambda3p(j)];

        % Negative flux
        Fhm(:,j) = 0.5 * (1/gamma) * A(j) *  rho(j) * [ ...
            2 * (gamma - 1) * lambda1m(j) + lambda2m(j) + lambda3m(j); ...
            2 * (gamma - 1) * u(j) * lambda1m(j) + (u(j) + a(j)) * lambda2m(j) + (u(j) - a(j)) * lambda3m(j); ...
            (gamma - 1) * u(j)^2 * lambda1m(j) + ...
            (H(j) + u(j)*a(j)) * lambda2m(j) + (H(j) - u(j)*a(j)) * lambda3m(j)];
    end

end
