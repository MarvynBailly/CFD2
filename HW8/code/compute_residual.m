function residual = compute_residual(Qh, Q, Fhp, Fhm, area, dx, dt, x)
    N = size(Qh, 2);
    residual = zeros(3, N);
    for j = 2:size(Qh, 2)-1
        residual(:, j) = -dt / dx * (Fhp(:, j) - Fhp(:, j-1) + Fhm(:, j+1) - Fhm(:, j));
        residual(2, j) = residual(2,j)  + dt / dx * (Q(3,j) * abs((0.5*(area(j+1)+area(j)) - 0.5*(area(j)+area(j-1)))));
        %residual(2, j) = residual(2,j)  + dt / dx * (Q(3,j) * (calcarea(j*dx + 0.5 *dx) - calcarea(j*dx - 0.5 *dx)));
    end
end


