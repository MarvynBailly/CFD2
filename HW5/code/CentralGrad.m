function grad = CentralGrad(f, direction, grid_size)
% compute_gradient_vectorized computes the gradient of a 2D array in the specified direction.
%
%   grad = compute_gradient_vectorized(f, direction, grid_size) returns the gradient of the
%   2D array f. The input 'direction' specifies which derivative to compute:
%       - 'horizontal' for differences along columns (x-direction)
%       - 'vertical' for differences along rows (y-direction)
%   grid_size is the spacing between grid points.
%
%   The function uses central differences for interior points and one-sided
%   differences at the boundaries.

    [nrows, ncols] = size(f);
    grad = zeros(nrows, ncols);

    switch lower(direction)
        case 'horizontal'
            % Left boundary (forward difference)
            grad(:, 1) = (f(:, 2) - f(:, 1)) / grid_size;
            % Interior points (central difference)
            grad(:, 2:end-1) = (f(:, 3:end) - f(:, 1:end-2)) / (2 * grid_size);
            % Right boundary (backward difference)
            grad(:, end) = (f(:, end) - f(:, end-1)) / grid_size;

        case 'vertical'
            % Top boundary (forward difference)
            grad(1, :) = (f(2, :) - f(1, :)) / grid_size;
            % Interior points (central difference)
            grad(2:end-1, :) = (f(3:end, :) - f(1:end-2, :)) / (2 * grid_size);
            % Bottom boundary (backward difference)
            grad(end, :) = (f(end, :) - f(end-1, :)) / grid_size;

        otherwise
            error('Invalid direction. Use ''horizontal'' or ''vertical''.');
    end
end
