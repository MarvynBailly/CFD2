%% grid_generation_main.m
% This template sets up the grid generation problem for an O-grid around a
% NACA00xx airfoil. The code is modularized into functions for:
%   - Initializing grid and parameters
%   - Defining airfoil and outer boundary conditions
%   - Creating an algebraic grid (initial guess)
%   - Setting elliptic boundary conditions (with and without control)
%   - Calling the solver (solver section is intentionally left empty)
%   - Plotting the resulting grid

clear; clc; close all;

%% PARAMETERS
% Grid dimensions and locations (as specified)
JMAX  = 217;       % Number of points in wrap-around (airfoil contour) direction
KMAX  = 71;        % Number of points in the normal direction
JLE   = 109;       % Index corresponding to the leading edge
ds    = 0.004;     % Spacing to first point off the airfoil
XSF   = 1.02;      % Stretching factor for the normal direction
XMAX  = 1.5;       % Outer boundary x-location (assumed; adjust as needed)
tau   = 0.12;      % Thickness parameter for the NACA00xx airfoil
omega = 1.5;

% Choose which method to run:
% 1 -> Algebraic Grid Generation (initial guess)
% 2 -> Elliptic grid generation, no control (P = Q = 0)
% 3 -> Elliptic grid generation, with control (P and Q computed via a control method)
method = 2;  % Change this value to 2 or 3 for the other cases

gd.x = zeros(JMAX, KMAX);
gd.y = zeros(JMAX, KMAX);

%% APPLY BOUNDARY CONDITIONS
% Generate and apply the airfoil boundary and outer boundaries based on the problem statement.
% This routine fills in the boundary values for x and y using the given formulas.
params.JMAX = JMAX;
params.KMAX = KMAX;
params.JLE  = JLE;
params.ds   = ds;
params.XSF  = XSF;
params.XMAX = XMAX;
params.tau  = tau;
params.dxi  = 1 / (JMAX-1);
params.deta = 1 / (KMAX-1);
params.omega = omega;
params.max_iter = 2000;
params.tol = 1e-6;
gd = set_boundaries(gd, params);
[T1, T2] = TFI(flip([gd.x(:,1)';gd.y(:,1)'],2), [gd.x(1,:);gd.y(1,:)], ...
           [gd.x(:,end)';gd.y(:,end)'], [gd.x(end,:);gd.y(end,:)]);
gd.x = T1;
gd.y = T2;

%% (OPTIONAL) SOLVER: ELLIPTIC GRID GENERATION
% For methods 2 and 3, solve the elliptic grid generation equations.
if method >= 2
    % Initialize convergence criteria, relaxation factors, etc.
    % NOTE: The solver routine is left empty. Insert your iterative scheme here.
    gd = solve_elliptic(gd, params, method);
end

%% PLOT THE GRID
plot_grid(gd);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTION DEFINITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gd = set_boundaries(gd, params)
    
    JMAX = params.JMAX;
    KMAX = params.KMAX;
    JLE  = params.JLE;
    ds   = params.ds;
    XSF  = params.XSF;
    XMAX = params.XMAX;
    tau  = params.tau;
    
    %% Airfoil Boundary (j-index = 1)
    % Cosine stretching for the x-locations on the airfoil.
    % Here, xi varies from 0 (leading edge) to 1 (trailing edge).
    % Airfoil x locations: x(j,1)
    for j = 1:JMAX
        % Using cosine clustering (adjust formula if necessary)
        theta = pi * (JLE - j) / (JLE - 1);  % For j = 1,...,JLE
        % For j beyond the leading edge, symmetry is assumed.
        % Note: Adjust the clustering for the lower and upper surfaces as needed.
        gd.x(j,1) = 0.5 - 0.5*cos(theta);
        %if j <= JLE
        %    gd.x(j,1) = 0.5 - 0.5*cos(theta);
        %else
        %    gd.x(j,1) = 0.5 - 0.5*cos(pi - theta);
        %end
        
        % Compute the corresponding y values using the analytical expression
        % for a NACA00xx airfoil.
        % f(x) = 5*tau*(0.2969*sqrt(x) - 0.1260*x - 0.3516*x^2 + 0.2843*x^3 - 0.1015*x^4)
        f = 5*tau*(0.2969*sqrt(gd.x(j,1)) - 0.1260*gd.x(j,1) ...
                  - 0.3516*gd.x(j,1)^2 + 0.2843*gd.x(j,1)^3 - 0.1015*gd.x(j,1)^4);
        % For the upper surface use +f, for the lower surface use -f.
        if j <= JLE
            gd.y(j,1) = f;
        else
            gd.y(j,1) = -f;
        end
    end
    % Force trailing-edge closure: set y at trailing edge (x=1.0) to zero.
    % Assume trailing edge corresponds to j = JLE and j = JLE+1.
    %gd.y(JLE,1)   = 0.0;
    gd.y(1,1)   = 0.0;
    %gd.y(JLE+1,1) = 0.0;
    gd.y(end,1)   = 0.0;

    gd.y(2,1) = 0.5*gd.y(2,1) + 0.25*(gd.y(3,1)+gd.y(1,1));
    
    % Optionally smooth the points near the trailing-edge (for j = JMAX-? region)
    % Example smoothing for adjacent points:
    if JMAX > (JLE+2)
        gd.y(JMAX-1,1) = 0.5*gd.y(JMAX-1,1) + 0.25*gd.y(JMAX,1) + 0.25*gd.y(JMAX-2,1);
    end
    
    %% Normal (K-direction) Boundaries: Rearward from the airfoil
    % j = 1: boundary along the airfoil surface is already set.
    % For k = 2:KMAX, set the grid by geometric stretching starting from the airfoil.
    % The first normal offset is given by ds.
    % For the rear boundary, a stretching factor XSF is applied.
    for k = 2:KMAX
        % For the first grid point along the normal from the airfoil
        if k == 2
            % At the airfoil surface, use ds to define the first offset.
            gd.x(1,k) = gd.x(1,1) + ds;  % For example, adjust along x-direction.
            gd.y(1,k) = 0;       % y remains the same (can be modified if needed)
        else
            % For k >= 3, apply the stretching factor XSF
            % Here we assume a one-dimensional stretching along the normal direction.
            gd.x(1,k) = gd.x(1,k-1) + (gd.x(1,k-1) - gd.x(1,k-2)) * XSF;
            gd.y(1,k) = 0.0;  % Normal line remains horizontal (adjust if necessary)
        end
        
        % Similarly, for the far-boundary (j = JMAX) use symmetry or define by
        % the outer circular boundary.
        gd.x(JMAX,k) = gd.x(1,k);
        gd.y(JMAX,k) = gd.y(1,k);
    end
    XMAX = gd.x(1,end); 
    
    %% Outer Boundary (k = KMAX) along a circle
    % Points along the outer boundary are placed along a circle centered at the
    % origin, starting at (XMAX, 0) and going clockwise.
    for j = 1:JMAX
        % Angle definition for equal spacing along the outer boundary.
        % This simple example assumes uniform angular spacing.
        %theta = (j-1)/(JMAX-1)*pi; % Change range as needed for full circle
        theta = (JLE-j)/(JLE-1)*pi;
        gd.x(j,KMAX) = -XMAX * cos(theta);
        gd.y(j,KMAX) = -XMAX * sin(theta);
    end
end

function gd = solve_elliptic(gd, params, method)
    res_list = [];
    Dx = zeros(size(gd.x));
    Dy = zeros(size(gd.y));
    for iter=1:params.max_iter
        dxdxi = CentralGrad(gd.x, 'horizontal', params.dxi);
        dxdeta = CentralGrad(gd.x, 'vertical', params.deta);
        dydxi = CentralGrad(gd.y, 'horizontal', params.dxi);
        dydeta = CentralGrad(gd.y, 'vertical', params.deta);
    
        A1 = dxdeta.^2 + dydeta.^2;
        A2 = dxdeta.*dxdxi + dydeta.*dydxi;
        A3 = dxdxi.^2 + dydxi.^2;
    
        % x-component
        res_x = A1(2:end-1,2:end-1).*CentralD2(gd.x, 'horizontal', params) ...
                -2*A2(2:end-1,2:end-1).*CentralD2(gd.x, 'mixed', params) ...
                +A3(2:end-1,2:end-1).*CentralD2(gd.x, 'vertical', params);
    
        for k=2:params.KMAX-1
            f_x = -res_x(k-1,:) + 2*A2(k,2:end-1).*(Dx(k-1,1:end-2)-Dx(k-1,3:end)) / (4*params.dxi*params.deta) ...
                  -A3(k,2:end-1).*Dx(k-1,2:end-1) / params.deta^2;
    
            % Update Dirichlet's BC
            %f_x(1) = f_x(1) - A1(k,2) / params.dxi.^2;
            %f_x(end) = f_x(end) - A1(k,end-1) / params.dxi.^2;
    
            a = A1(k,3:end)/params.dxi.^2;
            b = -2*A1(k,2:end-1)/params.dxi.^2 - 2*A3(k,2:end-1)/params.deta.^2;
            c = A1(k,1:end-2)/params.dxi.^2;
    
            A = spdiags([a', b', c'], [-1, 0, 1], params.JMAX-2, params.JMAX-2);
            Dx(k,2:end-1) = (A\f_x')';
        end
    
        % y-component
        res_y = A1(2:end-1,2:end-1).*CentralD2(gd.y, 'horizontal', params) ...
                -2*A2(2:end-1,2:end-1).*CentralD2(gd.y, 'mixed', params) ...
                +A3(2:end-1,2:end-1).*CentralD2(gd.y, 'vertical', params);
        for k=2:params.KMAX-1
            f_y = -res_y(k-1,:) + 2*A2(k,2:end-1).*(Dy(k-1,1:end-2)-Dy(k-1,3:end)) / (4*params.dxi*params.deta) ...
                  -A3(k,2:end-1).*Dy(k-1,2:end-1) / params.deta^2;
            % Update Dirichlet's BC
            %f_y(1) = f_y(1) - A1(k,2) / params.dxi.^2;
            %f_y(end) = f_y(end) - A1(k,end-1) / params.dxi.^2;
    
            a = A1(k,3:end)/params.dxi.^2;
            b = -2*A1(k,2:end-1)/params.dxi.^2 - 2*A3(k,2:end-1)/params.deta.^2;
            c = A1(k,1:end-2)/params.dxi.^2;
    
            A = spdiags([a', b', c'], [-1, 0, 1], params.JMAX-2, params.JMAX-2);
            Dy(k,2:end-1) = A\f_y';
        end
    
        gd.x = gd.x + params.omega*Dx;
        gd.y = gd.y + params.omega*Dy;
        res = norm([res_x, res_y], 'fro');
        res_list = [res_list, res];

        
    end
    plot(gd.x, gd.y, 'b.-'); hold on;
    plot(gd.x', gd.y', 'b.-');
    title('Full Grid');
    xlabel('x'); ylabel('y');
    axis equal; grid on;

    figure();
    semilogy(res_list); grid on
end


function plot_grid(gd)
    % Plot the grid in three different views:
    %   1. Full grid including outer boundaries.
    %   2. Zoomed view around the airfoil.
    %   3. Zoomed view around the trailing edge.
    %
    % Input:
    %   grid - structure with grid.x and grid.y coordinates.
    
    [JMAX, KMAX] = size(gd.x);
    
    % Full gd plot
    figure;
    plot(gd.x, gd.y, 'b.-'); hold on;
    plot(gd.x', gd.y', 'b.-');
    title('Full Grid');
    xlabel('x'); ylabel('y');
    axis equal; grid on;
    
    % Zoom around the airfoil (near j-index corresponding to the airfoil boundary)
    figure;
    % Choose indices near the airfoil surface
    airfoil_j = round(JMAX/2);
    plot(gd.x(airfoil_j,:), gd.y(airfoil_j,:), 'r.-');
    title('Zoomed View Around the Airfoil');
    xlabel('x'); ylabel('y');
    axis equal; grid on;
    
    % Zoom around the trailing-edge (assume trailing edge is near j = JLE)
    figure;
    plot(gd.x(params.JLE-5:params.JLE+5, :), gd.y(params.JLE-5:params.JLE+5, :), 'k.-');
    title('Zoomed View Around the Trailing Edge');
    xlabel('x'); ylabel('y');
    axis equal; grid on;
end

function [T1, T2] = TFI(bottom, right, top, left)
    Nx1 = length(bottom);
    Nx2 = length(right);

    rb = right(:,1);
    rt = right(:,end);
    lt = left(:,end);
    lb = left(:,1);
    
    s1 = linspace(0,1,Nx1);
    s2 = linspace(0,1,Nx2);
    [S1, S2] = meshgrid(s1, s2);
    
    T1 = (1-S2).*bottom(1,:) + S2.*top(1,:) + (1-S1).*right(1,:)' + S1.*left(1,:)' ...
    - (1-S1).*(1-S2)*rb(1) - S1.*S2*lt(1) - S1.*(1-S2)*lb(1) - (1-S1).*S2*rt(1);
    
    T2 = (1-S2).*bottom(2,:) + S2.*top(2,:) + (1-S1).*right(2,:)' + S1.*left(2,:)' ...
    - (1-S1).*(1-S2)*rb(2) - S1.*S2*lt(2) - S1.*(1-S2)*lb(2) - (1-S1).*S2*rt(2);
end
