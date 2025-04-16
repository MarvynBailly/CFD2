


function Qh = convert_to_conservative(Qh, rho, u, e, area)
    % convert to conservative variables
    Qh(1, :) = rho .* area;
    Qh(2, :) = rho .* u .* area;
    Qh(3, :) = e .* area;
  end
  
  function Q = convert_to_primitive(Q, Qh, gamma, area)
    % convert to primitive variables
    Q(1, :) = Qh(1, :) ./ area;
    Q(2, :) = Qh(2, :) ./ Qh(2, :);
    Q(3, :) = (gamma - 1) .* (Qh(3,:) ./ area - 0.5 .* Q(1,:) .* Q(2,:).^2);
  end
  
  function dt = compute_timestep(cfl, dx, u, c)
    % compute local timesteps
    dt_local = cfl * dx ./ (abs(u) + abs(c));
    dt = min(dt_local, [], "all");
  end