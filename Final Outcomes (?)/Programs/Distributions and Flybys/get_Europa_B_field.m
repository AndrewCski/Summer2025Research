function B = get_Europa_B_field(point, period_time)

    % Function for retrieving the magnetic field vector at a specific point
    % and specific time in the Europa system. "point" is a vector
    % describing the position in meters to be calculated for, "period_time"
    % is a time in seconds describing the point in Europa's synodic period 
    % to be calculated for.

    point = point(:);
    B0_vec = [75e-9; 200e-9; 0];
    omega = 2*pi/(11.23*3600); % synodic frequency (rad/s)
    B_prim = [B0_vec(1) * cos(omega * period_time); B0_vec(2) * -sin(omega * period_time); 0];
    mu0 = 4*pi*1e-7;           % vacuum permeability (H/m)
    R_E = 1560e3;              % Europa radius (m)

    M0 = -(4*pi/mu0) * B_prim .* (R_E^3)/2;
    M_real = real(M0);

    B_prim(3) = -385e-9;

    B_sec = (mu0/(4*pi)) * (3*dot(point, M_real) .* point - ...
    M_real .* norm(point).^2) ./ norm(point)^5;

    B = B_prim + B_sec; 
end