function N_c = get_electron_col_density(point, plume_theta, plume_phi, theta_p)

    % Function for retrieving the electron column density from a point
    % around Europa, returned in parts / m^2. "point" is the vector for the
    % position of the observer, "plume_theta" is the longitude of the base 
    % of the plume on Europa's surface, "plume_phi" is the plume's
    % longitude, and "theta_p" is the opening angle of the plume.

    R_E = 1560e3;     % Europa radius (m)

    normal_vec = point / norm(point);
    if sum(point.^2) >= (1.5*R_E)^2
        s_vals = linspace(R_E, norm(point), 700);
    else
        ds = 1e4;  % step size in meters
        s_vals = R_E:ds:norm(point);  % distance along line of sight
    end

    points = s_vals.' .* normal_vec;
    neutral_densities = get_neutral_vol_density(points(:,1), points(:,2), points(:,3), ...
        plume_theta, plume_phi, theta_p, "m");

    O_ion_rates = neutral_densities .* 1e-6;
    O2_ion_rates = neutral_densities .* 1e-7;
    taus = 225 + 1300 * (1 - exp((-vecnorm(points, 2, 2) + R_E)/600e3));

    electron_densities = (O_ion_rates + O2_ion_rates) .* taus;

    N_c = trapz(s_vals, electron_densities);
end