function N_c = get_neutral_col_density(point, plume_theta, plume_phi, theta_p)
    R_E = 1560e3;     % Europa radius (m)

    plume_phi = mod(plume_phi * pi / 180, 2*pi) - pi;
    plume_theta = plume_theta * pi / 180;

    normal_vec = point / norm(point);
    ds = 1e4;  % step size in meters
    s_vals = R_E:ds:norm(point);  % distance along line of sight
    points = s_vals.' .* normal_vec;
    densities = get_neutral_vol_density(points(:,1), points(:,2), points(:,3), ...
        plume_theta, plume_phi, theta_p);

    N_c = trapz(s_vals, densities);
end