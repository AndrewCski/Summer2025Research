function N_c = get_neutral_col_density(point, plume_theta, plume_phi, theta_p)
    R_E = 1560e3;     % Europa radius (m)

    normal_vec = point / norm(point);
    ds = 1e4;  % step size in meters
    s_vals = R_E:ds:norm(point);  % distance along line of sight
    points = s_vals.' .* normal_vec;
    densities = get_neutral_vol_density(points(:,1), points(:,2), points(:,3), ...
        plume_theta, plume_phi, theta_p, "m");

    N_c = trapz(s_vals, densities);
end