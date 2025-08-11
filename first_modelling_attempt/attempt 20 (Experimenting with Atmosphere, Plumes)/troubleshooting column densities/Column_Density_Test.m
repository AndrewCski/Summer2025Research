function N_n = get_neutral_vol_density(x, y, z, plume_info)
    R_E = 1560e3;
    n10 = 4e7; hs1 = 100e3;
    n20 = 1e6; hs2 = 500e3;

    r_mag = sqrt(x.^2 + y.^2 + z.^2);
    mask = r_mag >= R_E;
    
    % azimuthal angle, shifted to phi = 0 at -x
    phis = atan2(y, x) + pi;
    phis = mod(phis + pi, 2*pi) - pi;

    baseN = n10 .* exp(-(r_mag - R_E)/hs1) + n20 .* exp(-(r_mag - R_E)/hs2);

    halfPlane = abs(phis) <= pi/2;
    N_n = NaN(size(r_mag));
    N_n(mask) = baseN(mask);
    applyMask = mask & halfPlane;
    N_n(applyMask) = (1 + 2 .* cos(phis(applyMask))) .* ...
        (n10 .* exp(-(r_mag(applyMask) - R_E)/hs1) + ...
        n20 .* exp(-(r_mag(applyMask) - R_E)/hs2));
    N_n(N_n < 0) = 0;

end

function [N_c, points, densities] = get_neutral_col_density(point, plume_info)
    R_E = 1560e3;     % Europa radius (m)

    normal_vec = point / norm(point);
    ds = 1e4;  % step size in meters
    s_vals = R_E:ds:norm(point);  % distance along line of sight
    points = s_vals.' .* normal_vec;
    densities = get_neutral_vol_density(points(:,1), points(:,2), points(:,3), plume_info);

    N_c = trapz(s_vals, densities);
end

R_E = 1560e3;     % Europa radius (m)
obs_pos = [2.5*R_E, 0, 0];

[col, points, dens] = get_neutral_col_density(obs_pos, 0);
figure;
plot(points/R_E, dens)
xlabel('radius (R_E)')
ylabel('neutral density (cm^-3)')