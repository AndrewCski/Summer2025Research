function N_n = get_neutral_vol_density(x, y, z, plume_theta, plume_phi, theta_p, units)

    % Function for retrieving the neutral particle density at a point
    % around Europa. "x", "y", and "z"  are the
    % coordinates of the observer, "plume_theta" is the longitude of the base 
    % of the plume on Europa's surface, "plume_phi" is the plume's
    % longitude, and "theta_p" is the opening angle of the plume, and "units"
    % dictates whther the result is returned in parts / m^3 or per cm^3,
    % chosen by units = "cm", or units == "no exo" or anything else. no exo
    % is included here to calculate the effects of solely the plume.

    R_E = 1560e3;
    n10 = 4e7 * 1e6; hs1 = 100e3;
    n20 = 1e6 * 1e6; hs2 = 500e3;

    r_mag = sqrt(x.^2 + y.^2 + z.^2);
    mask = r_mag >= R_E;
    
    % azimuthal angle, shifted to phi = 0 at -x
    phis = atan2(y, x) + pi;
    phis = mod(phis + pi, 2*pi) - pi;

    baseN = n10 .* exp(-(r_mag - R_E)/hs1) + n20 .* exp(-(r_mag - R_E)/hs2);

    halfPlane = abs(phis) <= pi/2;
    N_n = zeros(size(r_mag));
    N_n(mask) = baseN(mask);
    applyMask = mask & halfPlane;
    N_n(applyMask) = (1 + 2 .* cos(phis(applyMask))) .* ...
        (n10 .* exp(-(r_mag(applyMask) - R_E)/hs1) + ...
        n20 .* exp(-(r_mag(applyMask) - R_E)/hs2));
    N_n(N_n < 0) = 0;

    if ~isnan(plume_theta) && ~isnan(plume_phi)
        plume_phi = deg2rad(plume_phi);
        plume_theta = plume_theta * pi / 180;
        
        plume_r = [R_E * cos(plume_theta) * sin(plume_phi), ...
            R_E * cos(plume_theta) * cos(plume_phi), R_E * sin(plume_theta)];
        
        Np0 = 2e9 * 1e6;
        Hp = 150e3;
        theta_p = theta_p * pi / 180;
        
        % unit surface normal
        r_s = plume_r / norm(plume_r);   % 1x3
        
        % vector from surface point to every grid point
        Vx = x - plume_r(1);
        Vy = y - plume_r(2);
        Vz = z - plume_r(3);
        
        % dot products with r_s 
        dotRV = Vx * r_s(1) + Vy * r_s(2) + Vz * r_s(3);
        
        % norm of V's
        normV = sqrt(Vx.^2 + Vy.^2 + Vz.^2);
        normV(normV < 1e-12) = 1e-12; 

        % polar angles
        cosVals = dotRV ./ normV;
        cosVals = min(max(cosVals, -1), 1);
        thetas = acos(cosVals); 

        plume_indices = thetas < pi/2;
        if units == "no exo"
            for i = 1:numel(N_n)
                N_n(i) = 0;
                if plume_indices(i)
                    N_n(i) = N_n(i) + Np0 * exp((-r_mag(i) + R_E)/Hp) ...
                        * exp(-(thetas(i)/theta_p)^2);
                end
            end
        else
            for i = 1:numel(N_n)
                if plume_indices(i)
                    N_n(i) = N_n(i) + Np0 * exp((-r_mag(i) + R_E)/Hp) ...
                        * exp(-(thetas(i)/theta_p)^2);
                end
            end
        end
    end
    if units == "cm"
        N_n = N_n ./ 1e6;
    end
end