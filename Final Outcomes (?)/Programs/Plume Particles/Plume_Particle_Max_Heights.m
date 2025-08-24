function [ ] = Plume_Particle_Max_Heights(r_cell, mass_range, speed_range, simulated_masses, ...
    simulated_speeds, mode, p_theta, p_phi)

    % Function for plotting the maximum altitudes achieved by particles of
    % each configuration. r_cell is each of the positions at each step of
    % each particle, mass_range is the masses used in the simulation,
    % speed_range is the initial speeds used in the simulation,
    % simulated_masses are the masses of each particles, simulated_speeds
    % are the inital speeds of each particle, mode is whether dust or ions
    % are simulated (plot relative to energy or radius), p_theta is the
    % latitude of the base of the plume on Europa's surface, and p_phi is
    % the plume's longitude.

R_E = 1560e3;  % Europa radius (m)

if lower(mode) == "dust" || lower(mode) == "ice"

    radius_range = (3 .* mass_range ./ (4 * pi * 920)).^(1/3);
    max_dists = zeros(length(radius_range), 1);
    simulated_radii = (3 .* simulated_masses ./ (4 * pi * 920)).^(1/3);

    for i = 1:length(radius_range)
        this_radius = radius_range(i);

        idx = find(abs(simulated_radii - this_radius) < 1e-12);

        if isempty(idx)
            max_dists(i) = NaN; % no particles of this radius
            continue;
        end

        dist_vals = zeros(length(idx),1);
        for j = 1:length(idx)
            traj = r_cell{idx(j)};  % 3 × y trajectory matrix
            dists = sqrt(sum(traj.^2,1));  % norm at each timestep
            dist_vals(j) = max(dists);     % max distance this particle achieved
        end

        max_dists(i) = max(dist_vals);
    end
    max_dists = max_dists - R_E;

    figure;
    loglog(radius_range, max_dists)
    grid on;
    xlabel('Particle Radii (m)')
    ylabel('Height from Europa Surface (m)')
    title(sprintf('Max Heights of Europa Plume Ice (θ = %d, Φ = %d)', p_theta, p_phi))
elseif lower(mode) == "ion" || lower(mode) == "ions"
    c = 299792458;                                        % speed of light (m/s)
    max_dists = zeros(length(speed_range), 1);
    for i = 1:length(speed_range)
        this_speed = speed_range(i);

        idx = find(abs(simulated_speeds - this_speed) < 1e-12);

        if isempty(idx)
            max_dists(i) = NaN; % no particles of this radius
            continue;
        end

        dist_vals = zeros(length(idx),1);
        for j = 1:length(idx)
            traj = r_cell{idx(j)};  % 3 × y trajectory matrix
            dists = sqrt(sum(traj.^2,1));  % norm at each timestep
            dist_vals(j) = max(dists);     % max distance this particle achieved  dist_vals(j)  
        end

        max_dists(i) = max(dist_vals);
    end      
    max_dists = max_dists - R_E;

    energy_range = 6.242e18 .* (mass_range .* c.^2 ./ sqrt(1 - (speed_range./c).^2) - ...
        mass_range .* c.^2);
    figure;
    plot(energy_range, max_dists)
    set(gca, 'XScale', 'log');
    grid on;
    xlabel('Initial Ion Energies (eV)')
    ylabel('Height from Europa Surface (m)')
    title(sprintf('Max Heights of Europa Plume Ions (θ = %d, Φ = %d)', p_theta, p_phi))
end