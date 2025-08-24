function [ ] = Particle_Flux_2D(r_cell, first_step, last_step, dims, resolution, ...
    sim_size, energy, element, time_step, start_time, body, update_gap)

    % Function for creating 2d slice plots depicting the normalized
    % rates at which particles simulated pass through binned regions around 
    % a moon. r_cell is the positions of the particles, first_step is the first
    % timestep to be represented here, last_step the last timestep, dims
    % should be a string stating which axes the plot (or plots) should be
    % plotted on, resolution is how many bins per axis should be created,
    % sim_size is the distance from the origin that the faces of the
    % simulation box reside at, energy is the energy of the simulated
    % particles, element is the type of particle simulated, time_step is
    % the time step used in the simulation, start_time is the point in the
    % synodic period that the simulation is started in, body is the moon
    % that is being used, and update_gap is the number of particles that
    % should be simulated per each update of plotting progress.

if body == "Europa"
    body_radius = 1560e3;
    radius_string = "(R_E)";
elseif body == "Ganymede"
    body_radius = 2631e3;
    radius_string = "(R_G)";
end

% define grid and scaling constants
sim_range = linspace(-sim_size, sim_size, resolution);
bin_scale = (resolution - 1) / (2 * sim_size);

% set which dimensions to use
flag = 0;
flux_Counts = zeros(resolution, resolution);

if dims == "XY"
    dim1 = 1; dim2 = 2;
elseif dims == "YX"
    dim1 = 2; dim2 = 1; 
elseif dims == "XZ"
    dim1 = 1; dim2 = 3;
elseif dims == "ZX"
    dim1 = 3; dim2 = 1;
elseif dims == "YZ"
    dim1 = 2; dim2 = 3;
elseif dims == "ZY"
    dim1 = 3; dim2 = 2;
elseif dims == "all" || dims == "All"
    flag = 1;
    flux_Counts = zeros(resolution, resolution, 3);
end

len = length(r_cell);

disp("2D Flux Loop:");

% Loop through each particle
for i = 1:len
    % extract positions across steps
    pos = r_cell{i}(:, first_step:min(last_step, size(r_cell{i},2)));
    if isempty(pos)
        continue; 
    end

    % skip points inside body
    mask = vecnorm(pos, 2, 1) >= body_radius;
    pos = pos(:, mask);
    if isempty(pos)
        continue; 
    end

    if flag
        % XY
        idx1 = floor((pos(1,:) + sim_size) * bin_scale) + 1;
        idx2 = floor((pos(2,:) + sim_size) * bin_scale) + 1;
        valid = idx1 >= 1 & idx1 <= resolution & idx2 >= 1 & idx2 <= resolution;
        lin_idx = sub2ind([resolution resolution], idx2(valid), idx1(valid));
        flux_Counts(:,:,1) = flux_Counts(:,:,1) + reshape(accumarray(lin_idx.',1,...
            [resolution*resolution 1]),resolution,resolution);

        % XZ
        idx1 = floor((pos(1,:) + sim_size) * bin_scale) + 1;
        idx2 = floor((pos(3,:) + sim_size) * bin_scale) + 1;
        valid = idx1 >= 1 & idx1 <= resolution & idx2 >= 1 & idx2 <= resolution;
        lin_idx = sub2ind([resolution resolution], idx2(valid), idx1(valid));
        flux_Counts(:,:,2) = flux_Counts(:,:,2) + reshape(accumarray(lin_idx.',1,...
            [resolution*resolution 1]),resolution,resolution);

        % YZ
        idx1 = floor((pos(2,:) + sim_size) * bin_scale) + 1;
        idx2 = floor((pos(3,:) + sim_size) * bin_scale) + 1;
        valid = idx1 >= 1 & idx1 <= resolution & idx2 >= 1 & idx2 <= resolution;
        lin_idx = sub2ind([resolution resolution], idx2(valid), idx1(valid));
        flux_Counts(:,:,3) = flux_Counts(:,:,3) + reshape(accumarray(lin_idx.',1,...
            [resolution*resolution 1]),resolution,resolution);

    else
        % chosen 2D slice
        val1 = pos(dim1,:);
        val2 = pos(dim2,:);
        idx1 = floor((val1 + sim_size) * bin_scale) + 1;
        idx2 = floor((val2 + sim_size) * bin_scale) + 1;
        valid = idx1 >= 1 & idx1 <= resolution & idx2 >= 1 & idx2 <= resolution;
        lin_idx = sub2ind([resolution resolution], idx2(valid), idx1(valid));
        flux_Counts = flux_Counts + reshape(accumarray(lin_idx.',1,[resolution*resolution 1]),...
            resolution,resolution);
    end

    if mod(i, update_gap) == 0
        fprintf('%d/%d\n', i, len);
    end

end

% normalize flux
flux = (flux_Counts - min(flux_Counts, [], 'all')) ./ ...
       (max(flux_Counts, [], 'all') - min(flux_Counts, [], 'all')) .* 100;

if flag
    for i = 1:3
        figure;
        imagesc(sim_range/body_radius, sim_range/body_radius, flux(:,:,i));
        axis xy;
        if i == 1
            xlabel(sprintf('X %s', radius_string));
            ylabel(sprintf('Y %s', radius_string));
        elseif i == 2
            xlabel(sprintf('X %s', radius_string));
            ylabel(sprintf('Z %s', radius_string));
        else
            xlabel(sprintf('Y %s', radius_string));
            ylabel(sprintf('Z %s', radius_string));
        end

        if energy < 1e6
            title(sprintf('Particle Traffic Rates Per Region: %s %.1f KeV, t = %.2f to %.2f hours', ...
                element, energy / 1000, (start_time + first_step * abs(time_step))/3600, ...
                (start_time + last_step * abs(time_step))/3600));
        else
            title(sprintf('Particle Traffic Rates Per Region: %s %.1f MeV, t = %.2f to %.2f hours', ...
                element, energy / (1e6), (start_time + first_step * abs(time_step))/3600, ...
                (start_time + last_step * abs(time_step))/3600));
        end
        cb = colorbar; clim([0 100])
    end
else
    figure;
    imagesc(sim_range/body_radius, sim_range/body_radius, flux);
    axis xy;
    xlabel([dims(1), ' ', radius_string]);
    ylabel([dims(2), ' ', radius_string]);
    if energy < 1e6
        title(sprintf('Normalized Particle Traffic Rates: %s %.1f KeV, t = %.2f to %.2f hours', ...
            element, energy / 1000, (start_time + first_step * abs(time_step))/3600, ...
            (start_time + last_step * abs(time_step))/3600));
    else
        title(sprintf('Normalized Particle Traffic Rates: %s %.1f MeV, t = %.2f to %.2f hours', ...
            element, energy / (1e6), (start_time + first_step * abs(time_step))/3600, ...
            (start_time + last_step * abs(time_step))/3600));
    end
    cb = colorbar; clim([0 100])
end