function [ ] = True_Europa_Particle_Flux(r_cell, first_step, last_step, dims, resolution, ...
    sim_size, energy, element, time_step, start_time)

R_E = 1560e3;  % Europa radius in meters

% define grid along each dimension
x_range = linspace(-sim_size, sim_size, resolution);
y_range = linspace(-sim_size, sim_size, resolution);
z_range = linspace(-sim_size, sim_size, resolution);

% Select which dimensions to plot (e.g., [1,2] for x vs y)
if dims == "XY"
    dim1 = 1;
    dim2 = 2;
elseif dims == "YX"
    dim1 = 2;
    dim2 = 1; 
elseif dims == "XZ"
    dim1 = 1;
    dim2 = 3;
elseif dims == "ZX"
    dim1 = 3;
    dim2 = 1;
elseif dims == "YZ"
    dim1 = 2;
    dim2 = 3;
elseif dims == "ZY"
    dim1 = 3;
    dim2 = 2;
end

% Create empty flux matrix
Flux_Counts = zeros(resolution, resolution);

% Loop through each particle
for i = 1:length(r_cell)
    % Store which bins this particle has entered to avoid double-counting

    % Loop through timesteps
    for j = first_step:last_step

        if j > length(r_cell{i})
            break;
        end

        % Get position of particle at timestep j
        pos = r_cell{i}(:, j);  % 3x1 vector

        % Skip if inside Europa
        if norm(pos) < R_E
            continue;
        end

        % Extract only the two dimensions of interest
        val1 = pos(dim1);
        val2 = pos(dim2);

        % Convert position to indices
        idx1 = find(val1 >= x_range, 1, 'last');
        idx2 = find(val2 >= y_range, 1, 'last');

        if isempty(idx1) || isempty(idx2) || idx1 >= resolution || idx2 >= resolution
            continue; % Outside bounds
        end

        % Avoid counting same particle multiple times in the same bin
        Flux_Counts(idx2, idx1) = Flux_Counts(idx2, idx1) + 1;
    end
end

% area_per_bin = (2*sim_size / resolution)^2;
% Flux = Flux_Counts / ((last_step - first_step + 1) * area_per_bin);
Flux = Flux_Counts / (last_step - first_step + 1); % I dont know if I like using area

% Plot
figure;
imagesc(x_range/R_E, y_range/R_E, Flux);
axis xy;
xlabel([dims(1), ' (R_E)']);
ylabel([dims(2), ' (R_E)']);
if energy < 1e6
    title(sprintf('Particle Flux for: %s %.1f KeV, from t = %.2f to %.2f hours', ...
        element, energy / 1000, (start_time + first_step * abs(time_step))/3600, ...
        (start_time + last_step * abs(time_step))/3600));
else
    title(sprintf('Particle Flux for: %s %.1f MeV, from t = %.2f to %.2f hours', ...
        element, energy / (1e6), (start_time + first_step * abs(time_step))/3600, ...
        (start_time + last_step * abs(time_step))/3600));
end
cb = colorbar;
cb.Label.String = 'Counts/Second';
end