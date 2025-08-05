function [ ] = True_Europa_Particle_Flux(r_cell, first_step, last_step, dims, resolution, ...
    sim_size, energy, element, time_step, start_time)

R_E = 1560e3;  % Europa radius in meters

% define grid along each dimension
sim_range = linspace(-sim_size, sim_size, resolution);

% Select which dimensions to plot (e.g., [1,2] for x vs y)
flag = 0;
Flux_Counts = zeros(resolution, resolution);

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
elseif dims == "all" || dims == "All"
    flag = 1;
    Flux_Counts = zeros(resolution, resolution, 3);
end

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
        if flag
            for k = 1:3
                if k == 1
                    val1 = pos(1);
                    val2 = pos(2);
                elseif k == 2
                    val1 = pos(1);
                    val2 = pos(3);
                else
                    val1 = pos(2);
                    val2 = pos(3);
                end

                % Convert position to indices
                idx1 = find(val1 >= sim_range, 1, 'last');
                idx2 = find(val2 >= sim_range, 1, 'last');
        
                if isempty(idx1) || isempty(idx2) || idx1 >= resolution || idx2 >= resolution
                    continue; % Outside bounds
                end
                Flux_Counts(idx2, idx1, k) = Flux_Counts(idx2, idx1, k) + 1;
            end
        else
            val1 = pos(dim1);
            val2 = pos(dim2);

            % Convert position to indices
            idx1 = find(val1 >= sim_range, 1, 'last');
            idx2 = find(val2 >= sim_range, 1, 'last');
    
            if isempty(idx1) || isempty(idx2) || idx1 >= resolution || idx2 >= resolution
                continue; % Outside bounds
            end
            Flux_Counts(idx2, idx1) = Flux_Counts(idx2, idx1) + 1;
        end
    end
end

area_per_bin = (2*sim_size / resolution)^2;
Flux = Flux_Counts ./ ((last_step - first_step + 1) * area_per_bin);

% Plot

if flag
    for i = 1:3
        figure;
        imagesc(sim_range/R_E, sim_range/R_E, Flux(:,:,i));
        axis xy;

        if i == 1
            xlabel('X (R_E)');
            ylabel('Y (R_E)');
        elseif i == 2
            xlabel('X (R_E)');
            ylabel('Z (R_E)');
        else
            xlabel('Y (R_E)');
            ylabel('Z (R_E)');
        end

        if energy < 1e6
            title(sprintf('Particle Flux (parts/m^3/s) for: %s %.1f KeV, from t = %.2f to %.2f hours', ...
                element, energy / 1000, (start_time + first_step * abs(time_step))/3600, ...
                (start_time + last_step * abs(time_step))/3600));
        else
            title(sprintf('Particle Flux (parts/m^3/s) for: %s %.1f MeV, from t = %.2f to %.2f hours', ...
                element, energy / (1e6), (start_time + first_step * abs(time_step))/3600, ...
                (start_time + last_step * abs(time_step))/3600));
        end
        cb = colorbar;
        cb.Label.String = 'Counts/Second';
    end
else
    figure;
    imagesc(sim_range/R_E, sim_range/R_E, Flux);
    axis xy;
    xlabel([dims(1), ' (R_E)']);
    ylabel([dims(2), ' (R_E)']);
    if energy < 1e6
        title(sprintf('Particle Flux (parts/m^3/s) for: %s %.1f KeV, from t = %.2f to %.2f hours', ...
            element, energy / 1000, (start_time + first_step * abs(time_step))/3600, ...
            (start_time + last_step * abs(time_step))/3600));
    else
        title(sprintf('Particle Flux (parts/m^3/s) for: %s %.1f MeV, from t = %.2f to %.2f hours', ...
            element, energy / (1e6), (start_time + first_step * abs(time_step))/3600, ...
            (start_time + last_step * abs(time_step))/3600));
    end
    cb = colorbar;
    cb.Label.String = 'Counts/Second';
end