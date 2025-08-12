function [ ] = Particle_Flux_Iso(r_cell, first_step, last_step, resolution, ...
    sim_size, energy, element, time_step, start_time, body, bar)

if body == "Europa"
    body_radius = 1560e3;
    radius_string = "(R_E)";
elseif body == "Ganymede"
    body_radius = 2631e3;
    radius_string = "(R_G)";
end

% define grid along each dimension
sim_range = linspace(-sim_size, sim_size, resolution);
flux_Counts = zeros(resolution, resolution, resolution);

% loop through each particle
for i = 1:length(r_cell)
    % store which bins this particle has entered to avoid double-counting
    % loop through timesteps
    for j = first_step:last_step
        if j > size(r_cell{i}, 2)
            break;
        end
        pos = r_cell{i}(:, j);

        % skip if inside moon
        if norm(pos) < body_radius
            continue;
        end

        val1 = pos(1);
        val2 = pos(2);
        val3 = pos(3);
        idx1 = discretize(val1, sim_range);
        idx2 = discretize(val2, sim_range);
        idx3 = discretize(val3, sim_range);

        if isempty(idx1) || isempty(idx2) || isempty(idx3) || ...
            idx1 > resolution || idx2 > resolution || idx3 > resolution
            continue; % outside bounds
        end
        flux_Counts(idx1, idx2, idx3) = flux_Counts(idx1, idx2, idx3) + 1;
    end
end

vol_per_bin = (2*sim_size / resolution)^3;
flux = flux_Counts ./ ((last_step - first_step + 1) * vol_per_bin);
sprintf("max flux: %d", max(flux(:)))
sprintf("min flux: %d", min(flux(:)))

figure;

% create sphere
[XS, YS, ZS] = sphere(50); 

% plot isosurface
[X, Y, Z] = meshgrid(sim_range/body_radius, sim_range/body_radius, sim_range/body_radius);
p1 = patch(isosurface(X, Y, Z, flux, bar));
p1.FaceColor = 'red';
p1.EdgeColor = 'none';
p1.FaceAlpha = 0.75; 

hold on;

% plot moon
p2 = surf(XS, YS, ZS);
set(p2, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.25);
axis equal;
xlabel(sprintf('x %s', radius_string));
ylabel(sprintf('y %s', radius_string));
zlabel(sprintf('z %s', radius_string));
if energy < 1e6
    title(sprintf('Particle Flux above %.d m^{-3}s^{-1} for: %s %.1f KeV, from t = %.2f to %.2f hours', ...
        bar, element, energy / 1000, (start_time + first_step * abs(time_step))/3600, ...
        (start_time + last_step * abs(time_step))/3600));
else
    title(sprintf('Particle Flux above %.d m^{-3}s^{-1} for: %s %.1f MeV, from t = %.2f to %.2f hours', ...
        bar, element, energy / (1e6), (start_time + first_step * abs(time_step))/3600, ...
        (start_time + last_step * abs(time_step))/3600));
end
hold off;