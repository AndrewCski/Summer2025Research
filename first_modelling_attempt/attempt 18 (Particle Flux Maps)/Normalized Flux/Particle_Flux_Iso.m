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

bin_scale = (resolution - 1) / (2 * sim_size);
len = length(r_cell);
disp("Iso Flux Loop:");

% loop through each particle
for i = 1:len
    % grab particle trajectory over relevant steps
    pos = r_cell{i}(:, first_step:min(last_step, size(r_cell{i},2)));
    % skip if empty
    if isempty(pos)
        continue; 
    end

    % remove positions inside moon
    mask = vecnorm(pos, 2, 1) >= body_radius;
    pos = pos(:, mask);
    if isempty(pos)
        continue; 
    end

    % compute indices directly (faster than discretize)
    idx = floor((pos + sim_size) * bin_scale) + 1;

    % filter valid indices
    valid = all(idx >= 1 & idx <= resolution, 1);
    idx = idx(:, valid);
    if isempty(idx)
        continue; 
    end

    % accumulate into flux_Counts
    lin_idx = sub2ind([resolution resolution resolution], ...
                      idx(1,:), idx(2,:), idx(3,:));
    flux_Counts = flux_Counts + reshape(accumarray(lin_idx.', 1, ...
                       [numel(flux_Counts) 1]), size(flux_Counts));

    if mod(i, 10000) == 0
        fprintf('%d/%d\n', i, len);
    end

end

flux = (flux_Counts - min(flux_Counts, [], 'all'))./(max(flux_Counts, [], 'all') - ...
    min(flux_Counts, [], 'all')) .* 100;

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
    title(sprintf('Normalized Particle Traffic Rates above %.d/100: %s %.1f KeV, t = %.2f to %.2f hours', ...
        bar, element, energy / 1000, (start_time + first_step * abs(time_step))/3600, ...
        (start_time + last_step * abs(time_step))/3600));
else
    title(sprintf('Normalized Particle Traffic Rates above %.d/100: %s %.1f MeV, t = %.2f to %.2f hours', ...
        bar, element, energy / (1e6), (start_time + first_step * abs(time_step))/3600, ...
        (start_time + last_step * abs(time_step))/3600));
end
hold off;