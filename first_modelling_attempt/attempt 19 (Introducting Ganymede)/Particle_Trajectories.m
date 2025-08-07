function [ ] = Particle_Trajectories(r_cell, mode, parameter, gif, body)

if body == "Europa"
    body_radius = 1560e3;
    title_string = "Particle Trajectories in Europa Magnetic Environment";
elseif body == "Ganymede"
    body_radius = 2631e3;
    title_string = "Particle Trajectories in Ganymede Magnetic Environment";
end
 
if gif
    clf;  % clear figure content
end

figure;
hold on
xlabel('x'); ylabel('y'); zlabel('z')
title(title_string)

% adding moon sphere
[xs, ys, zs] = sphere(50);
surf(body_radius * xs, body_radius * ys, body_radius * zs, ...
    'FaceColor', [0.2 0.5 1], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.5);

axis equal
grid on
view(3)
if gif
    camzoom(2)
end

N_particles = numel(r_cell);

if mode == "energies"
    log_energies = log10(parameter);
    minE = min(log_energies);
    maxE = max(log_energies);
    norm_energies = (log_energies - minE) / (maxE - minE);
    cmap = jet(256);

    for i = 1:N_particles
        r_i = r_cell{i};  % 3 x num_steps
        cidx = round(1 + norm_energies(i) * (size(cmap,1)-1));
        color = cmap(cidx, :);
        plot3(r_i(1, :), r_i(2, :), r_i(3, :), 'Color', color, 'LineWidth', 1.5)
        plot3(r_i(1, length(r_i)), r_i(2, length(r_i)), r_i(3, length(r_i)), 'ro')
    end
    if length(parameter) > 1
        colormap(cmap)
        cb = colorbar;
        cb.Label.String = 'Particle Energy (eV)';
        clim([min(parameter) max(parameter)])
        set(gca,'ColorScale','log')
    end
    hold off

elseif mode == "elements"
    elements = string(parameter);
    unique_elements = unique(elements, 'stable');
    cmap = lines(numel(unique_elements));
    plot_handles = gobjects(numel(unique_elements),1);
    legend_added = false(numel(unique_elements),1);

    for i = 1:N_particles
        r_i = r_cell{i};
        this_element = elements(i);
        idx = find(unique_elements == this_element, 1);
        color = cmap(idx, :);

        h = plot3(r_i(1, :), r_i(2, :), r_i(3, :), ...
                  'Color', color, 'LineWidth', 1.5);
        plot3(r_i(1, size(r_i, 2)), r_i(2, size(r_i, 2)), r_i(3, size(r_i, 2)), 'ro')

        if ~legend_added(idx)
            plot_handles(idx) = h;
            legend_added(idx) = true;
        end
    end

else
    for i = 1:N_particles
        if length(r_cell{i}) > 40000
            r_i = r_cell{i};
            plot3(r_i(1, :), r_i(2, :), r_i(3, :), 'LineWidth', 1.5)
            plot3(r_i(1, length(r_i)), r_i(2, length(r_i)), r_i(3, length(r_i)), 'ro')
        end
    end
end