function [ ] = Europa_Particle_Trajectories_Opt(r_cell, mode, parameter, gif)

R_E = 1560e3;               % Europa radius in meters

if gif
    clf;  % clear figure content
end

figure;
hold on
xlabel('x'); ylabel('y'); zlabel('z')
title('Particle Trajectories in Europa Magnetic Environment')

% adding Europa sphere
[xs, ys, zs] = sphere(50);
surf(R_E * xs, R_E * ys, R_E * zs, ...
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
    end

    colormap(cmap)
    cb = colorbar;
    cb.Label.String = 'Particle Energy (eV)';
    clim([min(parameter) max(parameter)])
    set(gca,'ColorScale','log')
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

        if ~legend_added(idx)
            plot_handles(idx) = h;
            legend_added(idx) = true;
        end
    end

    legend(plot_handles, unique_elements, 'Location', 'best');

else
    for i = 1:N_particles
        r_i = r_cell{i};
        plot3(r_i(1, :), r_i(2, :), r_i(3, :), 'LineWidth', 1.5)
    end
end