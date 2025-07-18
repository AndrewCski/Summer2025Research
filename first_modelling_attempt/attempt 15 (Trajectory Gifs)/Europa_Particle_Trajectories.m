function [ ] = Europa_Particle_Trajectories(r, num_recorded, mode, parameter)

% This function file is used for mapping the trajectories of various
% particles around Europa.
%
% The two parameters that must be passed are the positions of the particles
% at every timestep (r), and the number of recorded timesteps per particles
% (num recorded). From there, there are three modes for this function.
% 'energies' makes the plot disciminate particles via their starting
% energies, 'elements' does the same but for the elements of the particles,
% and putting nothing results in no identification of particles. the
% parameter named 'parameter' should be the energies, element names, or
% nothing, respectively.

R_E = 1560e3;               % Europa radius in meters

clf;  % clear figure content
figure;
hold on
xlabel('x'); ylabel('y'); zlabel('z')
title('Particle Trajectories in Europa Magnetic Environment')

% adding Europa sphere
[xs, ys, zs] = sphere(50);  % higher resolution sphere
surf(R_E * xs, R_E * ys, R_E * zs, ...
    'FaceColor', [0.2 0.5 1], ...     % easy-ish to see blue
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.5);               % semi-transparent

% set view and axis
axis equal
grid on
view(3)
camzoom(2)

if mode == "energies"
    log_energies = log10(parameter);
    minE = min(log_energies);
    maxE = max(log_energies);
    norm_energies = (log_energies - minE) / (maxE - minE);
    cmap = jet(256);
    
    for i = 1:size(r, 3)
        cidx = round(1 + norm_energies(i) * (size(cmap,1)-1));
        color = cmap(cidx, :);
        plot3(r(1, 1:num_recorded(i), i), r(2, 1:num_recorded(i), i), ...
            r(3, 1:num_recorded(i), i), 'Color', color, 'LineWidth', 1.5)
    end
    
    colormap(cmap)
    cb = colorbar;
    cb.Label.String = 'Particle Energy (eV)';
    clim([min(parameter) max(parameter)])
    set(gca,'ColorScale','log')
    hold off
elseif mode == "elements"
    % Make sure 'parameter' is a cell array of strings or categorical
    elements = string(parameter);  % or use cellstr if parameter is cell
    
    % get unique elements and assign each one a color
    unique_elements = unique(elements, 'stable');  % keep order if needed
    cmap = lines(numel(unique_elements));
    
    % store plot handles for legend
    plot_handles = gobjects(numel(unique_elements),1);
    legend_added = false(numel(unique_elements),1);
    
    for i = 1:size(r, 3)
        this_element = elements(i);
        % find index of this element in unique list
        idx = find(unique_elements == this_element, 1);
        color = cmap(idx, :);
        
        h = plot3(r(1, 1:num_recorded(i), i), ...
                  r(2, 1:num_recorded(i), i), ...
                  r(3, 1:num_recorded(i), i), ...
                  'Color', color, 'LineWidth', 1.5);
        
        % save the first plotted handle for each unique element
        if ~legend_added(idx)
            plot_handles(idx) = h;
            legend_added(idx) = true;
        end
    end
    legend(plot_handles, unique_elements, 'Location', 'best')  
else
    for i = 1:size(r, 3)
        plot3(r(1, 1:num_recorded(i), i), r(2, 1:num_recorded(i), i), ...
            r(3, 1:num_recorded(i), i), 'LineWidth', 1.5)
    end
end