function [ ] = Europa_Particle_Trajectories(r, num_recorded, energies)

% This function file is used for mapping the trajectories of various
% particles around Europa.
%
% The two parameters that must be passed are the positions of the particles
% at every timestep (r), and the number of recorded timesteps per particles
% (num recorded). Other parameters, such as energy as of now, can be chosen
% to be included or not. Pass a [] instead of a variable when calling the
% function to not utilize particle energies, or whatever other variable
% eventually.

R_E = 1560e3;               % Europa radius in meters
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

if energies
    log_energies = log10(energies);
    minE = min(log_energies);
    maxE = max(log_energies);
    norm_energies = (log_energies - minE) / (maxE - minE);
    cmap = jet(256);
    
    for i = 1:size(r, 3)
        cidx = round(1 + norm_energies(i) * (size(cmap,1)-1));
        color = cmap(cidx, :);
        plot3(r(1, 1:num_recorded(i), i), r(2, 1:num_recorded(i), i), ...
            r(3, 1:num_recorded(i), i), 'Color', color)
    end
    
    colormap(cmap)
    cb = colorbar;
    cb.Label.String = 'Particle Energy [eV]';
    clim([min(energies) max(energies)])
    set(gca,'ColorScale','log')
    hold off
else
    for i = 1:size(r, 3)
        plot3(r(1, 1:num_recorded(i), i), r(2, 1:num_recorded(i), i), ...
            r(3, 1:num_recorded(i), i))
    end
end
