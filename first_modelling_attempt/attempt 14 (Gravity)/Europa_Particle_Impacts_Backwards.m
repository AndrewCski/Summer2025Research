function [ ] = Europa_Particle_Impacts_Backwards(escaped, init_r, energies)

% This function file is used for mapping the impact sites of various
% energies of particles across the surface of Europa, when tracing
% particles back in time. 
%
% The three parameters to pass are an array of 1's and 0's, where a 1
% indicates that a corresponding particle came from beyond Europa and a 0
% indicates otherwise, an array of the initial positions of the particles
% in the simulation, and an array of the starting energies of each of the
% particles.


R_E = 1560e3;                          % Europa radius in meters

figure;
hold on
xlabel('x'); ylabel('y'); zlabel('z')
title('Particle Landing Sites')
[xs, ys, zs] = sphere(50);             % higher resolution sphere
surf(R_E * xs, R_E * ys, R_E * zs, ...
    'FaceColor', [0.2 0.5 1], ...      % easy-ish to see blue
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.5);                 % semi-transparent

% set view and axis
axis equal
grid on
view(3)


for i = 1:size(escaped, 1)
    if escaped(i) == 1
        scatter3(init_r(i, 1), init_r(i, 2), init_r(i, 3), ...
            [], energies(i), 'filled');
    end
end

Emin = min(energies(escaped == 1));
Emax = max(energies(escaped == 1));

if Emin ~= Emax
    cb = colorbar;
    cb.Label.String = 'Particle Energy (eV)';
    clim([Emin Emax])
    set(gca,'ColorScale','log')
end

hold off

% plot 2D projections

escaped_indices = escaped == 1;
escaped_r = init_r(escaped_indices, :);
escaped_energies = energies(escaped_indices);

% get unique energies
unique_energies = unique(escaped_energies);
num_energies = length(unique_energies);

for k = 1:num_energies
    E = unique_energies(k);

    % get indices for this energy
    idx = escaped_energies == E;
    coords = escaped_r(idx, :);
    E_coords = coords ./ R_E;
    
    figure;
    
    % x-y projection
    subplot(1, 3, 1);
    scatter(E_coords(:,1), E_coords(:,2), 10, E * ones(size(E_coords,1),1), 'filled');
    xlabel('x (R_E)'); ylabel('y (R_E)');

    % x-z projection
    subplot(1, 3, 2);
    scatter(E_coords(:,1), E_coords(:,3), 10, E * ones(size(E_coords,1),1), 'filled');
    xlabel('x (R_E)'); ylabel('z (R_E)');

    % y-z projection
    subplot(1, 3, 3);
    scatter(E_coords(:,2), E_coords(:,3), 10, E * ones(size(E_coords,1),1), 'filled');
    xlabel('y (R_E)'); ylabel('z (R_E)');

    sgtitle(['Particle Impact Coords, E = ', num2str(E, '%.2e'), ' eV'])
end