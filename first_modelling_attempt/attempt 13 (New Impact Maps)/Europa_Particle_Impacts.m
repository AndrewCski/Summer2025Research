function [ ] = Europa_Particle_Impacts(escaped, init_r, energies)

% This function file is used for mapping the impact sites of various
% energies of particles across the surface of Europa
%
% The two parameters that must be passed are the impact coordinates of
% the particles, as well as their energies


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
cb = colorbar;
cb.Label.String = 'Particle Energy (eV)';
clim([min(energies) max(energies)])
set(gca,'ColorScale','log')
hold off