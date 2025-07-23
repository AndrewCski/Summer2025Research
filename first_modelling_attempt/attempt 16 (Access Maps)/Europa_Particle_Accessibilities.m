function [ ] = Europa_Particle_Accessibilities(theta_list, phi_list, escaped, energy, element, grid_len)

% fff

thetas = -90:grid_len:90;
phis = 0:grid_len:359;
Z = zeros(length(phis), length(thetas));  % rows = y (phi), cols = x (theta)

for i = 1:length(thetas)
    for j = 1:length(phis)
        count = 0;
        total = 0;
        for k = 1:length(theta_list)
            if theta_list(k) == thetas(i) && phi_list(k) == phis(j)
                if escaped(k)
                    count = count + 1;
                end
                total = total + 1;
            end
        end
        if total > 0
            Z(j, i) = count / total * 100;  % phi is row, theta is col
        else
            Z(j, i) = NaN;
        end
    end
end

% figure;
% clf;
% imagesc(phis, thetas, Z);                      % #5 !!!
% set(gca,'YDir','normal'); % I think this is supposed to be Ydir???
% colorbar;
% ylabel('Latitude (Degrees)')                                                   % #6 !!!
% xlabel('Longitude (Degrees)')                                                   % #7 !!!
% xlim([0, 360])
% ylim([-90, 90])
% ax=xticklabels;
% xticklabels(flip(ax))
% % title(sprintf('Ion Access Map for:', element, energy));
% title(sprintf('Ion Access Map for: %s %.1f keV', element, energy/1000))

% Create finer grids
fine_theta = linspace(min(thetas), max(thetas), 200);
fine_phi   = linspace(min(phis), max(phis), 200);
[fine_theta_grid, fine_phi_grid] = meshgrid(fine_theta, fine_phi);
[theta_grid, phi_grid] = meshgrid(thetas, phis);

% Interpolate Z to the finer grid
Z_fine = interp2(theta_grid, phi_grid, Z, fine_theta_grid, fine_phi_grid, 'linear');

% Plot the smoothed heatmap
figure;
imagesc(fine_phi, fine_theta, Z_fine);
set(gca, 'YDir', 'normal');
axis([min(phis), max(phis), min(thetas), max(thetas)]);
cb = colorbar;
cb.Label.String = 'Local Accessibility (%)';
clim([0 100]);
ylabel('Latitude (Degrees)');
xlabel('Longitude (Degrees)');
ax=xticklabels;
xticklabels(flip(ax))
if energy < 1e6
    title(sprintf('Ion Access Map for: %s %.1f KeV', element, energy / 1000));
else
    title(sprintf('Ion Access Map for: %s %.1f MeV', element, energy / (1e6)));
end