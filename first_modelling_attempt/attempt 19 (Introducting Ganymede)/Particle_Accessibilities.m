function [ ] = Particle_Accessibilities(theta_list, phi_list, escaped, energy, element, grid_len)

% fff

thetas = -90:grid_len:90;
phis = 0:grid_len:359;
Z = zeros(length(thetas), length(phis));  % rows = theta (lat), cols = phi (long)

% Loop order flipped accordingly:
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
            Z(i, j) = count / total * 100;  % now: row = lat, col = long
        else
            Z(i, j) = NaN;
        end
    end
end

i_north = find(thetas == 90);
i_south = find(thetas == -90);

% compute mean accessibility across longitudes at each pole
north_val = mean(Z(i_north, :), 'omitnan');
south_val = mean(Z(i_south, :), 'omitnan');

% assign same value to entire row for θ = ±90°
Z(i_north, :) = north_val;
Z(i_south, :) = south_val;

figure;
clf;
imagesc(phis, thetas, Z);              
set(gca,'YDir','normal'); 
ylabel('Latitude (Degrees)')                              
xlabel('Longitude (Degrees)')           
ax=xticklabels;
xticklabels(flip(ax))
colormap(jet);
cb = colorbar;
cb.Label.String = 'Local Accessibility (%)';
clim([0 100]);
if energy < 1e6
    title(sprintf('Ion Access Map for: %s %.1f KeV', element, energy / 1000));
else
    title(sprintf('Ion Access Map for: %s %.1f MeV', element, energy / (1e6)));
end

% % create finer grids
% fine_theta = linspace(min(thetas), max(thetas), 200);
% fine_phi   = linspace(min(phis), max(phis), 200);
% [fine_theta_grid, fine_phi_grid] = meshgrid(fine_theta, fine_phi);
% [theta_grid, phi_grid] = meshgrid(thetas, phis);
% 
% % interpolate Z to the finer grid
% Z_fine = interp2(theta_grid, phi_grid, Z, fine_theta_grid, fine_phi_grid, 'linear');
% 
% % plot the smoothed heatmap
% figure;
% imagesc(fine_phi, fine_theta, Z_fine);
% set(gca, 'YDir', 'normal');
% axis([min(phis), max(phis), min(thetas), max(thetas)]);
% cb = colorbar;
% cb.Label.String = 'Local Accessibility (%)';
% clim([0 100]);
% ylabel('Latitude (Degrees)');
% xlabel('Longitude (Degrees)');
% ax=xticklabels;
% xticklabels(flip(ax))
% if energy < 1e6
%     title(sprintf('Ion Access Map for: %s %.1f KeV', element, energy / 1000));
% else
%     title(sprintf('Ion Access Map for: %s %.1f MeV', element, energy / (1e6)));
% end