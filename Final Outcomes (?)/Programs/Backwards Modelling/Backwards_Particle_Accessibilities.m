function [ ] = Backwards_Particle_Accessibilities(theta_list, phi_list, escaped, energy, element, grid_len)

    % Function for plotting the accessibilities of each point on the
    % surface of a moon when conducting backwards modelling of surface
    % impacting particles. theta_list is the list of latitudes of initial 
    % positions of particles on Europa's surface, phi_list is their
    % longitudes, escaped is an array indicating whether each particle
    % reimpacted onto Europa's surface or came from beyond the moon, energy
    % is the impact energies of particles simulated, element is the type of
    % particle simulated, and grid_len is the width in latitude/longitude
    % of each bin of impact positions of particles.

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

% assign same value to entire row for theta = ±90°
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