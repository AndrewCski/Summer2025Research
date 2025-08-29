function [ ] = forward_particle_impact_accessibilities(impact_coords, group_width, q_over_m, moon, ...
    element, mode, normalized, fitted)

    % Function for finding how often particle impact certain parts of
    % Europa's/Ganymede's surface. Can be used for simulating plume particles
    % or Jovian system Io dust. impact_coords should be the two last
    % recorded locations of each particle, group_width is the
    % latitude/longitude width of the bins depicting particle impacts on
    % the surface, q_over_m is the charge/mass ratio of particles
    % simulated, moon is a string telling what Jovian moon is being used, 
    % element is what particle type is being used, mode is whether the
    % results should be displaying the element or q/m simulated, and
    % normalized is whether the results should be normalized or not.

if strcmp(moon, "Europa")
    moon_rad = 1560e3;
    moon_string = 'Europa';
elseif strcmp(moon, "Ganymede")
    moon_rad = 2631e3;
    moon_string = 'Ganymede';
end
thetas = -90:group_width:90;
phis = 0:group_width:360;
Z = zeros(length(thetas), length(phis));  % rows = theta (lat), cols = phi (long)

for i = 1:size(impact_coords, 1)
    if ~all(impact_coords(i,1:6) == 0)
        P1 = impact_coords(i,1:3);
        P2 = impact_coords(i,4:6);    
    
        d = P2 - P1;
        A = dot(d, d);
        B = 2 * dot(P1, d);
        C = dot(P1, P1) - moon_rad^2;
    
        t_vals = roots([A, B, C]);
        t_surface = t_vals(t_vals >= 0 & t_vals <= 1);
        p_surface = P1 + t_surface * d;
    
        [phi, theta, ~] = cart2sph(p_surface(1), p_surface(2), p_surface(3));
        phi = rad2deg(phi);
        theta = rad2deg(theta);
        phi(phi < 0) = phi(phi < 0) + 360;
        phi = mod(phi - 90, 360);  % rotate reference
    
        idx_theta = discretize(theta, thetas);
        idx_phi = discretize(phi, phis);
        Z(idx_theta, idx_phi) = Z(idx_theta, idx_phi) + 1;
    end
end

if normalized
    Z = (Z - min(Z, [], 'all')) ./ ...
       (max(Z, [], 'all') - min(Z, [], 'all')) .* 100;
end

figure;
clf;
if fitted
    logical_array = Z > 0.25; % zooming in to where plume particles actually land with some frequency
    rows = any(logical_array, 2);
    cols = any(logical_array, 1);
    Z_2 = Z(rows, :);
    Z_3 = Z_2(:, cols);
    phi_2 = phis(cols);
    theta_2 = thetas(rows);
    imagesc(phi_2, theta_2, Z_3);  
else
    imagesc(phis, thetas, Z);              
end             
set(gca,'YDir','normal'); 
ylabel('Latitude (Degrees)')                              
xlabel('Longitude (Degrees)')           
ax=xticklabels;
xticklabels(flip(ax))
colormap(parula);
cb = colorbar;
if ~normalized
    cb.Label.String = '# of Impacts';
end
clim([0 max(Z(:))]);

if lower(mode) == "q/m" && all(q_over_m == q_over_m(1), 'all')
    title(sprintf('Impact Rates Per %s Surface Positions, q/m = %.3g', moon_string, q_over_m(1)));
elseif lower(mode) == "element" || lower(mode) == "elements"
    title(sprintf('%s Impact Rates Per %s Surface Positions', element, moon_string));
else
    title(sprintf('Impact Rates Per %s Surface Positions', moon_string));

end


