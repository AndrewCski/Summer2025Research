function [ ] = forward_particle_impact_accessibilities(impact_coords, group_width, q_over_m)

R_E = 1560e3;                                         % Europa radius (m)
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
        C = dot(P1, P1) - R_E^2;
    
        t_vals = roots([A, B, C]);
        t_surface = t_vals(t_vals >= 0 & t_vals <= 1);
        p_surface = P1 + t_surface * d;
    
        [phi, theta, ~] = cart2sph(p_surface(1), p_surface(2), p_surface(3));
        phi = rad2deg(phi);
        theta = rad2deg(theta);
        phi(phi < 0) = phi(phi < 0) + 360;
        phi = mod(phi + 90, 360);  % rotate reference
    
        idx_theta = discretize(theta, thetas);
        idx_phi = discretize(phi, phis);
        Z(idx_theta, idx_phi) = Z(idx_theta, idx_phi) + 1;
    end
end

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
cb.Label.String = '# of Impacts';
clim([0 max(Z(:))]);
title(sprintf('Dust Impact Rates Per Moon Positions, with q/m = %d', q_over_m))
