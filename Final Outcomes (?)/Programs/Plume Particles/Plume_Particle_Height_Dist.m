function [] = Plume_Particle_Height_Dist(r_cell, masses, element, p_theta, p_phi)

    % Function for plotting the maximum altitudes achieved by every particle. 
    % r_cell is each of the positions at each step of each particle, masses
    % is the masses of each particle, element is the type of particle
    % simulated, p_theta is the latitude of the base of the plume on
    % Europa's surface, and p_phi is the longitude of the plume.

    R_E = 1560e3;                                         % Europa radius (m)

    max_altitudes = zeros(length(masses),1);
    
    for p_index = 1:length(masses)
        h_max = max(sqrt(sum(r_cell{p_index}.^2,1)) - R_E);
        max_altitudes(p_index) = h_max;
    end
    
    % Make a histogram of all altitudes
    figure;
    histogram(max_altitudes, 50, 'Normalization', 'percentage'); % in km
    xlabel('Altitude (m)');
    ylabel('Percentage of Particles');
    title(sprintf('Height distribution of Europa Plume %s (θ = %d, Φ = %d)', element, p_theta, p_phi))
end