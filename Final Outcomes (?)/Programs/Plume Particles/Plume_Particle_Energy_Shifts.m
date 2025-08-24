function [] = Plume_Particle_Energy_Shifts(v_cell, masses, element, p_theta, p_phi)

    % Function for plotting the change in energy (in eV) of particles from
    % their starting states of leaving a plume and ending states of hitting
    % Europa's surface. v_cell is the velocities at each step of every
    % particles, masses are the masses of each particle, element is the
    % type of particle simulated, p_theta is the latitude of the plume's
    % base on Europa's surface, and p_phi is the plume's longitude.

    c = 299792458;                                        % speed of light (m/s)

    deltaE = zeros(length(masses),1);
    E_init = zeros(length(masses),1);
    E_final = zeros(length(masses),1);
    
    for p_index = 1:length(masses)
        v_traj = v_cell{p_index}; % 3 x numsteps
        
        v_init = norm(v_traj(:,1));
        v_end = norm(v_traj(:,end));
        
        gamma_init = 1/sqrt(1 - (v_init/c)^2);
        gamma_end  = 1/sqrt(1 - (v_end/c)^2);
        
        E_init(p_index) = (gamma_init - 1) * masses(p_index) * c^2;
        E_final(p_index) = (gamma_end  - 1) * masses(p_index) * c^2;
        
        deltaE(p_index) = E_final(p_index) - E_init(p_index);
    end
    
    figure;
    histogram(deltaE/1.602e-19, 50, 'Normalization', 'percentage'); % convert J -> eV
    xlabel('\DeltaE (eV)');
    ylabel('Percentage of Particles');
    title('Energy change distribution');
    title(sprintf('Energy Change distribution of Europa Plume %s (θ = %d, Φ = %d)', element, p_theta, p_phi))
end