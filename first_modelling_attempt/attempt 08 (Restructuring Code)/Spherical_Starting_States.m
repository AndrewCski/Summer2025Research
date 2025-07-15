function [init_r, init_v_vectors, parts_per_E, N_total] = Spherical_Starting_States(...
    theta_list, phi_list, alpha_list, beta_list, initial_velocities)

% Function for setting the initial velocities and positions of particles
% set to start on the surface of Europa. 
% 
% Takes lists of latitudinal and longitudinal positions, in addition to
% directions and pitch angles of impacts. Initial velocity magnitudes are
% also needed.
%
% Outputs computed initial positions and velocities of all particles, in
% addition to number of particles simulated for every starting energy and
% number of particles simulated overall.

num_points = length(theta_list) * length(phi_list);
num_directions = length(alpha_list) * length(beta_list);
parts_per_E = num_points * num_directions;
N_total = num_points * num_directions * length(initial_velocities);

R_E = 1560e3;               % Europa radius in meters
init_r = zeros(N_total, 3);
init_v_vectors = zeros(N_total, 3);

index = 1;
for h = 1: length(initial_velocities)
    for i = 1:length(theta_list)
        theta = deg2rad(theta_list(i));
    
        for j = 1:length(phi_list)
            phi = deg2rad(phi_list(j));
    
            % sphere position
            x = R_E * cos(theta) * cos(phi);
            y = R_E * cos(theta) * sin(phi);
            z = R_E * sin(theta);
            point = [x, y, z];
    
            % normal vector (outward)
            n = [x, y, z] / R_E;
    
            % local tangent basis
            e_theta = [-sin(theta)*cos(phi), -sin(theta)*sin(phi), cos(theta)];
            e_phi   = [-sin(phi), cos(phi), 0];
    
            % loop over incidence angles and directions
            for k = 1:length(alpha_list)
                alpha = deg2rad(alpha_list(k));
    
                for a = 1:length(beta_list)
                    beta = deg2rad(beta_list(a));
    
                    % tilted vector (inward direction)
                    v = cos(alpha)*(-n) + sin(alpha)*(cos(beta)*e_theta + sin(beta)*e_phi);
    
                    % storing positions and velocities
                    init_r(index, :) = point;
                    init_v_vectors(index, :) = v .* initial_velocities(h);
    
                    index = index + 1;
                end
            end
        end
    end
end
