function [init_r, init_v_vectors, N_total, particle_traits, elements] = Spherical_Starting_States(...
    theta_list, phi_list, alpha_list, beta_list, init_energies, masses, charges, elements, body_radius)

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

c = 299792458;                                                    % speed of light (m/s)
E_to_v = @(E, m) c .* sqrt(1 - 1 ./ ((E ./ (m .* c.^2) + 1).^2)); % get velocity from mass/energy
E_to_v_eV = @(E, m) E_to_v(E .* 1.60218 .* 10 .^ -19, m); 

num_points = length(theta_list) * length(phi_list);
num_directions = length(alpha_list) * length(beta_list);
N_total = num_points * num_directions * length(init_energies) * length(masses);

init_r = zeros(N_total, 3);
init_v_vectors = zeros(N_total, 3);

particle_traits = ones(N_total, 5);
element_array = strings(N_total, 1);

index = 1;
for g = 1:length(init_energies)
    for h = 1:length(masses)
        for i = 1:length(theta_list)
            theta = deg2rad(theta_list(i));
            if abs(theta) ~= pi/2                                 % if at poles, phi doesnt matter
                for j = 1:length(phi_list)
                    phi = deg2rad(phi_list(j));
            
                    % sphere position
                    x = body_radius * cos(theta) * sin(phi);
                    y = body_radius * cos(theta) * cos(phi);
                    z = body_radius * sin(theta);
                    point = [x, y, z];
            
                    % normal vector (outward)
                    n = [x, y, z] / body_radius;
            
                    % local tangent basis
                    e_theta = [-sin(theta)*cos(phi), -sin(theta)*sin(phi), cos(theta)];
                    e_phi   = [-sin(phi), cos(phi), 0];
            
                    % loop over incidence angles and directions
                    for k = 1:length(alpha_list)
                        alpha = deg2rad(alpha_list(k));
                        if abs(alpha) ~= 0    % if radially inward, beta doesnt matter
                            for a = 1:length(beta_list)
                                beta = deg2rad(beta_list(a));
                
                                % tilted vector (inward direction)
                                v = cos(alpha)*(-n) + sin(alpha)*(cos(beta)*e_theta + sin(beta)*e_phi);
                
                                % storing positions and velocities
                                init_r(index, :) = point;
                                init_v_vectors(index, :) = v .* E_to_v_eV(init_energies(g), masses(h));
        
                                particle_traits(index, 1) = init_energies(g);
                                particle_traits(index, 2) = masses(h);
                                particle_traits(index, 3) = charges(h);
                                particle_traits(index, 4) = theta_list(i);
                                particle_traits(index, 5) = phi_list(j);
                                element_array(index) = elements(h);
        
                                index = index + 1;
                            end
                        else
                            beta = 0;
            
                            % tilted vector (inward direction)
                            v = cos(alpha)*(-n) + sin(alpha)*(cos(beta)*e_theta + sin(beta)*e_phi);
            
                            % storing positions and velocities
                            init_r(index, :) = point;
                            init_v_vectors(index, :) = v .* E_to_v_eV(init_energies(g), masses(h));
    
                            particle_traits(index, 1) = init_energies(g);
                            particle_traits(index, 2) = masses(h);
                            particle_traits(index, 3) = charges(h);
                            particle_traits(index, 4) = theta_list(i);
                            particle_traits(index, 5) = phi_list(j);
                            element_array(index) = elements(h);
    
                            index = index + 1;
                        end
                    end
                end
            else
                phi = 0;
        
                % sphere position
                x = body_radius * cos(theta) * cos(phi);
                y = body_radius * cos(theta) * sin(phi);
                z = body_radius * sin(theta);
                point = [x, y, z];
        
                % normal vector (outward)
                n = [x, y, z] / body_radius;
        
                % local tangent basis
                e_theta = [-sin(theta)*cos(phi), -sin(theta)*sin(phi), cos(theta)];
                e_phi   = [-sin(phi), cos(phi), 0];
        
                % loop over incidence angles and directions
                for k = 1:length(alpha_list)
                    alpha = deg2rad(alpha_list(k));
                    if abs(alpha) ~= 0
                        for a = 1:length(beta_list)
                            beta = deg2rad(beta_list(a));
            
                            % tilted vector (inward direction)
                            v = cos(alpha)*(-n) + sin(alpha)*(cos(beta)*e_theta + sin(beta)*e_phi);
            
                            % storing positions and velocities
                            init_r(index, :) = point;
                            init_v_vectors(index, :) = v .* E_to_v_eV(init_energies(g), masses(h));
    
                            particle_traits(index, 1) = init_energies(g);
                            particle_traits(index, 2) = masses(h);
                            particle_traits(index, 3) = charges(h);
                            particle_traits(index, 4) = theta_list(i);
                            particle_traits(index, 5) = 0;
                            element_array(index) = elements(h);
    
                            index = index + 1;
                        end
                    else
                        beta = 0;
        
                        % tilted vector (inward direction)
                        v = cos(alpha)*(-n) + sin(alpha)*(cos(beta)*e_theta + sin(beta)*e_phi);
        
                        % storing positions and velocities
                        init_r(index, :) = point;
                        init_v_vectors(index, :) = v .* E_to_v_eV(init_energies(g), masses(h));

                        particle_traits(index, 1) = init_energies(g);
                        particle_traits(index, 2) = masses(h);
                        particle_traits(index, 3) = charges(h);
                        particle_traits(index, 4) = theta_list(i);
                        particle_traits(index, 5) = 0;
                        element_array(index) = elements(h);

                        index = index + 1;
                    end
                end
            end
        end
    end
end

% need to get rid of rows that are not used

rows_to_delete = all(init_r == 0, 2);

init_r = init_r(~rows_to_delete, :);
init_v_vectors = init_v_vectors(~rows_to_delete, :);
particle_traits = particle_traits(~rows_to_delete, :);
element_array = element_array(~rows_to_delete, :);
elements = element_array;

N_total = length(elements);
