function [init_r, init_v_vectors, N_total, particle_traits] = Plume_Particle_Starting_States(plume_theta, ...
    plume_phi, theta_p, d_theta, direction_angle_gap, init_speeds, q_over_m, masses, mode)

    % Function for creating the starting states of plume particles.
    % plume_theta should be the latitude of the base of the plume on
    % Europa's surface, plume_phi should be its longitude, theta_p is the
    % plume's opening angle, d_theta is the gap between each starting
    % normal angle for the particles started between 0 and p_theta, 
    % direction_angle_gap is the gap between each starting direction angle
    % for the particles started between 0 and 359 degrees, init_speeds are
    % the initial velocity magnitude of the particles, q_over_m is the q/m
    % ratio of the particles, masses is the masses of the particles, and
    % mode is for whether ions or dust is being simulated (dust = mass/speed
    % coupled, ions = mass/speed unrelated)

R_E = 1560e3;                                         % Europa radius (m)

plume_phi = mod(90 - plume_phi, 360) * pi / 180;
plume_theta = plume_theta * pi / 180;
plume_r = [(R_E) * cos(plume_theta) * cos(plume_phi), ...
    (R_E) * cos(plume_theta) * sin(plume_phi), (R_E) * sin(plume_theta)];
e_theta = [-sin(plume_theta)*cos(plume_phi), -sin(plume_theta)*sin(plume_phi), cos(plume_theta)];
e_phi   = [-sin(plume_phi), cos(plume_phi), 0];

direction_angles = 0:direction_angle_gap:359;

% normal vector (outward)
n = plume_r / R_E;

thetas = 0:d_theta:theta_p;

N_total = length(thetas) * length(direction_angles) * length(init_speeds) * length(q_over_m) * ...
    length(masses);
init_r = zeros(N_total, 3);
init_r(:,1) = plume_r(1);
init_r(:,2) = plume_r(2);
init_r(:,3) = plume_r(3);

particle_traits = zeros(N_total, nargin - 4); % not counting plume coords or opening angle
init_v_vectors = zeros(N_total, 3);

index = 1;
if lower(mode) == "ion" || lower(mode) == "ions" % masses and init_speeds not coupled 
    for a = 1:length(thetas)
        theta = deg2rad(thetas(a));
        if abs(theta) ~= 0
            for b = 1:length(direction_angles)
                dir_angle = deg2rad(direction_angles(b));
                for c = 1:length(init_speeds)
                    for d = 1:length(q_over_m)
                        for i = 1:length(masses)
                            particle_traits(index, 1) = thetas(a);
                            particle_traits(index, 2) = direction_angles(b);
                            particle_traits(index, 3) = init_speeds(c);
                            particle_traits(index, 5) = masses(i);
                            particle_traits(index, 4) = q_over_m(d) .* particle_traits(index, 5);
                                    
                            % tilted vector (outward)
                            v = cos(theta)*(n) + sin(theta)*(cos(dir_angle)*e_theta + ...
                                sin(dir_angle)*e_phi);
                            init_v_vectors(index, :) = v .* init_speeds(c);
    
                            index = index + 1;
                        end
                    end
                end
            end
        else
            dir_angle = 0;
            for c = 1:length(init_speeds)
                for d = 1:length(q_over_m)
                    for i = 1:length(masses)
                        particle_traits(index, 1) = thetas(a);
                        particle_traits(index, 2) = 0;
                        particle_traits(index, 3) = init_speeds(c);
                        particle_traits(index, 5) = masses(i);
                        particle_traits(index, 4) = q_over_m(d) .* particle_traits(index, 5);
                                
                        % tilted vector (outward)
                        v = cos(theta)*(n) + sin(theta)*(cos(dir_angle)*e_theta + sin(dir_angle)*e_phi);
                        init_v_vectors(index, :) = v .* init_speeds(c);

                        index = index + 1;
                    end
                end
            end
        end
    end
elseif lower(mode) == "dust" || lower(mode) == "ice" % masses and init_speeds coupled
    for a = 1:length(thetas)
        theta = deg2rad(thetas(a));
        if abs(theta) ~= 0
            for b = 1:length(direction_angles)
                dir_angle = deg2rad(direction_angles(b));
                for c = 1:length(init_speeds)
                    for d = 1:length(q_over_m)
                        particle_traits(index, 1) = thetas(a);
                        particle_traits(index, 2) = direction_angles(b);
                        particle_traits(index, 3) = init_speeds(c);
                        particle_traits(index, 5) = masses(c);
                        particle_traits(index, 4) = q_over_m(d) .* particle_traits(index, 5);
                                
                        % tilted vector (outward)
                        v = cos(theta)*(n) + sin(theta)*(cos(dir_angle)*e_theta + sin(dir_angle)*e_phi);
                        init_v_vectors(index, :) = v .* init_speeds(c);

                        index = index + 1;
                    end
                end
            end
        else
            dir_angle = 0;
            for c = 1:length(init_speeds)
                for d = 1:length(q_over_m)
                    particle_traits(index, 1) = thetas(a);
                    particle_traits(index, 2) = 0;
                    particle_traits(index, 3) = init_speeds(c);
                    particle_traits(index, 5) = masses(c);
                    particle_traits(index, 4) = q_over_m(d) .* particle_traits(index, 5);
                            
                    % tilted vector (outward)
                    v = cos(theta)*(n) + sin(theta)*(cos(dir_angle)*e_theta + sin(dir_angle)*e_phi);
                    init_v_vectors(index, :) = v .* init_speeds(c);

                    index = index + 1;
                end
            end
        end
    end
end

rows_to_delete = all(particle_traits == 0, 2);

init_r = init_r(~rows_to_delete, :);
init_v_vectors = init_v_vectors(~rows_to_delete, :);
particle_traits = particle_traits(~rows_to_delete, :);
N_total = size(init_r, 1);