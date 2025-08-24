function [init_r, init_v_vectors, N_total, particle_traits] = Io_Dust_Starting_States(nums_per_dims, ...
    y_limits, z_limits, starting_x, thetas, direction_angles, init_speeds, q_over_m, radii, densities)

    % Function for setting the initial states of Io dust particles
    % approaching a Jovian moon. nums_per_dims is the number of particles
    % per dimension of the starting grid, y_limits is the maximum and
    % minimum y values of the grid, z_limits is the same but for z, starting_x
    % is the starting x value for the particles, thetas is the angles to be
    % taken by the particles' initial velocities relative to the normal
    % vector pointing in the +x direction, direction_angles is the
    % directions of the angles to be taken if they are at an angle to the
    % norm vector above, init_speeds is the intial speeds of the particles,
    % q_over_m is the q/m ratios of the particles, radii is the radii of
    % the grains, and densities is the densities of the grains.

N_total = nums_per_dims(1) * nums_per_dims(2) * length(thetas) * length(direction_angles) * ...
    length(init_speeds) * length(q_over_m) * length(radii) * length(densities);
init_r = zeros(N_total, 3);
init_r(:,1) = starting_x;
[y, z] = meshgrid(linspace(y_limits(1), y_limits(2), nums_per_dims(1)), ...
    linspace(z_limits(1), z_limits(2), nums_per_dims(2)));

particle_traits = zeros(N_total, nargin - 6); % not counting nums_per_dims, start coords, start velocities
% counting radii + densities as one

init_v_vectors = zeros(N_total, 3);
unit_n = [1, 0, 0];
e_1 = [0, 0, 1];
e_2 = [0, 1, 0];

index = 1;
for a = 1:numel(y)
    for b = 1:length(q_over_m)
        for c = 1:length(radii)
            for d = 1:length(densities)
                for i = 1:length(init_speeds)
                    for j = 1:length(thetas)
                        if thetas(j) ~= 0
                            theta = thetas(j) * pi / 180;
                            for k = 1:length(direction_angles)
                                dir_angle = direction_angles(k) * pi / 180;
    
                                init_r(index,2) = y(a);
                                init_r(index,3) = z(a);
                                particle_traits(index, 2) = 4/3 * pi * radii(c)^3 * densities(d);
                                particle_traits(index, 1) = q_over_m(b) .* particle_traits(index, 2);
                                v = cos(theta)*(unit_n) + sin(theta)*(cos(dir_angle)*e_1 + sin(dir_angle)*e_2);
                                init_v_vectors(index, :) = v .* init_speeds(i);
            
                                index = index + 1;
                            end
                        else
                            theta = 0;
                            dir_angle = 0;

                            init_r(index,2) = y(a);
                            init_r(index,3) = z(a);
                            particle_traits(index, 2) = 4/3 * pi * radii(c)^3 * densities(d);
                            particle_traits(index, 1) = q_over_m(b) .* particle_traits(index, 2);
                            v = cos(theta)*(unit_n) + sin(theta)*(cos(dir_angle)*e_1 + sin(dir_angle)*e_2);
                            init_v_vectors(index, :) = v .* init_speeds(i);
        
                            index = index + 1;
                        end
                    end
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