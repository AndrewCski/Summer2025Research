function [init_r, init_v_vectors, N_total, particle_traits] = Io_Dust_Starting_States(...
    nums_per_dims, y_limits, z_limits, starting_x, init_speeds, q_over_m, radii, densities)

N_total = nums_per_dims(1) * nums_per_dims(2) * length(q_over_m) * length(radii) * ...
    length(densities) * length(init_speeds);
init_r = zeros(N_total, 3);
init_r(:,1) = starting_x;
[y, z] = meshgrid(linspace(y_limits(1), y_limits(2), nums_per_dims(1)), ...
    linspace(z_limits(1), z_limits(2), nums_per_dims(2)));

particle_traits = zeros(N_total, nargin - 6); % not counting nums_per_dims, start coords, or start velocities.
% counting radii + densities as one

init_v_vectors = zeros(N_total, 3);
init_v_vectors(:,1) = init_speeds;

index = 1;
for i = 1:numel(y)
    for j = 1:length(q_over_m)
        for k = 1:length(radii)
            for l = 1:length(densities)
                for p = 1:length(init_speeds)
                    init_r(index,2) = y(i);
                    init_r(index,3) = z(i);
                    init_v_vectors(index,1) = init_speeds(p);
                    particle_traits(index, 2) = 4/3 * pi * radii(k)^3 * densities(l);
                    particle_traits(index, 1) = q_over_m(j) .* particle_traits(index, 2);

                    index = index + 1;
                end
            end
        end
    end
end

