% constants

omega = 2*pi/(11.23*3600);                            % synodic frequency (rad/s)
R_E = 1560e3;                                         % Europa radius (m)
c = 299792458;                                        % speed of light (m/s)

% simulation parameters

elements = "dust";

[init_r, init_v_vectors, N_total, particle_traits] = Io_Dust_Starting_States(...
    [5,5], [-1.25*R_E,1.25*R_E], [-1.25*R_E,1.25*R_E], -3*R_E, 100000, 5000, 10e-9, 2500);

% timestep = 2 * pi * m / (100 * q * 3.85 * 10^-7);     % 500 steps / gyration, 3.85.. ave B field
timestep = 0.01;
start = 0 * (2 * pi / omega)/4;                       % point in synodic period to start from

% numsteps = round(abs((2 * pi / omega) / timestep)); % steps for a synodic period
numsteps = 10000;                                     % doesnt tend to need above this number
sim_size = 4 * R_E;

B0_vec = [75 * 10 ^-9; 200 * 10 ^-9; 0]; % 75 nT azimuthal, 200 nT radial, using Zimmer coords

case_id = "Zimmer";
% case_id = "yeah";

time = linspace(start, start + numsteps * timestep, numsteps);

r_cell = cell(N_total, 1);           % position (m)
v_cell = cell(N_total, 1);           % velocity (m/s)

num_recorded = ones(N_total, 1);     % number of timesteps recorded for each particles
escaped = ones(N_total, 1);          % whether particles trace back to Europa or not
impact_coords = zeros(N_total, 6);   % coords before and after going through surface

%%

for p_index = 1:N_total
    % initialize per-particle arrays (overallocate to max possible)
    r_temp = zeros(3, numsteps);
    v_temp = zeros(3, numsteps);
    
    % initial state
    r_now = init_r(p_index, :)';
    v_now = init_v_vectors(p_index, :)';
    
    r_temp(:,1) = r_now;
    v_temp(:,1) = v_now;

    for u = 1:numsteps-1
        B_prim_now = real(B0_vec .* exp(1i * omega * time(u))); % B_z set in Boris loop
        
        [r_new, v_new] = new_relativistic_Boris(r_now, v_now, timestep, B_prim_now, ...
            0, case_id, m, q);
        
        r_now = r_new;
        v_now = v_new;

        r_temp(:, u+1) = r_now;
        v_temp(:, u+1) = v_now;
        
        num_recorded(p_index) = u + 1;

        r_norm2 = sum(r_now.^2);  % avoid sqrt for speed
        if r_norm2 <= R_E^2
            escaped(p_index) = 0;
            impact_coords(p_index, 1:3) = r_temp(:, u);
            impact_coords(p_index, 4:6) = r_temp(:, u+1);
            break;
        elseif r_norm2 > (sim_size)^2
        % elseif r_norm2 > (sim_size)^2 || r_now(1) > 1.5 * R_E % particles are probably not coming back...
            break;
        end
    end
    
    % save truncated histories
    r_cell{p_index} = r_temp(:, 1:num_recorded(p_index));
    v_cell{p_index} = v_temp(:, 1:num_recorded(p_index));
    
    if mod(p_index, 500) == 0
        fprintf('%d/%d\n', p_index, N_total);
    end
end

%%

Particle_Trajectories(r_cell, "idk", elements, 0, "Europa");
forward_particle_start_accessibilities(init_r, escaped, 0, q_over_m, "Europa")
% forward_particle_start_accessibilities(init_r, escaped, [2,2], q_over_m, "Europa")
forward_particle_impact_accessibilities(impact_coords, 5, q_over_m, "Europa")