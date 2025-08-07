% constants

omega = 2*pi/(10.53*3600);                            % synodic frequency (rad/s)
R_G = 2631e3;                                         % Ganymede radius (m)
c = 299792458;                                        % speed of light (m/s)

% simulation parameters

elements = "H+";
q = 1.602e-19;
m = 1.67262192e-27;

B0_vec = [-15 * 10 ^-9; 24 * 10 ^-9; 0]; % 75 nT azimuthal, 200 nT radial, using Zimmer coords

timestep = -2 * pi * m / (50 * q * norm(B0_vec));     % 50 steps / gyration
start = 0 * (2 * pi / omega)/4;                       % point in synodic period to start from
% start = 0;                                           % point in synodic period to start from

% numsteps = round(abs((2 * pi / omega) / timestep)); % steps for a synodic period
numsteps = 10000;                                     % doesnt tend to need above this number
sim_size = 5 * R_G;

init_energies = 1e5;
case_id = "Ganymede";
    
% calling function to set initial positions and velocities of particles
% around surface of Europa, and also get variables for later use

[init_r, init_v_vectors, N_total, particle_traits, elements] = Spherical_Starting_States(...
    -90:30:90, 0:30:359, 0:20:80, 0:60:359, init_energies, m, q, elements, R_G);

particle_traits(:, 5) = mod(particle_traits(:, 5) - 90, 360); % shifting to Nordheim coords

time = linspace(start, start + numsteps * timestep, numsteps);

r_cell = cell(N_total, 1);           % position (m)
v_cell = cell(N_total, 1);           % velocity (m/s)

num_recorded = ones(N_total, 1);     % number of timesteps recorded for each particles
escaped = ones(N_total, 1);          % whether particles trace back to Europa or not
impact_coords = zeros(N_total, 6);   % coords before and after going through surface
done = zeros(N_total, 1);            % whether each particle is done being modelled or not

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
        B_prim_now = real(B0_vec .* exp(0 * omega * time(u))); % B_z set in Boris loop
        
        [r_new, v_new] = Relativistic_Boris(r_now, v_now, timestep, B_prim_now, ...
            0, case_id, particle_traits(p_index, 2), particle_traits(p_index, 3));
        
        r_now = r_now + v_new * timestep;
        v_now = v_new;

        r_temp(:, u+1) = r_now;
        v_temp(:, u+1) = v_now;
        
        num_recorded(p_index) = u + 1;

        r_norm2 = sum(r_now.^2);  % avoid sqrt for speed
        if r_norm2 < R_G^2
            escaped(p_index) = 0;
            impact_coords(p_index, 1:3) = r_temp(:, u);
            impact_coords(p_index, 4:6) = r_temp(:, u+1);
            break
        elseif r_norm2 > (sim_size)^2
            break
        end
    end
    
    % save truncated histories
    r_cell{p_index} = r_temp(:, 1:num_recorded(p_index));
    v_cell{p_index} = v_temp(:, 1:num_recorded(p_index));
    
    if mod(p_index, 1000) == 0
        fprintf('%d/%d\n', p_index, N_total);
    end
end

%%

Particle_Accessibilities(particle_traits(:, 4), particle_traits(:, 5), ...
    escaped, init_energies, elements(1), 30);
Particle_Trajectories(r_cell, "elements", elements, 0, "Ganymede");
% Europa_Particle_Impacts_Backwards(escaped, init_r, particle_traits(:, 1));
% Europa_Particle_Energy_Escapes(escaped, particle_traits(:, 1));
% Europa_Surface_Field_Map(real((B0_vec) .* exp(1i * omega * time)),
% case_id, 1, timestep, start);
Particle_Flux(r_cell, 1, numsteps, 'all', 100, sim_size, init_energies, elements(1), ...
    timestep, start, "Ganymede");