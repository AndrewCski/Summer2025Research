    % Modelling particle impacts on Europa's surface!
    % To change settings of the simulation, elements/q/m dictate what
    % particles are modelled, init_energies dictates what impact energy the
    % particles have, and the parameters of Spherical_Starting_States will
    % dictate the starting positions and velocities of particles. I do not
    % recommend running Particle_Trajectories when more than a few thousand
    % particles are simulated, as MATLAB will struggle to run it.

% constants

omega = 2*pi/(11.23*3600);                            % synodic frequency (rad/s)
R_E = 1560e3;                                         % Europa radius (m)
c = 299792458;                                        % speed of light (m/s)

% simulation parameters

% elements = ["H+", "O++", "S+++"];
% q = [1.602e-19, 2 * 1.602e-19, 3 * 1.602e-19];      % charge (C)
% m = [1.67262192e-27, 2.656e-26, 5.325e-26];         % mass (kg)

% elements = "H+";
% q = 1.602e-19;
% m = 1.67262192e-27;

% elements = "O++";
% q = 2 * 1.602e-19;
% m = 2.656e-26;

elements = "S+++";
q = 3 * 1.602e-19;
m = 5.325e-26;

timestep = -2 * pi * m / (50 * q * 3.85 * 10^-7);     % 50 steps / gyration, 3.85.. is ave B field
start = 1 * (2 * pi / omega)/4;                       % point in synodic period to start from

numsteps = 10000;                                     % doesnt tend to need above this number
sim_size = 3 * R_E;

B0_vec = [75 * 10 ^-9; 200 * 10 ^-9; 0]; % 75 nT azimuthal, 200 nT radial, using Zimmer coords

init_energies = 1e7;
case_id = "Zimmer";
    
% calling function to set initial positions and velocities of particles
% around surface of Europa, and also get variables for later use

[init_r, init_v_vectors, N_total, particle_traits, elements] = Spherical_Starting_States(...
    -90:5:90, 0:5:359, 0:10:80, 0:60:359, init_energies, m, q, elements, R_E);

particle_traits(:, 5) = mod(particle_traits(:, 5) - 90, 360); % shifting to Nordheim coords

time = linspace(start, start + numsteps * timestep, numsteps);

r_cell = cell(N_total, 1);           % position (m)
v_cell = cell(N_total, 1);           % velocity (m/s)

num_recorded = ones(N_total, 1);     % number of timesteps recorded for each particles
escaped = ones(N_total, 1);          % whether particles trace back to Europa or not
impact_coords = zeros(N_total, 6);   % coords before and after going through surface
done = zeros(N_total, 1);            % whether each particle is done being modelled or not

update_gap = 10000;

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
        % B_prim_now = real(B0_vec .* exp(1i * omega * time(u))); % B_z set in Boris loop
        B_prim_now = [B0_vec(1) * cos(omega * time(u)); B0_vec(2) * -sin(omega * time(u)); 0];

        [r_new, v_new] = new_relativistic_Boris(r_now, v_now, timestep, B_prim_now, ...
            true, false, case_id, particle_traits(p_index, 2), particle_traits(p_index, 3));
        
        r_now = r_new;
        v_now = v_new;

        r_temp(:, u+1) = r_now;
        v_temp(:, u+1) = v_now;
        
        num_recorded(p_index) = u + 1;

        r_norm2 = sum(r_now.^2);  % avoid sqrt for speed
        if r_norm2 < R_E^2
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
    
    if mod(p_index, update_gap) == 0
        fprintf('%d/%d\n', p_index, N_total);
    end
end

%%

Backwards_Particle_Accessibilities(particle_traits(:, 4), particle_traits(:, 5), ...
    escaped, init_energies, elements(1), 5);
% Particle_Trajectories(r_cell, "elements", elements, 0, "Europa");
Particle_Flux_2D(r_cell, 1, numsteps, 'all', 100, sim_size, init_energies, elements(1), ...
    timestep, start, "Europa", update_gap);
Particle_Flux_Iso(r_cell, 1, numsteps, 100, sim_size, init_energies, elements(1), ...
    timestep, start, "Europa", [10, 5], update_gap);