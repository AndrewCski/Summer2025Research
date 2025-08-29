    % Modelling Io dust passing by Ganymede!
    % To change settings of the simulation, the parameters of Io_Dust_Starting_States 
    % will dictate the starting positions and velocities of particles, and
    % start will choose the point in the synodic period to simulate in. 

% constants

omega = 2*pi/(11.23*3600);                            % synodic frequency (rad/s)
R_G = 2631e3;                                         % Ganymede radius (m)
c = 299792458;                                        % speed of light (m/s)

% simulation parameters

elements = "dust";

[init_r, init_v_vectors, N_total, particle_traits] = Io_Dust_Starting_States([30,30], ...
    [-1.25*R_G,1.25*R_G], [-1.25*R_G,1.25*R_G], -3*R_G, 0, 0, 200000, 1000, 10e-9, 2500);

timestep = 0.01;
start = 1 * (2 * pi / omega)/4;                       % point in synodic period to start from

numsteps = 20000;                                     % doesnt tend to need above this number
sim_size = 5 * R_G;

B0_vec = [15 * 10 ^-9; 80 * 10 ^-9; 0];  

case_id = "Ganymede";

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
        B_prim_now = [B0_vec(1) * cos(omega * time(u)); B0_vec(2) * -sin(omega * time(u) - pi/4.5); 0];
        
        [r_new, v_new] = new_relativistic_Boris(r_now, v_now, timestep, B_prim_now, ...
            true, false, case_id, particle_traits(p_index, 2), particle_traits(p_index, 1));
        
        r_now = r_new;
        v_now = v_new;

        r_temp(:, u+1) = r_now;
        v_temp(:, u+1) = v_now;
        
        num_recorded(p_index) = u + 1;

        r_norm2 = sum(r_now.^2);  % avoid sqrt for speed
        if r_norm2 <= R_G^2
            escaped(p_index) = 0;
            impact_coords(p_index, 1:3) = r_temp(:, u);
            impact_coords(p_index, 4:6) = r_temp(:, u+1);
            break;
        elseif r_norm2 > (sim_size)^2 || r_now(1) > 1.5 * R_G % particles are probably not coming back...
            break;
        end
    end
    
    % save truncated histories
    r_cell{p_index} = r_temp(:, 1:num_recorded(p_index));
    v_cell{p_index} = v_temp(:, 1:num_recorded(p_index));
    
    if mod(p_index, 100) == 0
        fprintf('%d/%d\n', p_index, N_total);
    end
end

%%

Particle_Trajectories(r_cell, "idk", elements, 0, "Ganymede");
forward_particle_start_accessibilities(init_r, escaped, 0, particle_traits(p_index, 1)./ ...
    particle_traits(p_index, 2), "Ganymede")
forward_particle_impact_accessibilities(impact_coords, 5, particle_traits(p_index, 1)./ ...
    particle_traits(p_index, 2), "Ganymede", elements, "q/m", true, false)
