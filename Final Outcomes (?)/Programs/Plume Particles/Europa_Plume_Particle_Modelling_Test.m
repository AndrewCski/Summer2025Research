    % Modelling plume ejecta particle behavior!
    % To change settings, plume_theta and plume_phi dictate the latitude
    % and longitude of the plume's base, respectively, arguments to
    % Plume_Particle_Starting_States will determine the initial states of
    % particles, and the two sections below for dust and ions will dictate
    % which type of particle is simulated (just uncomment the one you want
    % and comment out the other).

% constants

omega = 2*pi/(11.23*3600);                            % synodic frequency (rad/s)
R_E = 1560e3;                                         % Europa radius (m)
c = 299792458;                                        % speed of light (m/s)

% simulation parameters

plume_theta = 45;
plume_phi = 180;

elements = "dust";

ice_velocities = [330, 325, 320, 310, 290, 265, 235, 190, 145, 110, 80, 55, 40, 25, ...
    15, 10, 5, 2, 1, 0.5, 0.1];
ice_radii = [1e-7, 1.5e-7, 2.5e-7, 4e-7, 6.5e-7, 1e-6, 1.5e-6, 2.5e-6, 4e-6, 6.5e-6, 1e-5, ...
    1.5e-5, 2.5e-5, 4e-5, 6.5e-5, 1e-4, 1.5e-4, 2.5e-4, 4e-4, 6.5e-4, 1e-3];
ice_masses = 920 .* 4 .* pi ./ 3 .* ice_radii.^3;  

% [init_r, init_v, N_total, particle_traits] = Plume_Particle_Starting_States(plume_theta, plume_phi, ...
%     15, 1, 15, ice_velocities, [0, 1, 10], ice_masses, "ice");
[init_r, init_v, N_total, particle_traits] = Plume_Particle_Starting_States(plume_theta, plume_phi, ...
    15, 7.5, 60, ice_velocities(1:2:length(ice_velocities)), [0, 10], ...
    ice_masses(1:2:length(ice_masses)), "ice");

timestep = 0.05;
numsteps = 15000;                                     % doesnt tend to need above this number
sim_size = 1.5 * R_E;

% elements = "H+";
% particle_masses = 1.67262192e-27;  
% particle_charges = 1.602e-19;
% 
% particle_speeds = linspace(1e2,1e4,30);
% 
% [init_r, init_v, N_total, particle_traits] = Plume_Particle_Starting_States(plume_theta, plume_phi, ...
%     15, 7.5, 30, particle_speeds, particle_charges./particle_masses, particle_masses, "ions");
% % [init_r, init_v, N_total, particle_traits] = Plume_Particle_Starting_States(plume_theta, ...
% % plume_phi, 15, 15, 45, particle_speeds, particle_charges./particle_masses, particle_masses, "ions");
% 
% timestep = 2 * pi * particle_masses / (100 * particle_charges * 3.85 * 10^-7);
% numsteps = 25000;                                     % doesnt tend to need above this number
% sim_size = 2.5 * R_E;

start = 0 * (2 * pi / omega)/4;                       % point in synodic period to start from

B0_vec = [75 * 10 ^-9; 200 * 10 ^-9; 0]; % 75 nT azimuthal, 200 nT radial, using EPhiO

case_id = "Zimmer";

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
    v_now = init_v(p_index, :)';
    
    r_temp(:,1) = r_now;
    v_temp(:,1) = v_now;

    for u = 1:numsteps-1
        B_prim_now = [B0_vec(1) * cos(omega * time(u)); B0_vec(2) * -sin(omega * time(u)); 0];
        
        [r_new, v_new] = new_relativistic_Boris(r_now, v_now, timestep, B_prim_now, ...
            true, true, case_id, particle_traits(p_index, 5), particle_traits(p_index, 4));
        
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
            break;
        elseif r_norm2 > (sim_size)^2
            break;
        end
    end
    
    % save truncated histories
    r_cell{p_index} = r_temp(:, 1:num_recorded(p_index));
    v_cell{p_index} = v_temp(:, 1:num_recorded(p_index));
    
    if mod(p_index, 5000) == 0
        fprintf('%d/%d\n', p_index, N_total);
    end
end
clear('r_temp');
clear('v_temp');

%%

Particle_Trajectories(r_cell, "idk", elements, 0, "Europa");

Plume_Particle_Energy_Shifts(v_cell, particle_traits(:, 5), elements(1), plume_theta, plume_phi);

Plume_Particle_Height_Dist(r_cell, particle_traits(:, 5), elements(1), plume_theta, plume_phi);

% Plume_Particle_Max_Heights(r_cell, particle_masses, particle_speeds, particle_traits(:, 5), ...
%     particle_traits(:, 3), "ions", plume_theta, plume_phi);

Plume_Particle_Max_Heights(r_cell, ice_masses, ice_velocities, particle_traits(:, 5), ...
    particle_traits(:, 3), "ice", plume_theta, plume_phi)