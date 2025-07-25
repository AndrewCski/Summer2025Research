% constants

% q = [1.602e-19, 2 * 1.602e-19, 3 * 1.602e-19];      % charge (C)
% m = [1.67262192e-27, 2.656e-26, 5.325e-26];         % mass (kg)
% elements = ["H+", "O++", "S+++"];
q = 1.602e-19;              % charge (C)
m = 1.67262192e-27;         % mass (kg)
elements = "H+";

omega = 2*pi/(11.23*3600);               % synodic frequency (rad/s)
R_E = 1560e3;                            % Europa radius (m)
c = 299792458;                           % speed of light (m/s)

% simulation parameters
% numsteps = round(abs(period / timestep));   % number of time steps
numsteps = 3000;                             % doesnt tend to need above this number
max_steps_per_particle = numsteps;
period = 2 * pi / omega;                      % modelling a synodic period

% 50 timesteps per gyration, 3.65.. is average primary B field
timestep = -2 * pi * m / (50 * q * 3.659 * 10^-7); 
start = 0;
time = linspace(start, start + numsteps * timestep, numsteps);
%time = linspace(0, period, numsteps);

B0_vec = [75 * 10 ^-9; 200 * 10 ^-9; 0]; % 75 nT azimuthal, 200 nT radial, using Zimmer coords
% B0_vec = [0; 0; 0];

init_energies = 1e3;
casenum = 0;
    
% calling function to set initial positions and velocities of particles
% around surface of Europa, and also get variables for later use

[init_r, init_v_vectors, N_total, particle_traits, elements] = Spherical_Starting_States(...
    -90:45:90, 0:45:359, 0:20:80, 0:90:359, init_energies, m, q, elements);
% B_prim = real(B0_vec .* exp(1i * omega * time));
% [init_r, init_v_vectors, N_total, particle_traits, elements] = New_Spherical_Starting_States(...
%     B_prim(:,1), casenum, 90:-90:-90, 0:90:359, [0:40:80 100:40:180], 0, init_energies, ...
%     m, q, elements);

particle_traits(:, 5) = mod(particle_traits(:, 5) - 90, 360); % shifting to nordheim coords

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
        B_prim_now = real(B0_vec .* exp(1i * omega * time(u))); % z comp counted in Boris loop
        
        [r_new, v_new] = Relativistic_Boris(r_now, v_now, timestep, B_now, ...
            0, 0, casenum, particle_traits(p_index, 2), particle_traits(p_index, 3));
        
        r_now = r_now + v_new * timestep;
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
        elseif r_norm2 > (3*R_E)^2
            break
        end
    end
    
    % save truncated histories
    r_cell{p_index} = r_temp(:, 1:num_recorded(p_index));
    v_cell{p_index} = v_temp(:, 1:num_recorded(p_index));
    
    if mod(p_index, 200) == 0
        disp(p_index)
    end
end

%%

Europa_Particle_Accessibilities(particle_traits(:, 4), particle_traits(:, 5), ...
    escaped, init_energies, elements(1), 45);
Europa_Particle_Trajectories_Opt(r_cell, "elements", elements, 0);
% Europa_Particle_Impacts_Backwards(escaped, init_r, particle_traits(:, 1));
% Europa_Particle_Energy_Escapes(escaped, particle_traits(:, 1));