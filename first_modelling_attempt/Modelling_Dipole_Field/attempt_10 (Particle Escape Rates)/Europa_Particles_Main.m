% constants

q = 1.602e-19;              % charge of proton (C)
m = 1.67262192e-27;         % mass of proton (kg)
omega = 2*pi/(11.23*3600);  % synodic frequency (rad/s)
R_E = 1560e3;               % Europa radius (m)
c = 299792458;              % speed of light (m/s)

% simulation parameters
% numsteps = round(abs(period / timestep));   % number of time steps
numsteps = 15000;                             % doesnt tend to need above this number
% period = 2 * pi / omega;                    % modelling a synodic period

% 50 timesteps per gyration, 3.65.. is average primary B field
timestep = -2 * pi * m / (50 * q * 3.659 * 10^-7); 
time = linspace(0, numsteps * timestep, numsteps);
%time = linspace(0, period, numsteps);

B0_vec = [75 * 10 ^-9; 200 * 10 ^-9; 0];
% 75 nT azimuthal, 200 nT radial, using Zimmer coords
% (-350 nT constant added to North-South (z) in upcoming loop)
B_prim = real(B0_vec .* exp(1i * omega * time));
B_prim(3,:) = -350 * 10^-9;

% functions to create velocities for particles of desired energies 
E_to_v = @(E, m) sqrt((E.^2./(m.^2 .* c.^2) + 2 .* E ./ m)./...
    (1 + (E.^2./(m.^2 .* c.^2) + 2 .* E ./ m)./c.^2));
E_to_v_eV = @(E, m) E_to_v(E .* 1.60218 .* 10 .^ -19, m); % convert eV to J
init_energies = [100, 1000, 10000, 100000];
initial_velocities = E_to_v_eV(init_energies, m);

% calling function to set initial positions and velocities of particles
% around surface of Europa, and also get variables for later use
[init_r, init_v_vectors, parts_per_E, N_total] = Spherical_Starting_States(90:-30:-90, ...
    0:60:359, 0, 0, initial_velocities);

r = zeros(3, numsteps, N_total);     % position (m)
v = zeros(3, numsteps, N_total);     % velocity (m/s)

r(:, 1, :) = transpose(init_r);          % final position
v(:, 1, :) = transpose(init_v_vectors);  % final velocity

num_recorded = ones(N_total, 1);     % number of timesteps recorded for each particles
escaped = ones(N_total, 1);          % whether particles trace back to Europa or not
impact_coords = zeros(N_total, 6);   % coords before and after going through surface

for p_index = 1:N_total
    for u = 1:numsteps-1
        r_now = r(:,u,p_index);
        v_now = v(:,u,p_index);
        [r_new, v_new] = Relativistic_Boris(r_now,v_now,timestep,B_prim(:,u),[], m, q);

        % update position and velocity
        r(:,u+1,p_index) = r_now + v_new * timestep;
        v(:,u+1,p_index) = v_new;

        % track number of steps
        num_recorded(p_index) = num_recorded(p_index) + 1;

        % stop if particle hits Europa or escapes
        r_norm = norm(r(:,u+1,p_index));
        if r_norm < R_E
            escaped(p_index) = 0;
            impact_coords(p_index, 1:3) = r(:,u,p_index);
            impact_coords(p_index, 4:6) = r(:,u+1,p_index);
            break
        elseif r_norm > 10 * R_E
            break
        end
    end
end

energies = zeros(N_total, 1);
for i = 1:N_total
    % figure out which initial velocity this particle has:
    h = ceil(i / (parts_per_E));
    energies(i) = init_energies(h); % in eV
end

Europa_Particle_Trajectories(r, num_recorded, energies);
Europa_Particle_Impacts(impact_coords, energies);
Europa_Particle_Energy_Escapes(escaped, energies);