% constants

% q = [1.602e-19, 2 * 1.602e-19, 3 * 1.602e-19];      % charge (C)
% m = [1.67262192e-27, 2.656e-26, 5.325e-26];         % mass (kg)
% elements = ["H+", "O++", "S+++"];
q = 1.602e-19;      % charge (C)
m = 1.67262192e-27;         % mass (kg)
elements = "H+";

omega = 2*pi/(11.23*3600);               % synodic frequency (rad/s)
R_E = 1560e3;                            % Europa radius (m)
c = 299792458;                           % speed of light (m/s)

% simulation parameters
% numsteps = round(abs(period / timestep));   % number of time steps
numsteps = 3000;                              % doesnt tend to need above this number
% period = 2 * pi / omega;                    % modelling a synodic period

% 50 timesteps per gyration, 3.65.. is average primary B field
timestep = -2 * pi * m / (50 * q * 3.659 * 10^-7); 
time = linspace(0, numsteps * timestep, numsteps);
%time = linspace(0, period, numsteps);

B0_vec = [75 * 10 ^-9; 200 * 10 ^-9; 0];
% 75 nT azimuthal, 200 nT radial, using Zimmer coords
B_prim = real(B0_vec .* exp(1i * omega * time));
B_prim(3,:) = -350 * 10^-9;

init_energies = 1e7;

% calling function to set initial positions and velocities of particles
% around surface of Europa, and also get variables for later use

[init_r, init_v_vectors, N_total, particle_traits, elements] = Spherical_Starting_States(...
    90:-5:-90, 0:5:359, 0:20:80, 0:120:359, init_energies, m, q, elements);

r = zeros(3, numsteps, N_total);     % position (m)
v = zeros(3, numsteps, N_total);     % velocity (m/s)

r(:, 1, :) = transpose(init_r);          % final position
v(:, 1, :) = transpose(init_v_vectors);  % final velocity

num_recorded = ones(N_total, 1);     % number of timesteps recorded for each particles
escaped = ones(N_total, 1);          % whether particles trace back to Europa or not
impact_coords = zeros(N_total, 6);   % coords before and after going through surface
done = zeros(N_total, 1);            % whether each particle is done being modelled or not

for u = 1:numsteps-1
    for p_index = 1:N_total
        if ~done(p_index)
            r_now = r(:,u,p_index);
            v_now = v(:,u,p_index);
            [r_new, v_new] = Relativistic_Boris(r_now,v_now,timestep,B_prim(:,u),...
                1, 1, particle_traits(p_index, 2), particle_traits(p_index, 3));
    
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
                done(p_index) = 1;
            elseif r_norm > 8 * R_E
                done(p_index) = 1;
            end
        end
    end
    if all(done == 1)
        break
    end
    if mod(u,500) == 0
        disp(u)
    end
end

%%

Europa_Particle_Accessibilities(particle_traits(:, 4), particle_traits(:, 5), ...
    escaped, init_energies, elements(1), 10);
% Europa_Particle_Trajectories(r, num_recorded, "elements", elements, 0);
% Europa_Particle_Impacts_Backwards(escaped, init_r, particle_traits(:, 1));
% Europa_Particle_Energy_Escapes(escaped, particle_traits(:, 1));