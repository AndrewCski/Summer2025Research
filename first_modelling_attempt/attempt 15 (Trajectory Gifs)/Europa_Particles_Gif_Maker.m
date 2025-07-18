% constants

q = [1.602e-19, 2 * 1.602e-19, 3 * 1.602e-19];      % charge (C)
m = [1.67262192e-27, 2.656e-26, 5.325e-26];         % mass (kg)
elements = ["H+", "O++", "S+++"];

omega = 2*pi/(11.23*3600);               % synodic frequency (rad/s)
R_E = 1560e3;                            % Europa radius (m)
c = 299792458;                           % speed of light (m/s)

% simulation parameters
numsteps = 6000;                              % doesnt tend to need above this number

% 50 timesteps per gyration, 3.65.. is average primary B field
timestep = -2 * pi * m / (100 * q * 3.659 * 10^-7); 
time = linspace(0, numsteps * timestep, numsteps);

B0_vec = [75 * 10 ^-9; 200 * 10 ^-9; 0];
% 75 nT azimuthal, 200 nT radial, using Zimmer coords
B_prim = real(B0_vec .* exp(1i * omega * time));
B_prim(3,:) = -350 * 10^-9;

init_energies = 1000000;

% calling function to set initial positions and velocities of particles
% around surface of Europa, and also get variables for later use
[init_r, init_v_vectors, N_total, particle_traits, elements] = Spherical_Starting_States(90:-90:-90, ...
    0:120:359, 0, 0, init_energies, m, q, elements);

r = zeros(3, numsteps, N_total);     % position (m)
v = zeros(3, numsteps, N_total);     % velocity (m/s)

r(:, 1, :) = transpose(init_r);          % final position
v(:, 1, :) = transpose(init_v_vectors);  % final velocity

num_recorded = ones(N_total, 1);     % number of timesteps recorded for each particles
done = zeros(N_total, 1);          % whether each particle is done being modelled or not

filename = 'Particle_Position_Evolution_Test.gif';

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
            if r_norm < R_E || r_norm > 10 * R_E
                done(p_index) = 1;
            end
        end
    end
    if mod(u,5) == 0
        Europa_Particle_Trajectories(r, num_recorded, "elements", elements); % plot results
    
        % capture the frame
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
    
        % write to GIF
        if u == 20  % first frame
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.3);
        else
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.3);
        end
        close(gcf);  % close the figure to avoid many open windows    
    end
    if all(done == 1)
        break
    end
end