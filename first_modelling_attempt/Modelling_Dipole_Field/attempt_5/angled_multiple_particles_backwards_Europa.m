% constants

q = 1.602e-19;              % charge of proton (C)
m = 1.67262192e-27;         % mass of proton (kg)

omega = 2*pi/(11.23*3600);  % Synodic frequency in rad/s

% simulation parameters
numsteps = 1000000;                      % Number of time steps
period = 2 * pi / omega;        % modelling 1/100 of a synodic period
time = linspace(period/5, period/5.5, numsteps);
timestep = -period / numsteps / 20;

B0_vec = [75 * 10 ^-9; 200 * 10 ^-9; 0];  % 75 nT azimuthal, 200 nT radial, using Zimmer coords
% (-350 nT constant added to North-South, represented by z, in upcoming loop)
B_prim = real(B0_vec .* exp(1i * omega * time));
B_prim(3,:) = -350 * 10^-9;

A = 1;
phi = 0;

mu0 = 4*pi*1e-7;  % vacuum permeability
R_E = 1560e3;     % Europa radius in meters

% induced dipole moment magnitude
M0 = -(4*pi/mu0) * A * exp(1i*phi) * B_prim .* (R_E^3)/2;
M_real = real(M0);


% define function to compute secondary magnetic field at a point
r_hat = @(r) r / norm(r);
B_sec = @(r, M) (mu0/(4*pi)) * (3*dot(r_hat(r), M)*r_hat(r) - M) / norm(r)^3;

E = [0; 0; 0]; % not using E-field yet, left here for boris method computation

% grid parameters
theta_list = 90:-45:-90;   % latitude
phi_list   = 0:120:359;     % longitude

alpha_list = 0:40:80;     % incidence angles relative to normal (in degrees)
beta_list  = 0:120:359;    % directions in tangent plane (in degrees)

% counting total vectors produced
num_points = length(theta_list) * length(phi_list);
num_directions = length(alpha_list) * length(beta_list);
N_total = num_points * num_directions;

positions = zeros(N_total, 3);
velocities = zeros(N_total, 3);

index = 1;

for i = 1:length(theta_list)
    theta = deg2rad(theta_list(i));

    for j = 1:length(phi_list)
        phi = deg2rad(phi_list(j));

        % sphere position
        x = R_E * cos(theta) * cos(phi);
        y = R_E * cos(theta) * sin(phi);
        z = R_E * sin(theta);
        point = [x, y, z];

        % normal vector (outward)
        n = [x, y, z] / R_E;

        % local tangent basis
        e_theta = [-sin(theta)*cos(phi), -sin(theta)*sin(phi), cos(theta)];
        e_phi   = [-sin(phi), cos(phi), 0];

        % loop over incidence angles and directions
        for k = 1:length(alpha_list)
            alpha = deg2rad(alpha_list(k));

            for a = 1:length(beta_list)
                beta = deg2rad(beta_list(a));

                % tilted vector (inward direction)
                v = cos(alpha)*(-n) + sin(alpha)*(cos(beta)*e_theta + sin(beta)*e_phi);

                % storing positions and velocities
                positions(index, :) = point;
                velocities(index, :) = v;

                index = index + 1;
            end
        end
    end
end

r = zeros(3, numsteps, N_total);     % position (m)
v = zeros(3, numsteps, N_total);     % velocity (m/s)

r(:, 1, :) = transpose(positions);  % final position
v(:, 1, :) = transpose(velocities) .* 10000000;   % final velocity

num_recorded = ones(N_total);

for p = 1:N_total
    for u = 1:numsteps-1
    
        B = B_sec(r(:,u,p), M_real(:,u)) + B_prim(:,u); 

        t = (q * B / m) * (0.5 * timestep);
        s = 2 * t / (1 + dot(t, t));
    
        v_minus = v(:,u,p) + 0.5 .* q ./ m .* timestep * E;
        
        % rotation
        v_prime = v_minus + cross(v_minus, t);
        v_plus  = v_minus + cross(v_prime, s);
    
        v_plus = v_plus + 0.5 * (q/m) * E * timestep;
        
        % advance position using v_plus (same as v at half-step)
        r(:,u+1,p) = r(:,u,p) + v_plus * timestep;
        
        % store updated velocity
        v(:,u+1,p) = v_plus;
        num_recorded(p) = num_recorded(p) + 1;

        % cut simulation if hitting Europa or flying away
        if norm(r(:,u+1,p)) < R_E || norm(r(:,u+1,p)) > 10 * R_E 
            break
        end
    end
end

% plotting
figure;
hold on
xlabel('x'); ylabel('y'); zlabel('z')
title('Particle Trajectories in Europa Magnetic Environment')

% adding Europa sphere
[xs, ys, zs] = sphere(50);  % higher resolution sphere
surf(R_E * xs, R_E * ys, R_E * zs, ...
    'FaceColor', [0.2 0.5 1], ...     % easy-ish to see blue
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.5);               % semi-transparent

% set view and axis
axis equal
grid on
view(3)

for i = 1:N_total
    plot3(r(1, 1:num_recorded(i), i), r(2, 1:num_recorded(i), i), r(3, 1:num_recorded(i), i))
end
hold off