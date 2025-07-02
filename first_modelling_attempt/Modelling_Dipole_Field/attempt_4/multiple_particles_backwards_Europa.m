% constants

q = 1.602e-19;              % charge of proton (C)
m = 1.67262192e-27;         % mass of proton (kg)

omega = 2*pi/(11.23*3600);  % Synodic frequency in rad/s

% simulation parameters
numsteps = 100000;                      % Number of time steps
period = 2 * pi / (2000 * omega);        % modelling 1/100 of a synodic period
time = linspace(0, period, numsteps);
timestep = -period / numsteps;

B0_vec = [75 * 10 ^-9; 200 * 10 ^-9; 0];  % 75 nT azimuthal, 200 nT radial, using Zimmer coords
% (-350 nT constant added to North-South, represented by z, in upcoming loop)
B_prim = real(B0_vec .* exp(1i * omega * time));
B_prim(3,:) = -350 * 10^-19;

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

dTheta = 30; % degree step for latitude
dPhi = 60;   % degree step for longitude

% Create grid
theta_deg = 90:-dTheta:-90;  % latitude from +90 to -90
phi_deg = 0:dPhi:355;        % longitude from 0 to 355
% theta_deg = 90:-dTheta:80;  % testing lower number
% phi_deg = 0:dPhi:10;        % testing lower number

% Initialize arrays
num_points = length(theta_deg) * length(phi_deg);
points = zeros(num_points, 3);
vectors = zeros(num_points, 3);

index = 1;
for i = 1:length(theta_deg)
    theta = deg2rad(theta_deg(i));
    for j = 1:length(phi_deg)
        phi = deg2rad(phi_deg(j));
        
        % Cartesian coordinates
        x = R_E * cos(theta) * cos(phi);
        y = R_E * cos(theta) * sin(phi);
        z = R_E * sin(theta);
        
        % Position
        points(index, :) = [x, y, z];
        
        % Radially inward unit vector
        vec = -[x, y, z] / R_E; % since ||r|| = R
        
        vectors(index, :) = vec;
        
        index = index + 1;
    end
end

r = zeros(3, numsteps, num_points);     % position (m)
v = zeros(3, numsteps, num_points);     % velocity (m/s)

r(:, 1,:) = transpose(points);  % final position
v(:, 1, :) = transpose(vectors) .* 10000000;   % final velocity

num_recorded = ones(num_points);

for j = 1:num_points
    for n = 1:numsteps-1
    
        % adding constant z-dir B field at end of B expression here
        B = B_sec(r(:,n,j), M_real(:,n)) + B_prim(:,n); 
    
        t = (q * B / m) * (0.5 * timestep);
        s = 2 * t / (1 + dot(t, t));
    
        v_minus = v(:,n,j) + 0.5 .* q ./ m .* timestep * E;
        
        % rotation
        v_prime = v_minus + cross(v_minus, t);
        v_plus  = v_minus + cross(v_prime, s);
    
        v_plus = v_plus + 0.5 * (q/m) * E * timestep;
        
        % advance position using v_plus (same as v at half-step)
        r(:,n+1,j) = r(:,n,j) + v_plus * timestep;
        
        % store updated velocity
        v(:,n+1,j) = v_plus;

        num_recorded(j) = num_recorded(j) + 1;

        % cut simulation if hitting Europa or flying away
        if norm(r(:,n+1,j)) < R_E || norm(r(:,n+1,j)) > 10 * R_E 
            break
        end
    end
end

% plotting
figure;
hold on
xlabel('x'); ylabel('y'); zlabel('z')
title('Particle Trajectory in Europa Magnetic Environment')

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

for i = 1:num_points
    plot3(r(1, 1:num_recorded(i), i), r(2, 1:num_recorded(i), i), r(3, 1:num_recorded(i), i))
end
hold off