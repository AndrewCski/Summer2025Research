% constants

q = 1.602e-19;              % charge of proton (C)
m = 1.67262192e-27;         % mass of proton (kg)

omega = 2*pi/(11.23*3600);  % Synodic frequency in rad/s

% simulation parameters
numsteps = 100000;                      % Number of time steps
period = 2 * pi / (100 * omega);        % modelling 1/100 of a synodic period
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

r = zeros(3, numsteps);     % position (m)
v = zeros(3, numsteps);     % velocity (m/s)

r(:,1) = [R_E; 0; 0];  % final position
v(:,1) = [-9500000; 0; 0];   % final velocity

for n = 1:numsteps-1

    % adding constant z-dir B field at end of B expression here
    B = B_sec(r(:,n), M_real(:,n)) + B_prim(:,n); 

    t = (q * B / m) * (0.5 * timestep);
    s = 2 * t / (1 + dot(t, t));

    v_minus = v(:,n) + 0.5 .* q ./ m .* timestep * E;
    
    % rotation
    v_prime = v_minus + cross(v_minus, t);
    v_plus  = v_minus + cross(v_prime, s);

    v_plus = v_plus + 0.5 * (q/m) * E * timestep;
    
    % advance position using v_plus (same as v at half-step)
    r(:,n+1) = r(:,n) + v_plus * timestep;
    
    % store updated velocity
    v(:,n+1) = v_plus;

    if norm(r(:,n+1)) < R_E
        break
    end
end

% plotting
figure;
plot3(r(1, 1:n+1), r(2, 1:n+1), r(3, 1:n+1))
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
hold off

energies = .5 .* m .* vecnorm(v(:, 1:n+1)).^2 * 6.242 * 10 ^ 18;
times = linspace(0, timestep * n+1, n+1);

figure;
plot(times, energies);
xlabel('time (s)');
ylabel('energy (eV)');
title('Proton Energy during flyby around Europa');
axis equal;
grid on;