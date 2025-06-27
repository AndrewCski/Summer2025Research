% Constants
q = 1.602e-19;              % Charge of proton (C)
m = 1.67262192e-27;         % Mass of proton (kg)
B = [0; 0; 31000 * 10^-9];  % Uniform magnetic field (T)
E = [0; .31; 0];

% Simulation parameters
timestep = 1e-5;                 % Time step (s)
numsteps = 10000;            % Number of time steps

% Initial conditions
r = zeros(3, numsteps);     % Position (m)
v = zeros(3, numsteps);     % Velocity (m/s)

v(:,1) = [1e5; 0; 0];       % Initial velocity
r(:,1) = [0; 6371000 * 8; 0];  % Initial position (8 earth radii)

t = (q * B / m) * (0.5 * timestep);       % t = (q B dt)/(2m)
s = 2 * t / (1 + dot(t, t));

for n = 1:numsteps-1
    v_minus = v(:,n) + 0.5 .* q ./ m .* timestep * E;
    
    % Rotation
    v_prime = v_minus + cross(v_minus, t);
    v_plus  = v_minus + cross(v_prime, s);

    v_plus = v_plus + 0.5 * (q/m) * E * timestep;
    
    % Advance position using v_plus (same as v at half-step)
    r(:,n+1) = r(:,n) + v_plus * timestep;
    
    % Store updated velocity
    v(:,n+1) = v_plus;
end

% Plotting
figure;
plot(r(1,:), r(2,:));
xlabel('x (m)');
ylabel('y (m)');
title('Proton Gyromotion via Boris Method');
axis equal;
grid on;

energies = .5 .* m .* vecnorm(v).^2;
times = linspace(0, timestep * numsteps, numsteps);

figure;
plot(times, energies);
xlabel('time (s)');
ylabel('energy (joules)');
title('Proton Gyromotion energy');
axis equal;
grid on;