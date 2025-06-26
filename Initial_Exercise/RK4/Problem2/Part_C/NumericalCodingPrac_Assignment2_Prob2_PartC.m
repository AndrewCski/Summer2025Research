%% Part C

q = 1.602 * 10 ^ -19;
m = 1.67262192 * 10 ^ -27;
B = [0; 0; 31000 * 10 ^ -9];
v0 = 9500000;
RE = 6371000;
y0 = RE * 8;
B0 = 31000 * 10 ^-9;
w = q / m * B(3);

timestep = 0.01;
numsteps = 10000;

radius = @(r) sqrt(r(1).^2 + r(2).^2 + r(3).^2);

% [x; y; z; vx; vy; vz]
state = zeros(6, numsteps);
state(:,1) = [y0; 0; 0; v0; v0; 0];  

% define function to compute dy/dt
Bfunc = @(r) (B0 * RE^3 / radius(r)^5) * ...
    [3 * r(1) * r(3);
     3 * r(2) * r(3);
     3 * r(3)^2 - radius(r)^2];
derivative = @(r, v) [v; q/m * cross(v, Bfunc(r))];

energies = .5 * m .* (state(4)^2 + state(5)^2 + state(6)^2);

for i = 1:numsteps-1
    % current position and velocity
    r = state(1:3,i);
    v = state(4:6,i);

    k1 = derivative(r, v);
    k2 = derivative(r + 0.5*timestep*k1(1:3), v + 0.5*timestep*k1(4:6));
    k3 = derivative(r + 0.5*timestep*k2(1:3), v + 0.5*timestep*k2(4:6));
    k4 = derivative(r + timestep*k3(1:3), v + timestep*k3(4:6));

    state(:,i+1) = state(:,i) + timestep/6 * (k1 + 2*k2 + 2*k3 + k4);
    energies(i + 1) = .5 * m .* (state(4,i+1)^2 + state(5,i+1)^2 + state(6,i+1)^2);
end

times = linspace(0, timestep * numsteps, numsteps);

timestep_str = strrep(sprintf('timestep_%g', timestep), '.', 'p'); % e.g., 'timestep0p001'
title_suffix = sprintf(' (Timestep = %.4f)', timestep);     % e.g., ' (Timestep = 0.001)'

figure(10)
plot3(state(1, :), state(2, :), state(3, :))
hold on
xlabel('x'); ylabel('y'); zlabel('z')
title(['Particle Trajectory in Dipole Magnetic Field' title_suffix])


% adding Earth sphere
[xs, ys, zs] = sphere(50);  % higher resolution sphere
surf(RE * xs, RE * ys, RE * zs, ...
    'FaceColor', [0.2 0.5 1], ...     % earth-like blue
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.5);               % semi-transparent

% set view and axis
axis equal
grid on
view(3)
saveas(gcf,['Varying_Bfield_xyz_vxvy' timestep_str '.jpg']);


figure(11)
plot(times, energies)
xlabel('time')
ylabel('energy')
title(['particle energy in varying Bfield' title_suffix])
saveas(gcf,['Varying_Bfield_energy_vxvy' timestep_str '.jpg']);


E0 = .5 * m * v0^2; % initial energy
perc_error = abs(energies - E0)/E0;

figure(12)
plot(times, perc_error)
xlabel('time')
ylabel('percentage error')
title(['varying Bfield energy error' title_suffix])
saveas(gcf,['Varying_Bfield_perc_error_vxvy' timestep_str '.jpg']);
