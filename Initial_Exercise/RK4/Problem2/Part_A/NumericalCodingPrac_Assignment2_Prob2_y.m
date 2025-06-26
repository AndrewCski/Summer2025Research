%% #2: Add E field in y direction...

q = 1.602 * 10 ^ -19;
m = 1.67262192 * 10 ^ -27;
B = [0; 0; 31000 * 10 ^ -9];
E = [0; 0.31; 0];
v0 = 100000;
y0 = 6371000 * 8;
w = q / m * B(3);

timestep = 0.00001;
numsteps = round(10 * 2 * pi / (w * timestep)); % model 5 cycles

% [x; y; z; vx; vy; vz]
state = zeros(6, numsteps);
state(:,1) = [0; y0; 0; v0; 0; 0];  

% define function to compute dr/dt and d^2r/dt^2
derivative = @(r, v) [v; q/m * (E + cross(v, B))];

energies = zeros(1, numsteps);
energies(1) = 0.5 * m * sum(state(4:6,1).^2);

for i = 1:numsteps-1
    % current position and velocity
    r = state(1:3,i);
    v = state(4:6,i);

    k1 = derivative(r, v);
    k2 = derivative(r + 0.5*timestep*k1(1:3), v + 0.5*timestep*k1(4:6));
    k3 = derivative(r + 0.5*timestep*k2(1:3), v + 0.5*timestep*k2(4:6));
    k4 = derivative(r + timestep*k3(1:3), v + timestep*k3(4:6));

    state(:,i+1) = state(:,i) + timestep/6 * (k1 + 2*k2 + 2*k3 + k4);
    energies(i + 1) = 0.5 * m * sum(state(4:6,i+1).^2);
end

times = linspace(0, timestep * numsteps, numsteps);

timestep_str = strrep(sprintf('timestep_%g', timestep), '.', 'p'); % e.g., 'timestep0p001'
title_suffix = sprintf(' (Timestep = %.4f)', timestep);     % e.g., ' (Timestep = 0.001)'

figure(4)
plot(state(1, :), state(2, :))
xlabel('x')
ylabel('y')
title(['particle in Y-dir Efield' title_suffix])
saveas(gcf,['10Cycle_YDir_Efield_xVSy' timestep_str '.jpg']);

figure(5)
plot(times, state(4, :))
xlabel('time')
ylabel('vx')
title(['x-velocity of particle in Y-dir Efield' title_suffix])
saveas(gcf,['10Cycle_YDir_Efield_vx' timestep_str '.jpg']);


figure(6)
plot(times, energies)
xlabel('time')
ylabel('energy')
title(['energy of particle in Y-dir Efield' title_suffix])
saveas(gcf,['10Cycle_YDir_Efield_energy' timestep_str '.jpg']);


% drift velocity
v_drift = cross(E, B) / norm(B)^2;

% compute exact vx at each time step
vx_exact = v_drift(1) + (v0 - v_drift(1)) * cos(w * times);

% Compute absolute or relative error
vx_numerical = state(4, :);
vx_perc_error = abs(vx_numerical - vx_exact) ./ vx_exact;

figure(7)
plot(times, vx_perc_error)
xlabel('time')
ylabel('percent error in vx')
title(['Numerical Error in vx vs Analytical Solution' title_suffix])
saveas(gcf,['10Cycle_YDir_Efield_vx_error' timestep_str '.jpg']);
