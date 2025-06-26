%% #1: Gyromotion and trying RK4

q = 1;
m = 1;
B = [0; 0; 1];
w = 1;
v0 = 1;
y0 = 1;

timestep = 0.001;
numsteps = round(5 * 2 * pi / w / timestep); % model 5 cycles

% [x; y; z; vx; vy; vz]
state = zeros(6, numsteps);
state(:,1) = [0; y0; 0; v0; 0; 0];  % position (0,1,0), velocity (1,0,0)

% define function to compute dy/dt
dydt = @(r, v) [v; q/m * cross(v, B)];

energies = .5 * m .* (state(4)^2 + state(5)^2 + state(6)^2);

for i = 1:numsteps-1
    % current position and velocity
    r = state(1:3,i);
    v = state(4:6,i);

    k1 = dydt(r, v);
    k2 = dydt(r + 0.5*timestep*k1(1:3), v + 0.5*timestep*k1(4:6));
    k3 = dydt(r + 0.5*timestep*k2(1:3), v + 0.5*timestep*k2(4:6));
    k4 = dydt(r + timestep*k3(1:3), v + timestep*k3(4:6));

    state(:,i+1) = state(:,i) + timestep/6 * (k1 + 2*k2 + 2*k3 + k4);
    energies(i + 1) = .5 * (state(4,i+1)^2 + state(5,i+1)^2 + state(6,i+1)^2);
end

times = linspace(0, timestep * numsteps, numsteps);

% radius over time
r_numerical = sqrt(state(1,:).^2 + state(2,:).^2);

% analytical solution
analytical_x = sin(w .* times) .* v0 ./ w;
analytical_y = cos(w .* times) .* v0 ./ w;
analytical_r = sqrt(analytical_x.^2 + analytical_y.^2);

timestep_str = strrep(sprintf('ts_%g', timestep), '.', 'p'); % e.g., 'dt0p001'
title_suffix = sprintf(' (Timestep = %.4f)', timestep);     % e.g., ' (Timestep = 0.001)'

figure(1)
hold on
plot(times, r_numerical)
plot(times, analytical_r, '--k')
legend('numerical r', 'analytical r')
xlabel('Time')
title(['Radial Distance vs Time' title_suffix])
saveas(gcf, ['Gyromotion_AnalVSNumer_r_' timestep_str '.jpg']);
hold off

figure(2)
hold on
plot(times, state(1, :))
plot(times, state(2, :))
plot(times, analytical_x)
plot(times, analytical_y)
legend('numerical x', 'numerical y', 'analytical x', 'analytical y')
xlabel('Time')
title(['X/Y vs Time' title_suffix])
saveas(gcf, ['Gyromotion_AnalVSNumer_xy_' timestep_str '.jpg']);
hold off

figure(3)
plot(times, energies)
xlabel('Time')
ylabel('Energy')
title(['Energy vs Time' title_suffix])
saveas(gcf, ['Gyromotion_energy_' timestep_str '.jpg']);

perc_error = abs(r_numerical - analytical_r) ./ analytical_r;
figure(4)
plot(times, perc_error)
xlabel('Time')
ylabel('Percentage Error')
title(['Radius Perc-Error Over Time' title_suffix])
saveas(gcf, ['Perc_Error_' timestep_str '.jpg']);
