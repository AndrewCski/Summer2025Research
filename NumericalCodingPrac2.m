%% #1: Gyromotion and trying RK4

numsteps = 100000;
timestep = 0.0001;

q = 1;
m = 1;
B = [0; 0; 1];
w = 1;
v0 = 1;
y0 = 1;

% [x; y; z; vx; vy; vz]
y = zeros(6, numsteps);
y(:,1) = [0; y0; 0; v0; 0; 0];  % position (0,1,0), velocity (1,0,0)

% define function to compute dy/dt
dydt = @(r, v) [v; q/m * cross(v, B)];

energies = .5 * m .* (y(4)^2 + y(5)^2 + y(6)^2);

for i = 1:numsteps-1
    % current position and velocity
    r = y(1:3,i);
    v = y(4:6,i);

    k1 = dydt(r, v);
    k2 = dydt(r + 0.5*timestep*k1(1:3), v + 0.5*timestep*k1(4:6));
    k3 = dydt(r + 0.5*timestep*k2(1:3), v + 0.5*timestep*k2(4:6));
    k4 = dydt(r + timestep*k3(1:3), v + timestep*k3(4:6));

    y(:,i+1) = y(:,i) + timestep/6 * (k1 + 2*k2 + 2*k3 + k4);
    energies(i + 1) = .5 * (y(4,i+1)^2 + y(5,i+1)^2 + y(6,i+1)^2);
end

times = linspace(0, timestep * numsteps, numsteps);

% radius over time
r_numerical = sqrt(y(1,:).^2 + y(2,:).^2);

% analytical solution
analytical_x = sin(w .* times) .* v0 ./ w;
analytical_y = cos(w .* times) .* v0 ./ w;
analytical_r = sqrt(analytical_x.^2 + analytical_y.^2);

figure(1)
hold on
plot(times, r_numerical)
plot(times, analytical_r, '--k')
legend('numerical r', 'analytical r')
xlabel('Time')
hold off

figure(2)

hold on
plot(times, y(1, :))
plot(times, y(2, :))
plot(times, analytical_x)
plot(times, analytical_y)
legend('numerical x', 'numerical y', 'analytical x', 'analytical y')
xlabel('Time')
hold off

figure(3)
plot(times, energies) % energy NOT conserved!

%% #2: Add E field...

numsteps = 1000;
timestep = 0.00001;

q = 1.602 * 10 ^ -19;
m = 1.67262192 * 10 ^ -27;
B = [0; 0; 31000 * 10 ^ -9];
v0 = 100000;
y0 = 6371000 * 8;

% doing y E field first

E = [0; 0.31; 0];

% [x; y; z; vx; vy; vz]
y = zeros(6, numsteps);
y(:,1) = [0; y0; 0; v0; 0; 0];  

% define function to compute dy/dt
dydt = @(r, v) [v; q/m * (E + cross(v, B))];

energies = .5 * m .* (y(4)^2 + y(5)^2 + y(6)^2);

for i = 1:numsteps-1
    % current position and velocity
    r = y(1:3,i);
    v = y(4:6,i);

    k1 = dydt(r, v);
    k2 = dydt(r + 0.5*timestep*k1(1:3), v + 0.5*timestep*k1(4:6));
    k3 = dydt(r + 0.5*timestep*k2(1:3), v + 0.5*timestep*k2(4:6));
    k4 = dydt(r + timestep*k3(1:3), v + timestep*k3(4:6));

    y(:,i+1) = y(:,i) + timestep/6 * (k1 + 2*k2 + 2*k3 + k4);
    energies(i + 1) = .5 * m .* (y(4,i+1)^2 + y(5,i+1)^2 + y(6,i+1)^2);
end

times = linspace(0, timestep * numsteps, numsteps);

figure(4)
plot(y(1, :), y(2, :))
xlabel('x')
ylabel('y')

figure(5)
plot(times, y(4, :))
xlabel('time')
ylabel('vx')

figure(6)
plot(times, energies)
xlabel('time')
ylabel('energy')

%% now trying x E field

E = [3.1; 0; 0];

% [x; y; z; vx; vy; vz]
y = zeros(6, numsteps);
y(:,1) = [0; y0; 0; v0; 0; 0];

% define function to compute dy/dt
dydt = @(r, v) [v; q/m * (E + cross(v, B))];

energies = .5 * m .* (y(4)^2 + y(5)^2 + y(6)^2);

for i = 1:numsteps-1
    % current position and velocity
    r = y(1:3,i);
    v = y(4:6,i);

    k1 = dydt(r, v);
    k2 = dydt(r + 0.5*timestep*k1(1:3), v + 0.5*timestep*k1(4:6));
    k3 = dydt(r + 0.5*timestep*k2(1:3), v + 0.5*timestep*k2(4:6));
    k4 = dydt(r + timestep*k3(1:3), v + timestep*k3(4:6));

    y(:,i+1) = y(:,i) + timestep/6 * (k1 + 2*k2 + 2*k3 + k4);
    energies(i + 1) = .5 * m .* (y(4,i+1)^2 + y(5,i+1)^2 + y(6,i+1)^2);
end

times = linspace(0, timestep * numsteps, numsteps);

figure(7)
plot(y(1, :), y(2, :))
xlabel('x')
ylabel('y')

figure(8)
plot(times, y(5, :))
xlabel('time')
ylabel('vy')

figure(9)
plot(times, energies)
xlabel('time')
ylabel('energy')

%% Part C

numsteps = 1000000;
timestep = 0.00001;

q = 1.602 * 10 ^ -19;
m = 1.67262192 * 10 ^ -27;
B = [0; 0; 31000 * 10 ^ -9];
v0 = 9500000;
RE = 6371000;
y0 = RE * 8;
B0 = 31000 * 10 ^-9;
E = [0; 0; 0];

radius = @(r) sqrt(r(1).^2 + r(2).^2 + r(3).^2);

% [x; y; z; vx; vy; vz]
y = zeros(6, numsteps);
y(:,1) = [y0; 0; 0; v0; 0; 0];  

% define function to compute dy/dt
Bfunc = @(r) [-3 * B0 * r(1) * r(2) / ((radius(r)/RE)^3 * radius(r)^2);
    -3 * B0 * r(2) * r(3) / ((radius(r)/RE)^3 * radius(r)^2);
    B0 * (r(1)^2 + r(2)^2 - 2 * r(3)^2) / ((radius(r)/RE)^3 * radius(r)^2)];
dydt = @(r, v) [v; q/m * (E + cross(v, Bfunc(r)))];

energies = .5 * m .* (y(4)^2 + y(5)^2 + y(6)^2);

for i = 1:numsteps-1
    % current position and velocity
    r = y(1:3,i);
    v = y(4:6,i);

    k1 = dydt(r, v);
    k2 = dydt(r + 0.5*timestep*k1(1:3), v + 0.5*timestep*k1(4:6));
    k3 = dydt(r + 0.5*timestep*k2(1:3), v + 0.5*timestep*k2(4:6));
    k4 = dydt(r + timestep*k3(1:3), v + timestep*k3(4:6));

    y(:,i+1) = y(:,i) + timestep/6 * (k1 + 2*k2 + 2*k3 + k4);
    energies(i + 1) = .5 * m .* (y(4,i+1)^2 + y(5,i+1)^2 + y(6,i+1)^2);
end

times = linspace(0, timestep * numsteps, numsteps);

figure(10)
plot3(y(1, :), y(2, :), y(3,:))
xlabel('x')
ylabel('y')
zlabel('z')

figure(11)
plot(times, energies)
xlabel('time')
ylabel('energy')