% constants

q = 1.602e-19;              % charge of proton (C)
m = 1.67262192e-27;         % mass of proton (kg)
c = 299792458;              % speed of light (m/s)

omega = 2*pi/(11.23*3600);  % Synodic frequency in rad/s

% simulation parameters
% numsteps = round(abs(period / timestep));   % number of time steps
numsteps = 10000;               % doesnt tend to get to above number
period = 2 * pi / omega;        % modelling a synodic period
time = linspace(0, period, numsteps);

% 100 timesteps per gyration, 3.65.. is computed average primary B field
timestep = -2 * pi * m / (50 * q * 3.659 * 10^-7); 

B0_vec = [75 * 10 ^-9; 200 * 10 ^-9; 0];
% 75 nT azimuthal, 200 nT radial, using Zimmer coords
% (-350 nT constant added to North-South (z) in upcoming loop)
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
phi_list   = 0:60:359;    % longitude

% alpha_list = 0;     % incidence angles relative to normal (in degrees)
% beta_list  = 0;    % directions in tangent plane (in degrees)
alpha_list = 0:40:80;     % incidence angles relative to normal (in degrees)
beta_list  = 0:60:359;    % directions in tangent plane (in degrees)

% creating protons of energies:
% 1 eV, 10 eV, 100 eV, 500 eV, 1keV, 10 keV, 100 keV, 1 MeV
% initial_velocities = [13841.22743, 43769.80428, 138412.2743, 309499.2542, 437698.0428,...
%     1384122.743, 4376980.428, 13841227.43];
initial_velocities = [437698.0428, 1384122.743, 4376980.428, 13841227.43];

% counting total vectors produced
num_points = length(theta_list) * length(phi_list);
num_directions = length(alpha_list) * length(beta_list);
N_total = num_points * num_directions * length(initial_velocities);

positions = zeros(N_total, 3);
velocities = zeros(N_total, 3);

index = 1;
for h = 1: length(initial_velocities)
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
                    velocities(index, :) = v .* initial_velocities(h);
    
                    index = index + 1;
                end
            end
        end
    end
end

r = zeros(3, numsteps, N_total);     % position (m)
v = zeros(3, numsteps, N_total);     % velocity (m/s)

r(:, 1, :) = transpose(positions);          % final position
v(:, 1, :) = transpose(velocities);   % final velocity

num_recorded = ones(N_total, 1); 
escaped = ones(N_total, 1); % tracking if each particle is traced back to Europa or not
impact_coords = zeros(N_total, 6); % coords of before and after going through surface

for p_index = 1:N_total
    for u = 1:numsteps-1
        r_now = r(:,u,p_index);

        % magnetic field at this position and time
        B = B_sec(r_now, M_real(:,u)) + B_prim(:,u);
        E_now = E;  % still zero for now

        % current velocity and momentum
        v_now = v(:,u,p_index);
        p_now = m * v_now / sqrt(1 - norm(v_now)^2 / c^2);  % relativistic momentum

        % step 1: half electric field push
        p_minus = p_now + (q * timestep / 2) * E_now;

        % step 2: magnetic rotation (Boris rotation in momentum space)
        gamma_minus = sqrt(1 + norm(p_minus)^2 / (m^2 * c^2));
        t = (q * timestep / (2 * m * gamma_minus)) * B;
        t_mag2 = dot(t, t);
        s = 2 * t / (1 + t_mag2);

        % p' = p_minus + p_minus x t
        p_prime = p_minus + cross(p_minus, t);
        % p_plus = p_minus + (p_prime x s)
        p_plus = p_minus + cross(p_prime, s);

        % step 3: second half electric field push
        p_new = p_plus + (q * timestep / 2) * E_now;

        % final momentum and velocity update
        gamma_new = sqrt(1 + norm(p_new)^2 / (m^2 * c^2));
        v_new = p_new / (gamma_new * m);

        % update position using updated velocity
        r(:,u+1,p_index) = r(:,u,p_index) + v_new * timestep;

        % store velocity
        v(:,u+1,p_index) = v_new;

        % track number of steps
        num_recorded(p_index) = num_recorded(p_index) + 1;

        % stop if particle hits Europa or escapes
        r_norm = norm(r(:,u+1,p_index));
        if r_norm < R_E || r_norm > 10 * R_E
            escaped(p_index) = 0;
            impact_coords(p_index, 1:3) = r(:,u,p_index);
            impact_coords(p_index, 4:6) = r(:,u+1,p_index);
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
    plot3(r(1, 1:num_recorded(i), i), r(2, 1:num_recorded(i), i), ...
        r(3, 1:num_recorded(i), i))
end
hold off



figure;
hold on
xlabel('x'); ylabel('y'); zlabel('z')
title('Particle Landing Sites')
surf(R_E * xs, R_E * ys, R_E * zs, ...
    'FaceColor', [0.2 0.5 1], ...     % easy-ish to see blue
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.5);               % semi-transparent

% set view and axis
axis equal
grid on
view(3)

interp_points = zeros(3, N_total);

for i = 1:N_total
    if escaped(i) == 0
        
        direction_vector = impact_coords(i,4:6) - impact_coords(i,1:3);

        % Solve the quadratic equation
        syms f
        equation = sum((impact_coords(i,1:3) + f * direction_vector).^2) == R_E^2;
        solutions_f = solve(equation, f);
        
        % Filter for valid solutions (0 <= t <= 1)
        valid_f = double(solutions_f(solutions_f >= 0 & solutions_f <= 1));
        
        % Retrieve the interpolated point(s)
        if ~isempty(valid_f)
            interpolated_point = impact_coords(i,1:3) + valid_f(1) * direction_vector; 
            % If multiple solutions, take the first one
            % plot3(interpolated_point(1), interpolated_point(2), interpolated_point(3),...
            %     'o', 'MarkerSize', 3, 'LineWidth', 3)
            interp_points(:, i) = interpolated_point;
        end
    end
end

for i = 1:N_total
    if norm(interp_points(:, i)) ~= 0
        scatter3(interp_points(1, i), interp_points(2, i), interp_points(3, i),...
            [],0.5 .* m .* initial_velocities(ceil(i/(num_points * ...
            num_directions))),'filled')
        disp(0.5 .* m .* initial_velocities(ceil(i/(num_points * num_directions)))...
            .^2 .* 6.242 .* 10.^18)
    end
end
colorbar
set(gca,'ColorScale','log')
hold off