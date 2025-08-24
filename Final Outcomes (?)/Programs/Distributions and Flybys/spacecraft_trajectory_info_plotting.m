    % visualizing spacecraft trajectories!
    % Comment and uncomment parts of this code to change axes used.

R_E = 1560e3;                                % Europa radius (m)
omega = 2*pi/(11.23*3600);                   % synodic frequency (rad/s)
period_time = 0 * (2 * pi / omega)/4;

num_spots = 10000;
trajectory1 = zeros(num_spots, 3);
trajectory2 = zeros(num_spots, 3);
trajectory3 = zeros(num_spots, 3);
trajectory4 = zeros(num_spots, 3);
trajectory5 = zeros(num_spots, 3);

trajectory1(:,1) = linspace(-3*R_E + 50000, 3*R_E - 50000, num_spots);
trajectory2(:,1) = linspace(-3*R_E + 50000, 3*R_E - 50000, num_spots);
trajectory3(:,1) = linspace(-3*R_E + 50000, 3*R_E - 50000, num_spots);
trajectory4(:,1) = linspace(-3*R_E + 50000, 3*R_E - 50000, num_spots);
trajectory5(:,1) = linspace(-3*R_E + 50000, 3*R_E - 50000, num_spots);

trajectory1(:,2) = -25e3 - R_E;
trajectory2(:,2) = -1.05*R_E;
trajectory3(:,2) = -1.25*R_E;
trajectory4(:,2) = -1.5*R_E;
trajectory5(:,2) = -2.5*R_E;

%%

n10 = 4e7;
hs1 = 100e3;
n20 = 1e6;
hs2 = 500e3;

sim_range = linspace(-3*R_E, 3*R_E, 1000);
% choose which plane you want
[Y, Z] = meshgrid(sim_range, sim_range);
X = zeros(size(Y));

% distance and mask
r_mag = sqrt(X.^2 + Y.^2 + Z.^2);
mask = (r_mag >= R_E);

% compute phi in (-pi, pi]
phis = atan2(Y, X) + pi;
phis = mod(phis + pi, 2*pi) - pi;

% base density everywhere 
baseN = n10 .* exp(-(r_mag - R_E)/hs1) + n20 .* exp(-(r_mag - R_E)/hs2);

% apply the angular enhancement for |phi| <= pi/2
halfPlane = abs(phis) <= pi/2;  % picks +/-90° around -x
N_n = NaN(size(r_mag));         % NaN inside the body automatically

% combined mask for valid points
applyMask = mask & halfPlane;
N_n(mask) = baseN(mask);                     % default where outside body
N_n(applyMask) = (1 + 2 .* cos(phis(applyMask))) .* baseN(applyMask);

% prevent negative values if (1+2*cos) can go negative
N_n(N_n < 0) = 0;

% plume time!

% rand_theta = randi([-90, 90]);  % NOTE: NEED THETA = 0 FOR VIEWING ON XY,
% WILL NEED TO CHANGE PHI INSTEAD FOR VIEWING ON YZ (phi = 0/180?)
rand_theta = 15;
rand_theta = deg2rad(rand_theta);
rand_phi = 180;
rand_phi = deg2rad(rand_phi);

plume_r = [R_E * cos(rand_theta) * sin(rand_phi), ...
    R_E * cos(rand_theta) * cos(rand_phi), R_E * sin(rand_theta)];

Np0 = 2e9;
Hp = 150e3;
thetap = 15 * pi / 180;

% unit surface normal
r_s = plume_r / norm(plume_r); 

% vector from surface point to every grid point
Vx = X - plume_r(1);
Vy = Y - plume_r(2);
Vz = Z - plume_r(3);

% dot product with r_s 
dotRV = Vx * r_s(1) + Vy * r_s(2) + Vz * r_s(3);

% norm of V 
normV = sqrt(Vx.^2 + Vy.^2 + Vz.^2);

% polar angles
thetas = acos(dotRV ./ normV);
plume_indices = thetas < pi/2;

for i = 1:numel(N_n)
    if plume_indices(i)
        N_n(i) = N_n(i) + Np0 * exp((-r_mag(i) + R_E)/Hp) ...
            * exp(-(thetas(i)/thetap)^2);
    end
end

figure;
hold on
imagesc(sim_range, sim_range, N_n);
set(gca,'YDir','normal');
set(gca, 'ColorScale', 'log');
xlabel('y (m)'); ylabel('z (m)');
% p1 = plot(trajectory1(:,1), trajectory1(:,2), 'g', 'DisplayName', 'dist = 25 km');
% p2 = plot(trajectory2(:,1), trajectory2(:,2), 'm', 'DisplayName', 'dist = 78 km');
% p3 = plot(trajectory3(:,1), trajectory3(:,2), 'w', 'DisplayName', 'dist = 0.25 R_E');
% p4 = plot(trajectory4(:,1), trajectory4(:,2), 'y', 'DisplayName', 'dist = 0.5 R_E');
% p5 = plot(trajectory5(:,1), trajectory5(:,2), 'c', 'DisplayName', 'dist = 1.5 R_E');
p1 = plot(trajectory1(1,2), trajectory1(1,3), 'xg', 'DisplayName', 'dist = 25 km');
p2 = plot(trajectory2(1,2), trajectory2(1,3), 'xm', 'DisplayName', 'dist = 78 km');
p3 = plot(trajectory3(1,2), trajectory3(1,3), 'xw', 'DisplayName', 'dist = 0.25 R_E');
p4 = plot(trajectory4(1,2), trajectory4(1,3), 'xy', 'DisplayName', 'dist = 0.5 R_E');
p5 = plot(trajectory5(1,2), trajectory5(1,3), 'xc', 'DisplayName', 'dist = 1.5 R_E');
legend([p1 p2 p3 p4 p5], 'Location', 'best')
axis equal;
hold off