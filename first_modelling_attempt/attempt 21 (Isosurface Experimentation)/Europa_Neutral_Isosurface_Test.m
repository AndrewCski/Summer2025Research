R_E = 1560e3;     % Europa radius (m)
n10 = 4e7;
hs1 = 100e3;
n20 = 1e6;
hs2 = 500e3;

sim_range = linspace(-5*R_E, 5*R_E, 100);
% choose which plane you want
[X, Y, Z] = meshgrid(sim_range, sim_range, sim_range);

% distance and mask
r_mag = sqrt(X.^2 + Y.^2 + Z.^2);
mask = (r_mag >= R_E);

% compute azimuth phi in (-pi, pi]
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
rand_theta = 0;
% rand_phi = randi([0, 359]);
rand_phi = 0;
rand_phi = mod(rand_phi, 2*pi) - pi;

plume_r = [R_E * cos(rand_theta) * cos(rand_phi), ...
    R_E * cos(rand_theta) * sin(rand_phi), R_E * sin(rand_theta)];

Np0 = 2e9;
Hp = 150e3;
thetap = 15 * pi / 180;

% unit surface normal
r_s = plume_r / norm(plume_r);   % 1x3

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
bar = 3e7;
p1 = patch(isosurface(X./R_E, Y./R_E, Z./R_E, N_n, bar));
p1.FaceColor = 'red';
p1.EdgeColor = 'none';
p1.FaceAlpha = 0.75; 

hold on;

% plot Europa
[XS, YS, ZS] = sphere(50); 
p2 = surf(XS, YS, ZS);
set(p2, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.25);
axis equal;
xlabel('x (R_E)');
ylabel('y (R_E)');
zlabel('z (R_E)');
title(sprintf('Europa Neutral Density Distribution Above %d cm^{-3}', bar));
hold off; 