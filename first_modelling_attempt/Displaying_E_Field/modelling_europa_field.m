% constants

q = 1.602e-19;              % charge of proton (C)
m = 1.67262192e-27;         % mass of proton (kg)

omega = 2*pi/(11.23*3600);  % Synodic frequency in rad/s

% simulation parameters
numsteps = 10000;               % Number of time steps
period = 2 * pi / omega;        % modelling a synodic period
time = linspace(0, period, numsteps);

B0_vec = [75 * 10 ^-9; 200 * 10 ^-9; 0];  % 75 nT azimuthal, 200 nT radial, using Zimmer coords
% (-350 nT constant added to North-South, represented by z, in upcoming loop)
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


% parameters for display
time_idx = 5000;   % index of the chosen timestep
x_plane = 0;    % x = 0 plane for the heatmap
y_plane = 0;    % y = 0 plane for the heatmap
z_plane = 0;    % z = 0 plane for the heatmap

% define grid 
x_range = linspace(-5*R_E, 5*R_E, 1000);  % adjust range and resolution
y_range = linspace(-5*R_E, 5*R_E, 1000);
z_range = linspace(-5*R_E, 5*R_E, 1000);

% CHANGE DEPENDING ON AXES USED...
% [X, Y] = meshgrid(x_range, y_range);
% Z = z_plane * ones(size(X));

[Y, Z] = meshgrid(x_range, z_range);
X = x_plane * ones(size(Y));

% B_prim at chosen time
Bprim_now = B_prim(:, time_idx);

% M_real at this time
M = M_real(:, time_idx);

% distance from origin for each grid point
r_mag = sqrt(X.^2 + Y.^2 + Z.^2);

% define exclusion zone
r_min = R_E;  % exclude inside this radius
mask = (r_mag >= r_min);

B_comp = NaN(size(X));  % NaN => will not plot
E_comp = NaN(size(X));  % NaN => will not plot

% loop over valid points only
for i = 1:numel(X)
    if ~mask(i)
        continue;
    end

    r = [X(i); Y(i); Z(i)];
    Bsec = B_sec(r, M);
    Btotal = Bprim_now + Bsec;
    Etotal = -cross([90000; 0; 0], Btotal); 
    B_comp(i) = norm(Btotal);  % x = 1, y = 2, z = 3, norm is total
    E_comp(i) = norm(Etotal);  % x = 1, y = 2, z = 3, norm is total
end

figure;
imagesc(y_range/1e3, z_range/1e3, E_comp);
set(gca,'YDir','normal'); % I think this is supposed to be Ydir???
colorbar;
xlabel('y (km)')
ylabel('z (km)')
title(sprintf('E (V m^-1) at time step %d/%d', time_idx, numsteps));