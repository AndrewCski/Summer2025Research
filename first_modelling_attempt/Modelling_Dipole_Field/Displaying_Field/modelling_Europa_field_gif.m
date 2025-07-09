% constants

q = 1.602e-19;              % charge of proton (C)
m = 1.67262192e-27;         % mass of proton (kg)

omega = 2*pi/(11.23*3600);  % synodic frequency in rad/s

% simulation parameters
numsteps = 10000;               % number of time steps
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
x_plane = 0;    % x = 0 plane for the heatmap
y_plane = 0;    % y = 0 plane for the heatmap
z_plane = 0;    % z = 0 plane for the heatmap

% define grid 
x_range = linspace(-5*R_E, 5*R_E, 1000);  % adjust range and resolution
y_range = linspace(-5*R_E, 5*R_E, 1000);
z_range = linspace(-5*R_E, 5*R_E, 1000);

% CHANGE DEPENDING ON AXES USED...
% [BX, BY] = meshgrid(x_range, y_range);
% BZ = z_plane * ones(size(BX));

[BY, BZ] = meshgrid(y_range, z_range);
BX = zeros(size(BY));


% distance from origin for each grid point
r_mag = sqrt(BX.^2 + BY.^2 + BZ.^2);

% define exclusion zone
r_min = R_E;  % exclude inside this radius
mask = (r_mag >= r_min);

filename = 'Bz_field_evolution_yz_plane.gif';

% loop over time indices
for t_idx = 1:100:numsteps   % e.g., take every 100th step for speed

    % B_prim and M_real at this time
    Bprim_now = B_prim(:, t_idx);
    M = M_real(:, t_idx);

    % compute field component
    B_comp = NaN(size(BX));  % store B of choice for this frame
    for i = 1:numel(BX)
        if ~mask(i)
            continue;
        end

        r = [BX(i); BY(i); BZ(i)];
        Bsec = B_sec(r, M);
        Btotal = Bprim_now + Bsec;
        B_comp(i) = Btotal(3);  % 1 = x, 2 = y, 3 = z
    end

    % plot this particular instance
    imagesc(y_range/1e3, z_range/1e3, B_comp*1e9);
    set(gca, 'YDir', 'normal');
    colorbar;
    clim([-500 500]);   % fixing color scale for consistency (changing colors suck!)
    xlabel('y (km)');
    ylabel('z (km)');
    title(sprintf('B_z (nT) at t = %.2f hours', time(t_idx)/3600));
    drawnow;

    % capture frame and write to gif
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    if t_idx == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
end