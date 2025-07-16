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
% [X, Y] = meshgrid(x_range, y_range);
% Z = z_plane * ones(size(X));

[X, Y] = meshgrid(y_range, z_range);
Z = zeros(size(X));


% distance from origin for each grid point
r_mag = sqrt(X.^2 + Y.^2 + Z.^2);

% define exclusion zone
r_min = R_E;  % exclude inside this radius
mask = (r_mag >= r_min);

filename = 'E_total_field_evolution_xy_plane.gif';

% loop over time indices
for t_idx = 1:100:numsteps   % e.g., take every 100th step for speed

    % B_prim and M_real at this time
    Bprim_now = B_prim(:, t_idx);
    M = M_real(:, t_idx);

    % compute field component
    B_comp = NaN(size(X));  % store B of choice for this frame
    E_comp = NaN(size(X));  % NaN => will not plot

    for i = 1:numel(X)
        if ~mask(i)
            continue;
        end

        r = [X(i); Y(i); Z(i)];
        Bsec = B_sec(r, M);
        Btotal = Bprim_now + Bsec;
        Etotal = -cross([90000; 0; 0], Btotal); 

        B_comp(i) = norm(Btotal);  % 1 = x, 2 = y, 3 = z, norm is total
        E_comp(i) = norm(Etotal);  % x = 1, y = 2, z = 3, norm is total

    end

    % plot this particular instance
    imagesc(x_range/1e3, y_range/1e3, E_comp);
    set(gca, 'YDir', 'normal');
    colorbar;
    clim([0 0.06]);   % fixing color scale for consistency (changing colors suck!)
    xlabel('x (km)');
    ylabel('y (km)');
    title(sprintf('Total E-field (Vm^-1) at t = %.2f hours', time(t_idx)/3600));
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