% This file is used for creating images of the magnetic field or electric 
% field around Europa during any point in its synodic period.
%
% To change the axes of the image, alter the lines labelled with #1, #2, etc.

% constant
omega = 2*pi/(10.53*3600);                            % synodic frequency (rad/s)

% simulation parameters
numsteps = 10000;               % Number of time steps
period = 2 * pi / omega;        % modelling a synodic period
time = linspace(0, period, numsteps);

B0_vec = [-15 * 10 ^-9; 24 * 10 ^-9; 0];  
B_prim = real(B0_vec .* exp(0 * omega * time));

mu0 = 4*pi*1e-7;  % vacuum permeability
R_G = 2631e3;     % Ganymede radius (m)

% induced dipole moment magnitude

M_dip = [-4.1; 9; -131] * 1e18;
M_ind = ((-2 * pi * 0.84 * R_G^3) / mu0) .* B_prim(:,1);
M_real = M_dip + M_ind;
B_prim(3,:) = -75 * 10^-9;

% define function to compute secondary magnetic field at a point
r_hat = @(r) r / norm(r);
B_sec = @(r, M) (mu0/(4*pi)) * (3*dot(r_hat(r), M)*r_hat(r) - M) / norm(r)^3;

% parameters for display
time_idx = 1;   % index of the chosen timestep
x_plane = 0;    % x = 0 plane for the heatmap
y_plane = 0;    % y = 0 plane for the heatmap
z_plane = 0;    % z = 0 plane for the heatmap

% define grid 
x_range = linspace(-5*R_G, 5*R_G, 1000);  % adjust range and resolution
y_range = linspace(-5*R_G, 5*R_G, 1000);
z_range = linspace(-5*R_G, 5*R_G, 1000);

% CHANGE DEPENDING ON AXES USED...
% [X, Y] = meshgrid(x_range, y_range);
% Z = z_plane * ones(size(X));

[X, Z] = meshgrid(x_range, z_range);                               % #1 !!!
Y = x_plane * ones(size(Y));                                       % #2 !!!

% B_prim at chosen time
Bprim_now = B_prim(:, time_idx);

% M_real at this time
M = M_real(:, time_idx);

% distance from origin for each grid point
r_mag = sqrt(X.^2 + Y.^2 + Z.^2);

% define exclusion zone
r_min = R_G;  % exclude inside this radius
mask = (r_mag >= r_min);

B_x = NaN(size(X));  % NaN -> will not plot       % all of these are #3 !!!
B_y = NaN(size(X));
B_z = NaN(size(X));
B_tot = NaN(size(X));
% E_x = NaN(size(X)); 
% E_y = NaN(size(X)); 
% E_z = NaN(size(X)); 
% E_tot = NaN(size(X)); 

% loop over valid points only
for i = 1:numel(X)
    if ~mask(i)
        continue;
    end

    r = [X(i); Y(i); Z(i)];
    Bsec = B_sec(r, M);
    Btotal = Bprim_now + Bsec;
    % Etotal = -cross([90000; 0; 0], Btotal);     % all of these are #4 !!!
    B_x(i) = Btotal(1);  % x = 1, y = 2, z = 3, norm is total
    B_y(i) = Btotal(2); 
    B_z(i) = Btotal(3); 
    % B_tot(i) = norm(Btotal); 
    % E_x(i) = Etotal(1);  % x = 1, y = 2, z = 3, norm is total
    % E_y(i) = Etotal(2); 
    % E_z(i) = Etotal(3); 
    % E_tot(i) = norm(Etotal); 
end

figure;
imagesc(x_range/1e3, z_range/1e3, B_x * 1e9);                      % #5 !!!
set(gca,'YDir','normal'); % I think this is supposed to be Ydir???
colormap jet;
cb = colorbar;
cb.Label.String = 'B_x (nT)';
clim([-300, 100]);
xlabel('x (km)')                                                   % #6 !!!
ylabel('z (km)')                                                   % #7 !!!
title('Ganymede Magnetic Environment');
hold on
str_slc = streamslice(X/1e3, Z/1e3, B_x, B_z);                     % #8 !!!
set(str_slc,'Color','b');
hold off