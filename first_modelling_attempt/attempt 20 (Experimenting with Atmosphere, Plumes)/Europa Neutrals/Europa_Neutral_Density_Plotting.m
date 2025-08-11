R_E = 1560e3;     % Europa radius (m)
n10 = 4e7;
hs1 = 100e3;
n20 = 1e6;
hs2 = 500e3;

x_plane = 0;    % slice plane coordinate
z_plane = 0;
sim_range = linspace(-5*R_E, 5*R_E, 1000);
% choose which plane you want
[Y, Z] = meshgrid(sim_range, sim_range);
X = z_plane * ones(size(Y));

% distance and mask
r_mag = sqrt(X.^2 + Y.^2 + Z.^2);
mask = (r_mag >= R_E);

% compute azimuth phi in (-pi, pi]
phis = atan2(Y, X) + pi;
phis = mod(phis + pi, 2*pi) - pi;

% base density everywhere
baseN = n10 .* exp(-(r_mag - R_E)/hs1) + n20 .* exp(-(r_mag - R_E)/hs2);

% apply the angular enhancement for |phi| <= pi/2
halfPlane = abs(phis) <= pi/2;   % picks +/-90° around +x (or around -x if you shifted)
N_n = NaN(size(r_mag));         % NaN inside the body automatically

% combined mask for valid points
applyMask = mask & halfPlane;
N_n(mask) = baseN(mask);                     % default where outside body
N_n(applyMask) = (1 + 2 .* cos(phis(applyMask))) .* baseN(applyMask);

% prevent negative values if (1+2*cos) can go negative
N_n(N_n < 0) = 0;

figure;
imagesc(sim_range/1e3, sim_range/1e3, N_n);
set(gca,'YDir','normal');
cb = colorbar; cb.Label.String = 'parts/cm^-3';
xlabel('y (km)'); ylabel('z (km)');

title('Europa Neutral Density Distribution');
