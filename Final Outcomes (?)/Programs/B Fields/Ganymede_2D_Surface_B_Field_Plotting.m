% program for creating a 2d map of Ganymede's B field at its surface, plotted
% respect to latitude and longitude (0 degrees is towards Jupiter, 90 degrees
% is towards Ganymede orbital direction, etc.). To change time in synodic
% period, change the factor in front of the "time" variable

omega = 2*pi/(10.53*3600);                            % synodic frequency (rad/s)
R_G = 2631e3;                                         % Ganymede radius (m)
mu0 = 4*pi*1e-7;           % vacuum permeability (H/m)

measure_gap = 1;

thetas = (-90 + measure_gap/2):measure_gap:(90 - measure_gap/2);
phis = (measure_gap/2):measure_gap:(360 - measure_gap/2);

time = 1/4 * (2 * pi / omega);                       % point in synodic period to start from

B0_vec = [15 * 10 ^-9; 80 * 10 ^-9; 0];
case_id = "Ganymede";
B_prim_now = [B0_vec(1) * cos(omega * time); B0_vec(2) * -sin(omega * time - pi/4.5); 0];
M_dip = [-4.1; 9; -131] * 1e18;
M_ind = ((-2 * pi * 0.84 * R_G^3) / mu0) .* B_prim_now;
M0 = M_dip + M_ind;
B_prim_now(3) = -75 * 10^-9;

B_x = zeros(length(thetas), length(phis));
B_y = zeros(length(thetas), length(phis));
B_z = zeros(length(thetas), length(phis));
B_tot = zeros(length(thetas), length(phis));

index = 1;
for i = 1:length(thetas)
    for j = 1:length(phis)
        x = R_G * cos(deg2rad(thetas(i))) * sin(deg2rad(phis(j)));
        y = R_G * cos(deg2rad(thetas(i))) * cos(deg2rad(phis(j)));
        z = R_G * sin(deg2rad(thetas(i)));
        r = [x; y; z];

        B_sec = (mu0/(4*pi)) .* (3.*dot(r, M0) .* r - ...
            M0 .* norm(r).^2) ./ norm(r)^5;
        B_x(i,j)   = B_prim_now(1) + B_sec(1);
        B_y(i,j)   = B_prim_now(2) + B_sec(2);
        B_z(i,j)   = B_prim_now(3) + B_sec(3);
        B_tot(i,j) = norm(B_prim_now + B_sec);

        index = index + 1;
    end
end

plot_thetas = -90:measure_gap:90;
plot_phis = 0:measure_gap:360;

figure;
clf;
imagesc(plot_phis, plot_thetas, B_x.*1e9);              
set(gca,'YDir','normal'); 
ylabel('Latitude (Degrees)')                              
xlabel('Longitude (Degrees)')           
ax=xticklabels;
xticklabels(flip(ax))
colormap(jet);
cb = colorbar;
cb.Label.String = 'B_x (nT)';
title(sprintf('Ganymede Surface B_x Map, at %.2f of Synodic Period', time / (2 * pi / omega)));

figure;
clf;
imagesc(plot_phis, plot_thetas, B_y.*1e9);              
set(gca,'YDir','normal'); 
ylabel('Latitude (Degrees)')                              
xlabel('Longitude (Degrees)')           
ax=xticklabels;
xticklabels(flip(ax))
colormap(jet);
cb = colorbar;
cb.Label.String = 'B_y (nT)';
title(sprintf('Ganymede Surface B_y Map, at %.2f of Synodic Period', time / (2 * pi / omega)));

figure;
clf;
imagesc(plot_phis, plot_thetas, B_z.*1e9);              
set(gca,'YDir','normal'); 
ylabel('Latitude (Degrees)')                              
xlabel('Longitude (Degrees)')           
ax=xticklabels;
xticklabels(flip(ax))
colormap(jet);
cb = colorbar;
cb.Label.String = 'B_z (nT)';
title(sprintf('Ganymede Surface B_z Map, at %.2f of Synodic Period', time / (2 * pi / omega)));

figure;
clf;
imagesc(plot_phis, plot_thetas, B_tot.*1e9);              
set(gca,'YDir','normal'); 
ylabel('Latitude (Degrees)')                              
xlabel('Longitude (Degrees)')           
ax=xticklabels;
xticklabels(flip(ax))
colormap(jet);
cb = colorbar;
cb.Label.String = 'B_{tot} (nT)';
title(sprintf('Ganymede Surface B_{tot} Map, at %.2f of Synodic Period', time / (2 * pi / omega)));