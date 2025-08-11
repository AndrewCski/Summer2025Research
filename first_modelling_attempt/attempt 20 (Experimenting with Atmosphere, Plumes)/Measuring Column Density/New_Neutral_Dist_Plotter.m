R_E = 1560e3;     % Europa radius (m)
n10 = 4e7;
hs1 = 100e3;
n20 = 1e6;
hs2 = 500e3;

sim_range = linspace(-5*R_E, 5*R_E, 2000);
% choose which plane you want
[X, Y] = meshgrid(sim_range, sim_range);
Z = zeros(size(Y));

N_n = get_neutral_vol_density(X, Y, Z, 0, 0, 15, "cm");

figure;
imagesc(sim_range/R_E, sim_range/R_E, N_n);
set(gca,'YDir','normal');
cb = colorbar; cb.Label.String = 'parts/cm^-3'; %set(gca, 'ColorScale', 'log'); clim([1e9 1e20]);
xlabel('z (R_E)'); ylabel('y (R_E)');
title('Europa Neutral Density Distribution');