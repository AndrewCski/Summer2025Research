R_E = 1560e3;                                % Europa radius (m)
num_spots = 10000;
test_line = zeros(num_spots, 3);
test_line(:,2) = linspace(R_E, R_E+1500e3, num_spots);

neutral_densities = get_neutral_vol_density(test_line(:,1), test_line(:,2), test_line(:,3), NaN, NaN, 15, "cm");
O_ion_rates = neutral_densities .* 1e-6;
O2_ion_rates = neutral_densities .* 1e-7;

paper_e_density = 3.25e3*(exp((-test_line(:,2) + R_E)/440e3) + exp((-test_line(:,2) + R_E)/240e3));

tau = 60 + 80*exp((-test_line(:,2) + R_E)/1000e3);

electron_densities = (O_ion_rates + O2_ion_rates) .* tau;

%%

figure;
plot(O_ion_rates + O2_ion_rates, (test_line(:,2) - R_E)./1e3)
xlabel('electron density rates (cm^{-3}s^{-1})')
ylabel('altitude (km)')
title('Electron production rates in Europa exosphere')

%%

figure;
plot(electron_densities, (test_line(:,2) - R_E)./1e3)
xlabel('electron densities (cm^{-3})')
ylabel('altitude (km)')
title('Electron Densities in Europa exosphere')

%%

figure;
plot(paper_e_density, (test_line(:,2) - R_E)./1e3)
xlabel('electron densities (cm^{-3})')
ylabel('altitude (km)')
title('Electron Densities in Europa exosphere')



