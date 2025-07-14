function [] = Europa_Particle_Energy_Escapes(escaped, energies)

% This function file is used for creating a plot depicting the rates at
% which various energies end up "escaping" Europa (I'm probably modelling 
% backwards!).
%
% The two parameters that must be passed are an 'escaped' array of 1's and
% 0's, depicting whether particles escaped Europa (1) or not (0), in
% addition to the energies of each of the particles (in eV)

[unique_energies, ~, idx] = unique(energies);
counts = accumarray(idx(:), escaped(:), [], @sum);

figure;
x_labels = string(round(unique_energies));
bar(x_labels, counts)
title('Escaped Particles Per Energy')
xlabel('Energy (eV)');
ylabel('Counts');
