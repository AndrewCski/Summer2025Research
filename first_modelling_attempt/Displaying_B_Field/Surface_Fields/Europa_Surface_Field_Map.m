function [] = Europa_Surface_Field_Map(B_prim, case_id, step_num, time_step, start)

% This function file is for generating 3d plots of the magnetic field
% magnitudes around the surface of Europa. 
% 
% It will need the primary B fields from Jupiter from the entirety of the 
% simulation (can just pass the B field for non-Zimmer cases), the id of 
% the case the simulation is being done in, the timestep at which the plot 
% will be generated and the timestep length being used, and the starting 
% time of the simulation.

% constants
mu0 = 4 * pi * 1e-7;        % vacuum permeability (H/m)
R_E = 1560e3;               % Europa radius (m)

A = 1;                      % amplitude
phi = 0;                    % phase angle
B_ext = [0; 0; -385e-9];        % external field (T)

% Define magnetic dipole moment
if case_id == "Zimmer"
    M0 = -(4*pi/mu0) * A * exp(1i*phi) * B_prim(:,step_num) .* (R_E^3)/2;
    M_real = real(M0);
elseif case_id == "Case 0" || case_id == "Case 1"
    M_real = B .* 0;        % No induced field for these cases
end

% Function for r_hat and induced field
r_hat = @(r) r / norm(r);
B_sec = @(r, M) (mu0/(4*pi)) * (3 * dot(r_hat(r), M) * r_hat(r) - M) / norm(r)^3;

% Create sphere mesh
n = 100; % resolution
[xs, ys, zs] = sphere(n);
xs = xs * R_E;
ys = ys * R_E;
zs = zs * R_E;

% Initialize magnetic field magnitude matrix
B_magnitude = zeros(size(xs));

% Loop over each surface point to compute total B-field magnitude
for i = 1:size(xs,1)
    for j = 1:size(xs,2)
        r = [xs(i,j); ys(i,j); zs(i,j)];
        B_total = B_prim(:,step_num) + B_ext;             % Primary field
        if case_id == "Zimmer"
            B_total = B_total + B_sec(r, M_real);  % Add induced field
        end
        B_magnitude(i,j) = norm(B_total); % Store field strength
    end
end

% Plotting
figure;
surf(xs, ys, zs, B_magnitude, ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.9);
axis equal;
xlabel('x'); ylabel('y'); zlabel('z');
% title(['Magnetic Field Strength on Europa Surface: ', case_id]);
title(sprintf('Europa Surface Magnetic Field Magnitudes at t = %.2f hours', ...
    (start + step_num * time_step)/3600));

colormap jet;
cb = colorbar;
cb.Label.String = 'B (nT)'; 
clim([min(B_magnitude(:)), max(B_magnitude(:))]);

end