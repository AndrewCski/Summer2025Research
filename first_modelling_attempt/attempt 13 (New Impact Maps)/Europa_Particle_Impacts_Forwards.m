function [ ] = Europa_Particle_Impacts_Forwards(impact_coords, energies)

% This function file is used for mapping the impact sites of various
% energies of particles across the surface of Europa, when tracing
% particles foward in time. 
%
% The two parameters that must be passed are an array consisting of the 
% impact coordinates of the modelled particles, as well as an array of the
% starting energies of the particles.


R_E = 1560e3;                          % Europa radius in meters

figure;
hold on
xlabel('x'); ylabel('y'); zlabel('z')
title('Particle Landing Sites')
[xs, ys, zs] = sphere(50);             % higher resolution sphere
surf(R_E * xs, R_E * ys, R_E * zs, ...
    'FaceColor', [0.2 0.5 1], ...      % easy-ish to see blue
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.5);                 % semi-transparent

% set view and axis
axis equal
grid on
view(3)

interp_points = zeros(3, size(impact_coords, 1));

for i = 1:size(impact_coords, 1)
    p0 = impact_coords(i,1:3);
    p1 = impact_coords(i,4:6);    
    if norm(p0) >= R_E && norm(p1) <= R_E
        d  = p1 - p0;

        a = dot(d,d);
        b = 2*dot(d,p0);
        c = dot(p0,p0) - R_E^2;

        discriminant = b^2 - 4*a*c;

        if discriminant >= 0
            sqrt_disc = sqrt(discriminant);
            f1 = (-b + sqrt_disc) / (2*a);
            f2 = (-b - sqrt_disc) / (2*a);
            f_valid = [f1 f2];
            f_valid = f_valid(f_valid >= 0 & f_valid <= 1);

            if ~isempty(f_valid)
                f = f_valid(1);        % pick the first valid fraction
                interpolated_point = p0 + f*d;
                interp_points(:, i) = interpolated_point;
            end
        end

    end
end

for i = 1:size(impact_coords, 1)
    if norm(interp_points(:, i)) ~= 0
        scatter3(interp_points(1,i), interp_points(2,i), interp_points(3,i), ...
            [], energies(i), 'filled');
    end
end
colorbar
clim([min(energies) max(energies)])
set(gca,'ColorScale','log')
hold off