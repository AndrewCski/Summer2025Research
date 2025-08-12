% testing out a spacecraft trajectory!

R_E = 1560e3;                                % Europa radius (m)
omega = 2*pi/(11.23*3600);                   % synodic frequency (rad/s)
period_time = 0 * (2 * pi / omega)/4;

num_spots = 10000;
spacecraft_time = linspace(0, (2 * pi / omega)/100000/3600, num_spots);
trajectory = zeros(num_spots, 3);
trajectory(:,2) = linspace(-5*R_E, 5*R_E, num_spots);
trajectory(:,1) = -1.05*R_E;

B_fields = zeros(num_spots, 3);
neutral_densities = zeros(num_spots, 1);
column_densities = zeros(num_spots, 1);

for i = 1:num_spots
    B_fields(i,:) = get_Europa_B_field(trajectory(i,:), period_time);
    neutral_densities(i) = get_neutral_vol_density(trajectory(i,1), trajectory(i,2), ...
        trajectory(i,3), 0, 0, 15, "m");
    column_densities(i) = get_neutral_col_density(trajectory(i,:), 0, 0, 15);
end

%%

figure;

% left axis (B-fields) 
ax1 = axes;
hold(ax1, 'on')
p1 = plot(ax1, spacecraft_time, B_fields(:,1)*1e9, '--', 'DisplayName', 'B_x');
p2 = plot(ax1, spacecraft_time, B_fields(:,2)*1e9, '--', 'DisplayName', 'B_y');
p3 = plot(ax1, spacecraft_time, B_fields(:,3)*1e9, '--', 'DisplayName', 'B_z');
p4 = plot(ax1, spacecraft_time, vecnorm(B_fields,2,2)*1e9, ':', 'DisplayName', 'Total B');
ylabel(ax1, 'B_E (nT)')
xlabel(ax1, 'Time (hr)')
ax1.YColor = 'k';  
ax1.Box = 'off';   % remove right y-axis ticks from ax1

% right axis (densities, log scale)
ax2 = axes;
hold(ax2, 'on')
p5 = plot(ax2, spacecraft_time, neutral_densities, '-', 'DisplayName', 'Neutral Densities (m^{-3})');
p6 = plot(ax2, spacecraft_time, column_densities.*1e-4, '-',  'DisplayName', 'Column Densities (cm^{-2})');
set(ax2, 'YScale', 'log')
ylabel(ax2, 'Densities')
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax2.XColor = 'k';
ax2.YColor = 'k';

% match positions & link x-axes
ax2.Position = ax1.Position;     % keep exact alignment
linkaxes([ax1 ax2], 'x')         % sync horizontal zoom/pan

legend([p1 p2 p3 p4 p5 p6], 'Location', 'best')
title(ax1, 'Europa Flyby')