    % testing out readings from spacecraft trajectories!
    % To change settings of the flyby, dist is the distance from the
    % surface of Europa at closest approach, and theta is the longitude of 
    % the base of the plume at Europa's surface.

R_E = 1560e3;                                % Europa radius (m)
omega = 2*pi/(11.23*3600);                   % synodic frequency (rad/s)
period_time = 0 * (2 * pi / omega)/4;

num_spots = 10000;
spacecraft_time = linspace(0, (2 * pi / omega)/1000/3600, num_spots);
trajectory = zeros(num_spots, 3);
trajectory(:,1) = linspace(-5*R_E, 5*R_E, num_spots);
% dist = -25e3 - R_E;
% dist = -1.05*R_E;
% dist = -1.25*R_E;
% dist = -1.5*R_E;
dist = -2.5*R_E;
trajectory(:,2) = dist;

B_fields = zeros(num_spots, 3);
neutral_densities = zeros(num_spots, 1);
column_densities = zeros(num_spots, 1);

plume_phis = 90:1:270;
plume_theta = 10;

top_m = 1e16; % for plotting later

filename = sprintf('SpaceCraft Measurements Theta = %d y = %.2f R_E.gif', plume_theta, dist./R_E);

for phi_idx = 1:length(plume_phis)
    clf;

    B_fields = zeros(num_spots, 3);
    neutral_densities = zeros(num_spots, 1);
    column_densities = zeros(num_spots, 1);

    for i = 1:num_spots
        B_fields(i,:) = get_Europa_B_field(trajectory(i,:), period_time);
        neutral_densities(i) = get_neutral_vol_density(trajectory(i,1), trajectory(i,2), ...
            trajectory(i,3), plume_theta, plume_phis(phi_idx), 15, "m");
        column_densities(i) = get_electron_col_density(trajectory(i,:), plume_theta, ...
            plume_phis(phi_idx), 15);
    end
    
    % left axis (B-fields) 
    ax1 = axes;
    hold(ax1, 'on')
    p1 = plot(ax1, spacecraft_time, B_fields(:,1)*1e9, '--', 'DisplayName', 'B_x');
    p2 = plot(ax1, spacecraft_time, B_fields(:,2)*1e9, '--', 'DisplayName', 'B_y');
    p3 = plot(ax1, spacecraft_time, B_fields(:,3)*1e9, '--', 'DisplayName', 'B_z');
    p4 = plot(ax1, spacecraft_time, vecnorm(B_fields,2,2)*1e9, ':', 'DisplayName', 'Total B');
    ylabel(ax1, 'B (nT)')
    xlabel(ax1, 'Time (hr)')
    ax1.YColor = 'k';  
    ax1.Box = 'off';   % remove right y-axis ticks from ax1
    
    % right axis (densities, log scale)
    ax2 = axes;
    hold(ax2, 'on')
    p5 = plot(ax2, spacecraft_time, neutral_densities, '-', 'DisplayName', 'Neutral Densities (m^{-3})');
    p6 = plot(ax2, spacecraft_time, column_densities, '-',  'DisplayName', ...
        'e^- Column Densities (m^{-2})');
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
    t = title(ax1, sprintf('Europa Flyby: Plume Phi = %d°, Plume Theta = %d°, Closest Approach %.0f km', ...
        plume_phis(phi_idx), plume_theta, (abs(dist)-R_E)./1000));
    currentPos = get(t, 'Position');
    set(t, 'Position', [currentPos(1) - 0.0007, currentPos(2) + 10, currentPos(3)]);
    ylim(ax2, [5e5 top_m])
    xlim(ax1, [min(spacecraft_time) max(spacecraft_time)])
    hold(ax1, 'off')
    hold(ax2, 'off')
    drawnow;

    % capture frame and write to gif
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    if phi_idx == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.15);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.15);
    end    
end