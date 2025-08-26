% program for creating a gif of Europa's B field at its surface, plotted
% respect to latitude and longitude (0 degrees is towards Jupiter, 90 degrees
% is towards Europa orbital direction, etc.). To change modelled field,
% change the B_x/whatever is on parts

omega = 2*pi/(11.23*3600);                            % synodic frequency (rad/s)
R_E = 1560e3;                                         % Europa radius (m)
mu0 = 4*pi*1e-7;           % vacuum permeability (H/m)

measure_gap = 1;

thetas = (-90 + measure_gap/2):measure_gap:(90 - measure_gap/2);
phis = (measure_gap/2):measure_gap:(360 - measure_gap/2);
plot_thetas = -90:measure_gap:90;
plot_phis = 0:measure_gap:360;

time_idx = 1;   % index of the chosen timestep
numsteps = 10000;               % Number of time steps
time = linspace(0, 11.23*3600, numsteps);

B0_vec = [75 * 10 ^-9; 200 * 10 ^-9; 0]; % 75 nT azimuthal, 200 nT radial, using Zimmer coords

B_x = zeros(length(thetas), length(phis));
B_y = zeros(length(thetas), length(phis));
B_z = zeros(length(thetas), length(phis));
B_tot = zeros(length(thetas), length(phis));

figure;

filename = 'Europa Surface Btot evolution.gif'; 

for t_idx = 1:49:numsteps   % take every 49th step for speed

    % B_prim and M at this time
    B_prim_now = [B0_vec(1) .* cos(omega .* time(t_idx)); B0_vec(2) .* -sin(omega .* time(t_idx)); 0];
    M = -(4*pi/mu0) .* B_prim_now .* (R_E^3)./2;
    B_prim_now(3) = -385e-9;

    index = 1;
    for i = 1:length(thetas)
        for j = 1:length(phis)
            x = R_E * cos(deg2rad(thetas(i))) * sin(deg2rad(phis(j)));
            y = R_E * cos(deg2rad(thetas(i))) * cos(deg2rad(phis(j)));
            z = R_E * sin(deg2rad(thetas(i)));
            r = [x; y; z];
    
            B_sec = (mu0/(4*pi)) .* (3.*dot(r, M) .* r - ...
                M .* norm(r).^2) ./ norm(r)^5;
            % B_x(i,j)   = B_prim_now(1) + B_sec(1);
            % B_y(i,j)   = B_prim_now(2) + B_sec(2);
            % B_z(i,j)   = B_prim_now(3) + B_sec(3);
            B_tot(i,j) = norm(B_prim_now + B_sec);
  
            index = index + 1;
        end
    end
    
    disp(max(B_tot(:) .* 1e9))

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
    clim([333 564])
    title(sprintf('Europa Surface B_{tot} at t = %.2f hours', time(t_idx)/3600));
    drawnow;

    % capture frame and write to gif
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    if t_idx == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.075);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.075);
    end
end