% program for creating a gif of Ganymede's B field at its surface, plotted
% respect to latitude and longitude (0 degrees is towards Jupiter, 90 degrees
% is towards Ganymede orbital direction, etc.). To change modelled field,
% change the B_x/whatever is on parts

omega = 2*pi/(10.53*3600);                            % synodic frequency (rad/s)
R_G = 2631e3;                                         % Ganymede radius (m)
mu0 = 4*pi*1e-7;                                      % vacuum permeability (H/m)

measure_gap = 1;

thetas = (-90 + measure_gap/2):measure_gap:(90 - measure_gap/2);
phis = (measure_gap/2):measure_gap:(360 - measure_gap/2);
plot_thetas = -90:measure_gap:90;
plot_phis = 0:measure_gap:360;

time_idx = 1;           % index of the chosen timestep
numsteps = 10000;       % number of time steps
time = linspace(0, 10.53*3600, numsteps);

B0_vec = [15 * 10 ^-9; 80 * 10 ^-9; 0]; 
M_dip = [-4.1; 9; -131] * 1e18;

B_x = zeros(length(thetas), length(phis));
B_y = zeros(length(thetas), length(phis));
B_z = zeros(length(thetas), length(phis));
B_tot = zeros(length(thetas), length(phis));

figure;

filename = 'Ganymede Surface Bz evolution.gif'; 

for t_idx = 1:49:numsteps   % take every 49th step for speed

    % B_prim and M at this time
    B_prim_now = [B0_vec(1) .* cos(omega .* time(t_idx)); B0_vec(2) .* -sin(omega .* ...
        time(t_idx) - pi/4.5); 0];
    M_ind = ((-2 * pi * 0.84 * R_G^3) / mu0) .* B_prim_now;
    M = M_dip + M_ind;
    B_prim_now(3) = -75 * 10^-9;

    index = 1;
    for i = 1:length(thetas)
        for j = 1:length(phis)
            x = R_G * cos(deg2rad(thetas(i))) * sin(deg2rad(phis(j)));
            y = R_G * cos(deg2rad(thetas(i))) * cos(deg2rad(phis(j)));
            z = R_G * sin(deg2rad(thetas(i)));
            r = [x; y; z];
    
            B_sec = (mu0/(4*pi)) .* (3.*dot(r, M) .* r - ...
                M .* norm(r).^2) ./ norm(r)^5;
            % B_x(i,j)   = B_prim_now(1) + B_sec(1);
            % B_y(i,j)   = B_prim_now(2) + B_sec(2);
            B_z(i,j)   = B_prim_now(3) + B_sec(3);
            % B_tot(i,j) = norm(B_prim_now + B_sec);
  
            index = index + 1;
        end
    end
    
    disp(max(B_z(:) .* 1e9))

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
    clim([-1.5211e+03  651.8049])
    title(sprintf('Ganymede Surface B_z at t = %.2f hours', time(t_idx)/3600));
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