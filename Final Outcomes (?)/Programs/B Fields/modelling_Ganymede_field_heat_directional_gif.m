
% This file is used for creating gifs of the evolution of the magnetic
% field or electric field around Ganymede during one synodic period. 
%
% To change the axes/B field component the gif, alter the lines labelled with #1, #2, etc.

omega = 2*pi/(10.53*3600);                % synodic frequency (rad/s)

% simulation parameters
numsteps = 10000;               % number of time steps
period = 2 * pi / omega;        % modelling a synodic period
time = linspace(0, period, numsteps);

B0_vec = [15 * 10 ^-9; 80 * 10 ^-9; 0]; % 15 nT azimuthal, 80 nT radial, using EPhiO
% (-75 nT constant added to North-South, represented by z, in upcoming loop)
B_prim = zeros(3, numsteps);
for u = 1:length(time)
    B_prim(:,u) = [B0_vec(1) .* cos(omega .* time(u)); B0_vec(2) .* -sin(omega .* time(u) - pi/4.5); 0];
end

mu0 = 4*pi*1e-7;  % vacuum permeability
R_G = 2631e3;     % Ganymede radius (m)

% induced dipole moment
M_ind = ((-2 * pi * 0.84 * R_G^3) / mu0) .* B_prim;
M_dip = [-4.1; 9; -131] * 1e18;
M_real = M_ind + M_dip;
B_prim(3,:) = -75 * 10^-9;

% function to compute dipole magnetic field at a point
B_sec = @(r, M) (mu0/(4*pi)) * (3*dot(r, M)*r - M*norm(r)^2) / norm(r)^5;

% define grid 
sim_range = linspace(-5*R_G, 5*R_G, 1000);  % adjust range and resolution

% CHANGE DEPENDING ON AXES USED...
[Y, Z] = meshgrid(sim_range, sim_range);                           % #1 !!!
X = zeros(size(Z));                                                % #2 !!!

r_mag = sqrt(X.^2 + Y.^2 + Z.^2);

% define exclusion zone
r_min = R_G;  % exclude inside this radius
mask = (r_mag >= r_min);

filename = 'Final Ganymede Bz on yz evolution.gif';        % #3 !!!

% loop over time indices
for t_idx = 1:50:numsteps   % take every 50th step for speed

    % B_prim and M_real at this time
    Bprim_now = B_prim(:, t_idx);
    M = M_real(:, t_idx);

    % compute field component                     % all of these are #4 !!!
    B_x = NaN(size(X));  % NaN -> will not plot
    B_y = NaN(size(X));
    B_z = NaN(size(X));
    % B_tot = NaN(size(X));
    % E_x = NaN(size(X)); 
    % E_y = NaN(size(X)); 
    % E_z = NaN(size(X)); 
    % E_tot = NaN(size(X)); 

    for i = 1:numel(X)
        if ~mask(i)
            continue;
        end

        r = [X(i); Y(i); Z(i)];
        Bsec = B_sec(r, M);
        Btotal = Bprim_now + Bsec;                                 % #5 !!!
        % Etotal = -cross([90000; 0; 0], Btotal); 
                                                  % all of these are #6 !!!
        B_x(i) = Btotal(1);  % x = 1, y = 2, z = 3, norm is total
        B_y(i) = Btotal(2); 
        B_z(i) = Btotal(3); 
        % B_tot(i) = norm(Btotal); 
    end

    % plot this particular instance
    imagesc(sim_range/1e3, sim_range/1e3, B_z * 1e9);                  % #7 !!!
    set(gca, 'YDir', 'normal');
    colorbar;
    clim([-1500 650]);   % fixing color scale, depends on what field is modeled       % #8 !!!
    xlabel('y (km)');                                              % #9 !!!
    ylabel('z (km)');                                             % #10 !!!
    title(sprintf('Ganymede Magnetic Environment: t = %.2f hrs', time(t_idx)/3600));
    cb = colorbar;
    cb.Label.String = 'B_z (nT)';                                 % #11 !!!
    cb.Label.Rotation = 0; % to rotate the text
    hold on
    str_slc = streamslice(Y/1e3, Z/1e3, B_y, B_z, 1.5);                % #12 !!!
    set(str_slc,'Color','r');
    hold off
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