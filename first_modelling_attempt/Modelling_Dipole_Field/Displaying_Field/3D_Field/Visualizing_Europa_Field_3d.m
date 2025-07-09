% constants
q = 1.602e-19;              
m = 1.67262192e-27;         

omega = 2*pi/(11.23*3600);  

numsteps = 10000;               
period = 2 * pi / omega;        
time = linspace(0, period, numsteps);

B0_vec = [75e-9; 200e-9; 0];  
B_prim = real(B0_vec .* exp(1i * omega * time));
B_prim(3,:) = -350e-9;

A = 1;
phi = 0;

mu0 = 4*pi*1e-7;  
R_E = 1560e3;     

M0 = -(4*pi/mu0) * A * exp(1i*phi) * B_prim .* (R_E^3)/2;
M_real = real(M0);

r_hat = @(r) r / norm(r);
B_sec = @(r, M) (mu0/(4*pi)) * (3*dot(r_hat(r), M)*r_hat(r) - M) / norm(r)^3;

time_idx = 1;   

[x, y, z] = meshgrid(linspace(-3, 3, 15));  % reduced grid for clarity
Bx = zeros(size(x));
By = zeros(size(y));
Bz = zeros(size(z));

Bprim_now = B_prim(:, time_idx);
M = M_real(:, time_idx);

for i = 1:numel(x)
    r_vec = [x(i); y(i); z(i)];  % position in R_E units
    r_norm = norm(r_vec);

    if r_norm < 1  % inside the sphere
        Btotal = [NaN; NaN; NaN];  % mark as NaN so quiver3 skips it
    else
        Bsec = B_sec(r_vec, M);
        Btotal = Bprim_now + Bsec;
    end

    Bx(i) = Btotal(1);
    By(i) = Btotal(2);
    Bz(i) = Btotal(3);
end

figure;

[xe, ye, ze] = sphere(30);
surf(xe, ye, ze, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', [0.5 0.5 0.5]);
hold on

quiver3(x, y, z, Bx, By, Bz, 2, 'b');
axis equal;
xlabel('x (RE)');
ylabel('y (RE)');
zlabel('z (RE)');
title('Magnetic Field of Europa in 3D');
grid on;
saveas(gcf,'Europa_Field_3D','jpg');
hold off