% testing plotting of B field

% constants
B0 = 1;             % base magnetic field strength
RE = 6.371e6;       % Earth's radius in meters
N = 10;             % Number of grid points per axis
range = 2;          % domain range in units of RE (from -2RE to 2RE)

% 3D grid in units of RE
[x, y, z] = meshgrid(linspace(-range, range, N));
r = sqrt(x.^2 + y.^2 + z.^2);

% avoid dividing by zero at origin
r(r == 0) = NaN;

% components of the magnetic field
Bx = -3*B0 .* x.*z ./ ( (r/RE).^3 .* r.^2 );
By = -3*B0 .* y.*z ./ ( (r/RE).^3 .* r.^2 );
Bz =  B0 .* (x.^2 + y.^2 - 2*z.^2) ./ ( (r/RE).^3 .* r.^2 );

% normalize the vectors for better visualization
L = sqrt(Bx.^2 + By.^2 + Bz.^2);
Bx = Bx ./ L;
By = By ./ L;
Bz = Bz ./ L;

figure;
quiver3(x, y, z, Bx, By, Bz, 0.5, 'b');  % 0.5 scales arrow size
axis equal;
xlabel('x (RE)');
ylabel('y (RE)');
zlabel('z (RE)');
title('Magnetic Field of a Dipole in 3D');
grid on;
saveas(gcf,'Planetary_Dipole_BField_Model_Attempt','jpg');