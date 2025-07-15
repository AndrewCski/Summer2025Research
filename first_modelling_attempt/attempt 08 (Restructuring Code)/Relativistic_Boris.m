function [r_new, v_new] = Relativistic_Boris(r,v,timestep,B_prim, E, m, q)

% Relativistic Boris Pusher Step. 
% 
% Takes current position, velocity, primary B-field, and E-field, alongside
% chosen timestep, particle mass, and particle charge as parameters. 
% Outputs new position and velocity of the particle.

c = 299792458;                          % speed of light (m/s)
A = 1;
phi = 0;

mu0 = 4*pi*1e-7;                        % vacuum permeability (H/m)
R_E = 1560e3;                           % Europa radius (m)

M0 = -(4*pi/mu0) * A * exp(1i*phi) * B_prim .* (R_E^3)/2;
M_real = real(M0);

r_hat = @(r) r / norm(r);
B_sec = @(r, M) (mu0/(4*pi)) * (3*dot(r_hat(r), M)*r_hat(r) - M) / norm(r)^3;


% magnetic field at this position and time
B = B_sec(r, M_real(:)) + B_prim;

% current velocity and momentum
p = m * v / sqrt(1 - norm(v)^2 / c^2);  % relativistic momentum

% step 1: half electric field push
p_minus = p + (q * timestep / 2) * E;

% step 2: magnetic rotation (Boris rotation in momentum space)
gamma_minus = sqrt(1 + norm(p_minus)^2 / (m^2 * c^2));
t = (q * timestep / (2 * m * gamma_minus)) * B;
t_mag2 = dot(t, t);
s = 2 * t / (1 + t_mag2);

% p' = p_minus + p_minus x t
p_prime = p_minus + cross(p_minus, t);
% p_plus = p_minus + (p_prime x s)
p_plus = p_minus + cross(p_prime, s);

% step 3: second half electric field push
p_new = p_plus + (q * timestep / 2) * E;

% final momentum and velocity update
gamma_new = sqrt(1 + norm(p_new)^2 / (m^2 * c^2));
v_new = p_new / (gamma_new * m);

% update position using updated velocity
r_new = r + v_new * timestep;
