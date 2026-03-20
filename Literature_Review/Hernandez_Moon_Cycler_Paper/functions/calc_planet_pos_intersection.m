function [r_b_1, r_b_2, r_sc_1, r_sc_2, tof1, tof2] = ...
    calc_planet_pos_intersection(r0_vec, a, e, mu, r_b, w0, t0, retrograde)

n_b  = mean_motion(r_b, mu);
n_sc = sqrt(mu/a^3);

% Initial SC angles
phi0   = atan2(r0_vec(2), r0_vec(1));   % inertial
theta0 = mod(phi0 - w0, 2*pi);          % true anomaly

% Solve for intersections r(theta)=r_b
cos_th = (a*(1 - e^2)/r_b - 1)/e;
cos_th = max(-1, min(1, cos_th));

theta1 = acos(cos_th);
theta2 = mod(2*pi - theta1, 2*pi);

% Inertial angles for SC intersections
phi1 = mod(w0 + theta1, 2*pi);
phi2 = mod(w0 + theta2, 2*pi);

% Time from theta0 -> theta_int (forward only)
tof1 = tof_between_thetas(theta0, theta1, e, n_sc, retrograde);
tof2 = tof_between_thetas(theta0, theta2, e, n_sc, retrograde);

% Absolute times (moons start at phi=0 when t=0)
t_abs1 = t0 + tof1;
t_abs2 = t0 + tof2;

% Positions at intersection
r_sc_1 = r_b * [cos(phi1); sin(phi1); 0];
r_sc_2 = r_b * [cos(phi2); sin(phi2); 0];

r_b_1  = r_b * [cos(n_b*t_abs1); sin(n_b*t_abs1); 0];
r_b_2  = r_b * [cos(n_b*t_abs2); sin(n_b*t_abs2); 0];

end