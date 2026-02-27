% Based on: https://doi.org/10.2514/3.25519

clear; clc; close all;
mu_sun = 132712000000;                    % [km^3 s^-2]

r_earth_sun = 149.6e6;                    % [km]
r_mars_sun = 227.9e6;                     % [km]

T_earth = 2*pi * sqrt(r_earth_sun^3 / mu_sun);
T_mars = 2*pi * sqrt(r_mars_sun^3 / mu_sun);

T_ratio = T_mars / T_earth;

T_transfer = T_ratio * 3600*24*365;

synodic_period = 1 / (1/T_earth - 1/T_mars);

angle_turn = mod(T_ratio * 360, 360);
if angle_turn > 180
    angle_turn = angle_turn - 360;
end

[a_transfer, e_transfer] = solve_transfer_geometry(r_earth_sun, T_transfer, mu_sun, 0);

eps_transfer = -mu_sun/(2*a_transfer);
h_transfer = sqrt(-1/(2*eps_transfer)*mu_sun^2*(1-e_transfer^2));

v0_transfer = mu_sun/h_transfer *[
    e_transfer*sin(0)
    1 + e_transfer*cos(0)
];

period = T_earth * (ceil(T_ratio));

num_steps = 300;
ts = linspace(0, synodic_period, num_steps);

thetas = linspace(0, 2*pi, 100);
r_transfer = a_transfer * (1 - e_transfer^2) ./ (1 + e_transfer * cos(thetas));

figure(1)

plot(r_earth_sun * cos(thetas), r_earth_sun * sin(thetas), '-b', 'LineWidth', 1.5); hold on
plot(r_mars_sun * cos(thetas), r_mars_sun * sin(thetas), '-r', 'LineWidth', 1.5);

plot(r_transfer .* cos(thetas), r_transfer .* sin(thetas), '-g', 'LineWidth', 1.5);

legend("Earth", "Mars", "Spacecraft")

title("Aldrin Cycler")

axis equal

function [a, e] = solve_transfer_geometry(r0, T, mu, theta0)
       
    a = ((T/(2*pi))^2 * mu)^(1/3);

    syms e_sym;
    eqn = r0 == a * (1 - e_sym^2) / (1 + e_sym * cosd(theta0));
    e_sol = double(solve(eqn, e_sym));
    
    e = e_sol(find(e_sol > 0, 1));
end