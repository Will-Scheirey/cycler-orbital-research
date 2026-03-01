clear; clc; close all;

jupiter_moons = get_jupiter_moon_params();
mu            = jupiter_moons.mu_primary;
r_jupiter     = 69911;

delta = deg2rad(5.2);
T_syn = day2sec(7.05); 

a_io = jupiter_moons.io.r_primary;
n_io = mean_motion(a_io, mu);
a_eu = a_io * ((8*pi + delta) / (4*pi + delta))^(2/3);
n_eu = mean_motion(a_eu, mu);
a_ga = a_io * ((8*pi + delta) / (2*pi + delta))^(2/3);
n_ga = mean_motion(a_ga, mu);

rp_min = r_jupiter; rp_max = a_io;
ra_min = a_ga;      ra_max = inf;

a_min = (rp_min + ra_min) / 2;
T_min = orbital_period(a_min, mu);
search_space = generate_search_space(T_syn, 5, T_min);

n_syn = 1;
n_rev = 2;

Ti = (n_syn * T_syn) / n_rev;
ai = (mu * (Ti/2*pi)^2) ^ (1/3);

e_min = max(1 - rp_max/ai, ra_min/ai - 1);
e_max = min(1 - rp_min/ai, ra_max/ai - 1);

ei = e_min;
r0 = a_eu;
phi0 = 4*pi/3;

[r0_vec, v0_vec1, w0_1] = calc_orbit(ai, ei, phi0, r0, mu, true);
[~,      v0_vec2, w0_2] = calc_orbit(ai, ei, phi0, r0, mu, false);

[r_io_1, r_io_2] = calc_planet_pos_intersection(r0_vec, v0_vec1, ai, ei, mu, a_io, w0_2);
[r_eu_1, r_eu_2] = calc_planet_pos_intersection(r0_vec, v0_vec2, ai, ei, mu, a_eu, w0_2);
[r_ga_1, r_ga_2] = calc_planet_pos_intersection(r0_vec, v0_vec2, ai, ei, mu, a_ga, w0_2);

figure(1)
clf

plot_orbit(r0_vec, v0_vec2, mu, 'DisplayName', 'Transfer', 'LineWidth', 2); hold on
plot_orbit_circular(a_io, mu, 'DisplayName', 'Io', 'LineWidth', 2);
plot_orbit_circular(a_eu, mu, 'DisplayName', 'Europa', 'LineWidth', 2);
plot_orbit_circular(a_ga, mu, 'DisplayName', 'Ganymede', 'LineWidth', 2);

plot(r0_vec(1), r0_vec(2), 'rx', 'DisplayName', 'Initial Position', 'MarkerSize', 20, 'LineWidth', 2)

plot(r_io_1(1), r_io_1(2), 'kx', 'HandleVisibility', 'off', 'MarkerSize', 20, 'LineWidth', 2)
plot(r_io_2(1), r_io_2(2), 'kx', 'HandleVisibility', 'off', 'MarkerSize', 20, 'LineWidth', 2)
plot(r_eu_1(1), r_eu_1(2), 'kx', 'HandleVisibility', 'off', 'MarkerSize', 20, 'LineWidth', 2)
plot(r_eu_2(1), r_eu_2(2), 'kx', 'HandleVisibility', 'off', 'MarkerSize', 20, 'LineWidth', 2)
plot(r_ga_1(1), r_ga_1(2), 'kx', 'HandleVisibility', 'off', 'MarkerSize', 20, 'LineWidth', 2)
plot(r_ga_2(1), r_ga_2(2), 'kx', 'HandleVisibility', 'off', 'MarkerSize', 20, 'LineWidth', 2)


axis equal
legend

function [r0_vec, v0_vec, w0] = calc_orbit(ai, ei, phi0, r0, mu, other_branch)
    if nargin < 6
        other_branch = false;
    end

    v0 = sqrt(mu * (2/r0 - 1/ai));
    epsilon = v0^2/2 - mu/r0;
    h = sqrt(-1/2 * (mu^2/epsilon) * (1 - ei^2));
    
    cos_theta = (ai*(1 - ei^2)/r0 - 1)/ei;
    cos_theta = max(-1, min(1, cos_theta));
    
    theta_i = acos(cos_theta);
    if other_branch, theta_i = 2*pi - theta_i; end

    w0 = phi0 - theta_i;
    
    vr = mu/h *      ei*sin(theta_i);
    vt = mu/h * (1 + ei*cos(theta_i));
    
    vx = vr*cos(phi0) - vt*sin(phi0);
    vy = vr*sin(phi0) + vt*cos(phi0);
    
    r0_vec = r0 * [cos(phi0 ); sin(phi0); 0];
    v0_vec = [vx; vy; 0];
end

function [r0_1, r0_2] = calc_planet_pos_intersection(r0, v0, a, e, mu, r_planet, w0)
    n = mean_motion(r_planet, mu);

    theta_int1 = get_intersection(r0, v0, mu, r_planet) - w0;
    tof_int1   = calc_tof(a, mu, e, theta_int1);

    theta_int2 = 2*pi - theta_int1;
    tof_int2   = calc_tof(a, mu, e, theta_int2);


    r0_1   = r_planet * [cos(n*tof_int1); sin(n*tof_int1); 0];
    r0_2   = r_planet * [cos(n*tof_int2); sin(n*tof_int2); 0];
end

function out = generate_search_space(T_syn, syn_per_max, T_min)
    out = zeros(syn_per_max, 1);

    for syn_per = 1:syn_per_max
        out(syn_per) = floor(T_syn * syn_per / T_min);
    end
end