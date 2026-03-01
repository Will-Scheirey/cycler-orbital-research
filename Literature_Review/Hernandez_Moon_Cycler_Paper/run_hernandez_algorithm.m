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

r_all = [a_io, a_eu, a_ga];

rp_min = r_jupiter; rp_max = a_io;
ra_min = a_ga;      ra_max = inf;

a_min = (rp_min + ra_min) / 2;
T_min = orbital_period(a_min, mu);

search_space = generate_search_space(T_syn, 5, T_min);

num_e   = 10;
num_phi = 10;
phi_all = linspace(0, 2*pi, num_phi);

n_syn = 1;
r0 = a_eu;

n_rev_max = search_space(n_syn);

for n_rev = 1:n_rev_max

    Ti = (n_syn * T_syn) / n_rev;
    ai = (mu * (Ti/(2*pi))^2) ^ (1/3);
    
    e_min = max(1 - rp_max/ai, ra_min/ai - 1);
    e_max = min(1 - rp_min/ai, ra_max/ai - 1);

    e_all   = linspace(e_min, e_max, num_e);

    for e_idx = 1:num_e
        e_i = e_all(e_idx);
        for phi_idx = 1:num_phi
            phi_i = phi_all(phi_idx);

            flybys_all = generate_guess_flybys(ai, e_i, phi_i, r0, mu, r_all);
        end
    end
end

figure(1)
clf
hold on

plot_orbit_circular(a_io, mu, 'DisplayName', 'Io', 'LineWidth', 2);
plot_orbit_circular(a_eu, mu, 'DisplayName', 'Europa', 'LineWidth', 2);
plot_orbit_circular(a_ga, mu, 'DisplayName', 'Ganymede', 'LineWidth', 2);

axis equal
legend

function results = generate_guess_flybys(ai, ei, phi0, r0, mu, r_all)
[r0_vec, v0_1, w0_1] = calc_orbit(ai, ei, phi0, r0, mu, true);
[~,  v0_2, w0_2]     = calc_orbit(ai, ei, phi0, r0, mu, false);

flybys_1 = eval_orbit(r0_vec, v0_1, ai, ei, w0_1, mu, r_all);
flybys_2 = eval_orbit(r0_vec, v0_2, ai, ei, w0_2, mu, r_all);

results = {flybys_1, flybys_2};

end

function flybys = eval_orbit(r0, v0, ai, ei, w0, mu, r_all)
num_b = length(r_all);
flybys = cell(num_b, 2);

for n = 1:num_b
    r_b = r_all(n);

    [r_b_1, r_b_2, r_sc_1, r_sc_2, tof1, tof2] = calc_planet_pos_intersection(r0, v0, ai, ei, mu, r_b, w0);

    dist1 = norm(r_b_1 - r_sc_1);
    dist2 = norm(r_b_2 - r_sc_2);

    flybys{n, 1} = {r_b_1, r_sc_1, dist1, tof1};
    flybys{n, 2} = {r_b_2, r_sc_2, dist2, tof2};
end
end

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

function [r_b_1, r_b_2, r_sc_1, r_sc_2, tof_int1, tof_int2] = calc_planet_pos_intersection(r0, v0, a, e, mu, r_b, w0)
n = mean_motion(r_b, mu);

theta_sc_int1 = get_intersection(r0, v0, mu, r_b) - w0;
tof_int1      = calc_tof(a, mu, e, theta_sc_int1);

theta_sc_int2 = 2*pi - theta_sc_int1;
tof_int2      = calc_tof(a, mu, e, theta_sc_int2);

r_b_1   = r_b * [cos(n*tof_int1);    sin(n*tof_int1);    0];
r_sc_1  = r_b * [cos(theta_sc_int1); sin(theta_sc_int1); 0];

r_b_2   = r_b * [cos(n*tof_int2); sin(n*tof_int2); 0];
r_sc_2  = r_b * [cos(theta_sc_int2); sin(theta_sc_int2); 0];
end

function out = generate_search_space(T_syn, syn_per_max, T_min)
out = zeros(syn_per_max, 1);

for syn_per = 1:syn_per_max
    out(syn_per) = floor(T_syn * syn_per / T_min);
end
end

function wrap2pi(ang)
    ang = mod(ang, 2*pi);
end