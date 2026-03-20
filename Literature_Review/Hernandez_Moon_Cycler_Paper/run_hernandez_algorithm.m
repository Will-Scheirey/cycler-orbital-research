clear; clc; close all;

jupiter_moons = get_jupiter_moon_params();
r_jupiter     = 69911;
mu            = jupiter_moons.mu_primary;

Delta = deg2rad(5.2);
T_syn = day2sec(7.05);

n_io = (8*pi + Delta) / T_syn;
n_eu = (4*pi + Delta) / T_syn;
n_ga = (2*pi + Delta) / T_syn;

a_io = (mu / n_io^2)^(1/3);
a_eu = (mu / n_eu^2)^(1/3);
a_ga = (mu / n_ga^2)^(1/3);

r_all = [a_io, a_eu, a_ga];
num_bodies = length(r_all);

rp_min = r_jupiter; rp_max = a_io;
ra_min = a_ga;      ra_max = inf;

a_min = (rp_min + ra_min) / 2;
T_min = orbital_period(a_min, mu);

search_space = generate_search_space(T_syn, 5, T_min);

num_e   = 50;
num_phi = 50;
phi_all = linspace(0, 4*pi, num_phi);

n_syn = 1;
r0 = a_eu;
n0 = n_eu;

n_rev_max = search_space(n_syn);

flybys_all_pro = cell(n_rev_max, num_e, num_phi);
flybys_all_ret = cell(n_rev_max, num_e, num_phi);

costs_all_pro  = cell(n_rev_max, num_e, num_phi, 2);
costs_all_ret  = cell(n_rev_max, num_e, num_phi, 2);

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

            % Guess:         [1          x 2]
            %         One for each direction (inbound/outbound)
            
            % Direction:     [1          x 4]

            %         Flybys, r0_vec, v0_vec, w0
            % Flybys:        [num_bodies x 2]
            %         Two intersections for each body

            % Intersections: [1          x 4]
            %         Body pos, sc pos, distance, tof | of/at intersection

            t0   = phi_i / n0;            % epoch offset: moons at anomaly 0 when phi0=0
            
            guess_flybys_pro = generate_guess_flybys(ai, e_i, phi_i, r0, mu, r_all, t0, false);
            flybys_all_pro{n_rev, e_idx, phi_idx} = guess_flybys_pro;
            costs_all_pro{n_rev, e_idx, phi_idx, 1} = cost_fcn(guess_flybys_pro, r_all, 1);
            costs_all_pro{n_rev, e_idx, phi_idx, 2} = cost_fcn(guess_flybys_pro, r_all, 2);

            guess_flybys_ret = generate_guess_flybys(ai, e_i, phi_i, r0, mu, r_all, t0, true);
            flybys_all_ret{n_rev, e_idx, phi_idx} = guess_flybys_ret;
            costs_all_ret{n_rev, e_idx, phi_idx, 1} = cost_fcn(guess_flybys_ret, r_all, 1);
            costs_all_ret{n_rev, e_idx, phi_idx, 2} = cost_fcn(guess_flybys_ret, r_all, 2);
        end
    end
end

%% PLOT

costs_pro = vertcat(costs_all_pro{:}) / r_jupiter;
costs_ret = vertcat(costs_all_ret{:}) / r_jupiter;

max_bin  = max(max(costs_pro), max(costs_ret));
num_bins = 100;
bins = linspace(0, max_bin, num_bins);

figure(1)
clf
histogram(costs_pro, bins, 'Normalization','percentage', 'LineStyle', 'none')
title("Prograde")
ylim([0, 300/num_bins])

figure(2)
clf
histogram(costs_ret, bins, 'Normalization','percentage', 'LineStyle', 'none')
title("Retrograde")
ylim([0, 300/num_bins])

figure(3)
clf
C = cell2mat(costs_all_ret);
[globalMin, linearIdx] = min(C, [], 'all');
[i, j, k, w] = ind2sub(size(C), linearIdx); % For a 3D array

guess  = flybys_all_ret{i, j, k};

phi_i = phi_all(k);
plot_guess(guess, a_io, a_eu, a_ga, mu, n0, r_all, phi_i, w)

function cost = cost_fcn(guess, r_all, direction)
    num_bodies = length(r_all);
    
    direction1   = guess{direction};
    guess_flybys = direction1{1};
    
    closest   = zeros(num_bodies, 1);
    for n = 1:num_bodies
        dist1 = guess_flybys{n, 1}{3};
        dist2 = guess_flybys{n, 2}{3};

        closest(n) = min(dist1, dist2);
    end

    cost = norm(closest, 2);
end

function plot_guess(guess, a_io, a_eu, a_ga, mu, n0, r_all, phi, w)
if nargin < 9
    w = 1;
end
num_bodies = length(r_all);
rev1   = guess{w};

guess_flybys = rev1{1};
guess_orbit  = rev1(2:end);

r0 = guess_orbit{1};
v0 = guess_orbit{2};

clf
hold on

plot_orbit_circular(a_io, mu, 'DisplayName', 'Io', 'LineWidth', 2);
plot_orbit_circular(a_eu, mu, 'DisplayName', 'Europa', 'LineWidth', 2);
plot_orbit_circular(a_ga, mu, 'DisplayName', 'Ganymede', 'LineWidth', 2);

plot_orbit(r0, v0, mu, 'DisplayName', 'Orbit', 'LineWidth', 1.5);
plot(r0(1), r0(2), 'kx', 'DisplayName', 'Initial Position', 'MarkerSize', 20, 'LineWidth', 2)
quiver(r0(1), r0(2), v0(1), v0(2), 3e4, 'k-', 'DisplayName', 'Initial Velocity', 'MarkerSize', 20, 'LineWidth', 2)

distances = zeros(num_bodies, 2);
for n = 1:num_bodies
    distances(n, 1) = guess_flybys{n, 1}{3};
    distances(n, 2) = guess_flybys{n, 2}{3};

    rb1 = guess_flybys{n, 1}{1};
    rb2 = guess_flybys{n, 2}{1};

    rsc1 = guess_flybys{n, 1}{2};
    rsc2 = guess_flybys{n, 2}{2};

    p = plot(rb1(1), rb1(2), 'rx', 'DisplayName', "Intersect 1", 'MarkerSize', 20, 'LineWidth', 2);
    if n ~= 1, p.HandleVisibility = 'off'; end

    p = plot(rb2(1), rb2(2), 'bx', 'DisplayName', "Intersect 2", 'MarkerSize', 20, 'LineWidth', 2);
    if n ~= 1, p.HandleVisibility = 'off'; end

    plot(rsc1(1), rsc1(2), 'rx', 'DisplayName', sprintf("SC at Body %d Intersect 1", n), 'MarkerSize', 20, 'LineWidth', 2, 'HandleVisibility', 'off')
    plot(rsc2(1), rsc2(2), 'bx', 'DisplayName', sprintf("SC at Body %d Intersect 2", n), 'MarkerSize', 20, 'LineWidth', 2, 'HandleVisibility', 'off')

    t0         = phi / n0;
    theta_body = mean_motion(r_all(n), mu) * t0;
    r0 = r_all(n) * [cos(theta_body); sin(theta_body)];

    p = plot(r0(1), r0(2), 'k.', 'DisplayName', 'Starting Position', 'MarkerSize', 20);
    if n ~= 1, p.HandleVisibility = 'off'; end
end

axis equal
legend
end
