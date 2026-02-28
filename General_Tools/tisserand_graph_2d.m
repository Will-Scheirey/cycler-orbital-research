clear; clc; close all

AU = 149.6e6;
r_jupiter = 69911;

% planets = get_simplified_planet_params();

bodies = get_jupiter_moon_params();
mu_primary = bodies.mu_primary;

num_v_inf = 10;
v_infs = linspace(1, 10, num_v_inf);

num_alphas = 50;
alphas = linspace(0, 180, num_alphas);

name_b = fieldnames(bodies);
num_bodies = length(name_b);

colors = turbo(num_bodies);

figure(1)
hold on

for idx = 1:length(name_b)
    name = name_b{idx};
    if strcmp(name, 'mu_primary')
        continue;
    end

    body = bodies.(name);

    points = generate_tisserand_points(v_infs, alphas, mu_primary, body.r_primary);
    plot_tisserand(points, body.mu, body.r, v_infs, alphas, colors(idx, :), name_b{idx}, r_jupiter);
end

legend
xlabel("Periapsis Radius (RJ)")
ylabel("Apoapsis Radius (RJ")
title("Tisserand Plots for Multiple Bodies")

xlim([1, 15]);
ylim([5, 25]);

function plot_tisserand(points, mu, rp_min, v_infs, alphas, color, name, reference_length)
line_width = 1;
marker_size = 8;

num_points = height(points);

for i = 1:num_points
    delta_max = delta_from_rp(mu, v_infs(i), rp_min);
    alpha_step = find_closest_idx(alphas, delta_max);
    curr_step = 0;

    arg1 = squeeze(points(i, :, 1)) / reference_length; % rp
    arg2 = squeeze(points(i, :, 2)); % Specific Energy
    arg3 = squeeze(points(i, :, 3)); % Period
    arg4 = squeeze(points(i, :, 4)); % Inclination
    arg5 = squeeze(points(i, :, 5)) / reference_length; % ra

    x_point = arg1;
    y_point = arg5;

    if i == 1
        plot(x_point, y_point, 'Color', color, 'LineStyle', '-', 'DisplayName', name, 'LineWidth', line_width);
    else
        plot(x_point, y_point, 'Color', color, 'LineStyle', '-','HandleVisibility', 'off', 'LineWidth', line_width);
    end

    if floor(alpha_step) > curr_step
        plot(x_point, y_point, 'Color', color, 'Marker', '.', 'MarkerSize', marker_size, 'HandleVisibility', 'off');
    end

    mid  = floor(length(x_point) / 2);
    text(x_point(mid) + 0.1, y_point(mid) - 0.1, sprintf("%d", i))
end

end

function delta = delta_from_rp(mu, v_inf, rp)
e = 1 + (rp*v_inf^2)/mu;
delta = 2*asind(1/e);
end

function rp = rp_from_delta(mu, v_inf, delta_deg)
e = 1/sind(delta_deg/2);
rp = (e - 1) * (mu / v_inf^2);
end

function points = generate_tisserand_points(v_infs, alphas, mu_2, r_1_2)
num_v_inf = length(v_infs);
num_alphas = length(alphas);

v1_2 = sqrt(mu_2 / r_1_2);

v_vec_1_2 = [0, v1_2, 0];

points = nan(num_v_inf, num_alphas, 5);

for i=1:num_v_inf
    v_inf_1 = v_infs(i);
    for j=1:num_alphas
        alpha = alphas(j);

        r_vec_2 = [r_1_2, 0, 0];
        v_inf_vec_1 = v_inf_1 * [sind(alpha), cosd(alpha), 0];

        v_vec_2 = v_vec_1_2 + v_inf_vec_1;
        v_2 = norm(v_vec_2);

        h_vec_2 = cross(r_vec_2, v_vec_2);
        h_2 = norm(h_vec_2);

        epsilon_2 = v_2^2/2 - mu_2 / r_1_2;
        a_2 = -mu_2 / (2*epsilon_2);

        e_2 = sqrt(1 - h_2^2/mu_2 * 1/a_2);

        if e_2 > 1
            continue
        end

        rp_2 = h_2^2/mu_2 * 1/(1 + e_2);
        ra_2 = h_2^2/mu_2 * 1/(1 - e_2);

        period = 2*pi*sqrt(a_2^3/mu_2);

        if imag(period) ~= 0 || period > 1e9
            period = 0;
        end

        inclination = acosd(h_vec_2(3) / h_2);

        points(i, j, :) = [rp_2, epsilon_2, period / 86400 / 365, inclination, ra_2];
    end
end
end

function idx = find_closest_idx(arr, target)
[~, idx] = min(abs(arr - target));
end