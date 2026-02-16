clear; clc; close all

[consts, mu_altaira, planet_data] = load_problem_data();

for i = 1:height(planet_data)
    plot_planet_orbit(mu_altaira, planet_data(i, :));
end

plot(200*consts.AU, 0, 'gx', 'MarkerSize', 30, 'DisplayName', 'SC', 'LineWidth', 3)
viscircles([0, 0; 0, 0], [0.05*consts.AU; 0.01*consts.AU], 'Color', 'r', 'LineWidth', 3, 'LineStyle', ':');

xlim([-1, 1] * 3.2e10)
ylim([-1, 1] * 3.2e10)

axis square;

legend
title("2D Orbits")

function plot_planet_orbit(mu_altaira, planet_data)
    eps = -mu_altaira/(2*planet_data.a_km);

    % planet_data.e^2

    h = sqrt(-1/2 * (mu_altaira^2/eps) * (1 - planet_data.e^2));

    thetas = linspace(0, 360, 1000) + planet_data.theta0_deg;

    r = h^2 ./ (mu_altaira * (1 + planet_data.e * cosd(thetas)));

    pos = [r .*cosd(thetas); r .*sind(thetas)];

    plot(pos(1, :), pos(2, :), 'DisplayName', planet_data.Name{1}, 'LineWidth', 2); hold on

    plot(pos(1,1), pos(2,1), '.k', 'MarkerSize', 20, 'HandleVisibility', 'off');
end