clear;
[consts, mu_altaira, planet_data] = load_problem_data();
clc;
close all;

num_planets = height(planet_data);


figure(1)
clc
for i = 1:num_planets
    a = planet_data.a_km(i);
    e = planet_data.e(i);
    inc = planet_data.i_deg(i);
    w = planet_data.w_deg(i);
    RAAN = planet_data.RAAN_deg(i);
    theta0 = planet_data.theta0_deg(i);

    if e == 0
        truelon = mod(RAAN + theta0, 360);
        [r, v] = keplerian2ijk(a, e, inc, RAAN, w, theta0, 'GravitationalParameter', mu_altaira, 'truelon', truelon);
    else
        lonper = mod(RAAN + w, 360);
        [r, v] = keplerian2ijk(a, e, inc, RAAN, w, theta0, 'GravitationalParameter', mu_altaira, 'lonper', lonper);
    end
    plot3(r(1), r(2), r(3), '.', 'MarkerSize', 40, 'Clipping', 'off', 'DisplayName', planet_data.Name{i}); hold on

    mag = norm(v);
    u = v / norm(v);
    vec = u * 1e10;

    quiver3(r(1), r(2), r(3), vec(1), vec(2), vec(3), 'Clipping', 'off', 'LineWidth', 2)

    z = linspace(0, 360, 1000);

    theta = mod(theta0 + z, 360);
    if e == 0
        truelon = mod(RAAN + theta, 360);
        [r, ~] = keplerian2ijk(a, e, inc, RAAN, w, theta, 'GravitationalParameter', mu_altaira, 'truelon', truelon);
    else
        lonper = mod(RAAN + w, 360);
        [r, ~] = keplerian2ijk(a, e, inc, RAAN, w, theta, 'GravitationalParameter', mu_altaira, 'lonper', lonper);
    end
    plot3(r(1, :), r(2, :), r(3, :), '.', 'MarkerSize', 5, 'Clipping', 'off'); hold on
end

plane = linspace(-1e11, 1e11, 100);

[x, y] = meshgrid(plane);

z = zeros(size(x, 1));

surf(x, y, z, 'LineStyle', 'none', 'FaceAlpha', 0.1)

ax = gca;
% ax.Projection = 'perspective';