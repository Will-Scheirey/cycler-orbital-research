clear; clc; close all

AU = 1.496e+8;
mu = 132712000000;

ve = sqrt(mu / AU);

v_inf = 10;

ve_vec = [ve, 0, 0];
v_inf_vec = v_inf * [1, 1, 1] / sqrt(3);

draw_v_inf_sphere(ve_vec, v_inf_vec)

function draw_v_inf_sphere(ve_vec, v_inf_vec, names)

ve = norm(ve_vec);
v_inf = norm(v_inf_vec);

radius_1 = ve;
radius_2 = v_inf;

num_points = 30;

azimuth   = linspace(0, 2*pi, num_points);
elevation = linspace(-pi/2, pi/2, num_points);

[az, el] = meshgrid(azimuth, elevation);

r1 = radius_1 * ones(size(az));
r2 = radius_2 * ones(size(az));

x_const = radius_1 - radius_2^2 / (2*radius_1);

[x1, y1, z1] = sph2cart(az, el, r1);
x1(x1 > x_const) = x_const;

[x2, y2, z2] = sph2cart(az, el, r2);

x2 = x2 + radius_1;
x2(x2 < x_const) = x_const;

surf(x1, y1, z1, 'FaceAlpha', 0.1, 'FaceColor', 'blue', 'LineStyle','none', 'DisplayName', 'Primary Velocity Sphere'); hold on
surf(x2, y2, z2, 'FaceAlpha', 0.1, 'FaceColor', 'black', 'LineStyle','none', 'DisplayName', 'V_\infty Sphere')

x_intersect = ones(1, num_points) * x_const;
r_i2 = radius_1^2 - x_const^2;
r_i = sqrt(r_i2);

z_intersect = linspace(-r_i, r_i, num_points);
y_intersect1 = sqrt(r_i2 - z_intersect.^2);
y_intersect2 = -sqrt(r_i2 - z_intersect.^2);

plot3(x_intersect, y_intersect1, z_intersect, '-r', 'LineWidth', 5, 'HandleVisibility', 'off')
plot3(x_intersect, y_intersect2, z_intersect, '-r', 'LineWidth', 5, 'DisplayName', 'Resonance')

% quiver3(0, 0, 0, ve_vec(1), ve_vec(2), ve_vec(3), 0, 'LineWidth', 3, 'Color', 'blue', 'DisplayName', 'Earth Velocity')
quiver3(ve_vec(1), ve_vec(2), ve_vec(3), v_inf_vec(1), v_inf_vec(2), v_inf_vec(3),  0, 'LineWidth', 3, 'Color', 'black', 'DisplayName', 'V_\infty')

xlabel('x')
ylabel('y')
zlabel('z')

axis equal
legend
end