clear; clc; close all
load_params

%% Run Algorithm
% solution_str = '2.5.1.+0';
solution_str = '1.0.1.-1'

str_parts = split(solution_str, '.');

p_targ = abs(str2double(str_parts(1)));
h_targ = str2double(str_parts(2));
s_targ = str2double(str_parts(3));

fast = contains(str_parts(4), '-');
long = contains(str_parts(1), '-');

i_targ = str2double(erase(str_parts(4), '-'));

[N_max, tof, rev_all] = russell_algorithm(tau, p_targ, h_targ, s_targ, mean_motion, r_b1, v_b1, r_b2, r0, v0, mu_b1, mu, rp_min, TR_min, AR_min, year2sec);

data = struct('p', p_targ, 'h', h_targ, 's', s_targ, 'N_max', N_max, 'tof', tof, 'rev_data', rev_all);

feasible_solutions = get_feasible_solutions({data}, true);

target_solutions = get_solution(p_targ, h_targ, s_targ, i_targ, feasible_solutions, fast, long);

theta_b1_turn = mean_motion * tof;

%% CALCULATE INTERSECTIONS
vd = target_solutions{1}.vd;

[theta_intersect, r_intersect] = get_intersection(r0, vd, mu, r_b2);

tof_intersect = tof_from_theta(r0, vd, mu, theta_intersect);

mars_theta0 = theta_intersect - mean_motion_2*tof_intersect;

r0_2 = r_b2 * [cos(mars_theta0); sin(mars_theta0); 0];
v0_2 = v_b2 * [-sin(mars_theta0); cos(mars_theta0); 0];

r_intersect_vec = r_intersect * [cos(theta_intersect); sin(theta_intersect); 0];

theta_end = p_targ * angle_turn;

r_end = r_b1 * [cos(theta_end); sin(theta_end); 0];

%% PLOT

figure(1)
clf

plot_orbit(r0, v0, mu, 'DisplayName', 'Earth', 'LineStyle', '--', 'LineWidth', 2); hold on
plot_orbit(r0_2, v0_2, mu, 'DisplayName', 'Mars', 'LineStyle', '--', 'LineWidth', 2);
hold on

lim = r_b2 * 1.1;

solution = target_solutions{1};
[orbits, pos_all] = plot_solution(solution, mu, vd, r_b1, v_b1, theta_b1_turn, lim/20);

plot_intersection(r0, lim/20, "0", 'b')
plot_intersection(r_intersect_vec, lim/20, "1", 'r')
plot_intersection(r_end, lim/20, "end")
plot_intersection(r0_2, lim/20, "Mars Start")

view(2)

title(sprintf("Solutions for %s", solution_str));
grid on
axis equal
xlabel("X");
ylabel("Y");
zlabel("Z");

figure(2)
clf

C = solution.rotm;

% Flatten the nested per-s storage into one list for plotting/debug
v_local_flat = {};
v_inf_flat   = {};

for s_idx = 1:solution.s
    v_local_s = solution.v_local{s_idx};
    v_inf_s   = solution.v_inf_minus{s_idx};

    for k = 1:length(v_local_s)
        v_local_flat{end+1,1} = v_local_s{k};
        v_inf_flat{end+1,1}   = v_inf_s{k};
    end
end

min_z = inf;
for n = 1:length(v_local_flat)
    v_local = v_local_flat{n};
    min_z = min(min_z, v_local(3));
end

draw_sphere(v_b1, [0, 0, 0], [-inf, inf], [-inf, inf], [-inf, inf], ...
    'FaceAlpha', 0.1, ...
    'FaceColor', 'blue', ...
    'LineStyle', 'none', ...
    'DisplayName', 'Primary Velocity Sphere');

v_inf_norm = norm(v_inf_flat{1});
draw_sphere(v_inf_norm, [0, 0, v_b1], [-inf, inf], [-inf, inf], [-inf, inf], ...
    'FaceAlpha', 0.1, ...
    'FaceColor', [0.2, 0.2, 0.2], ...
    'LineStyle', 'none', ...
    'DisplayName', 'V_\infty Sphere');

quiver3(0, 0, 0, 0, 0, v_b1, 'Color', 'blue', 'AutoScale', 'off', 'LineWidth', 3, 'DisplayName', 'Earth Velocity')

plot3(0, 0, 0, '.b', 'MarkerSize', 30, 'HandleVisibility', 'off')

hold on

num_vec = length(v_local_flat);
colors = turbo(num_vec);

num_print = min(length(pos_all), num_vec);

for n = 1:num_print
    r0_sc = pos_all{n}{1} / AU;
    v0_sc = pos_all{n}{2};
    tof   = pos_all{n}{3};

    v_inf = v_inf_flat{n};
    v     = v_local_flat{n};

    fprintf("Segment %d\n" + ...
        "\tr0:    (%0.2f, %0.2f, %0.2f) AU; mag %0.2f\n" + ...
        "\tv0:    (%0.2f, %0.2f, %0.2f) km/s; mag %0.2f\n" + ...
        "\tv_inf: (%0.2f, %0.2f, %0.2f) km/s; mag %0.2f\n" + ...
        "\tv_inf_FRAME: (%0.2f, %0.2f, %0.2f) km/s; mag %0.2f\n\n", ...
        n, ...
        r0_sc(1), r0_sc(2), r0_sc(3), norm(r0_sc), ...
        v0_sc(1), v0_sc(2), v0_sc(3), norm(v0_sc), ...
        v_inf(1), v_inf(2), v_inf(3), norm(v_inf), ...
        v(1), v(2), v(3), norm(v));

    quiver3(0, 0, v_b1, v(1), v(2), v(3), ...
        'LineWidth', 3, ...
        'Color', colors(n, :), ...
        'DisplayName', sprintf('V inf %d', n), ...
        'AutoScale', 'off');
    hold on
end

legend
axis equal
xlabel("X")
ylabel("Y")
zlabel("Z (Earth Velocity)")
title("Inertial Frame")

num_orbits = length(pos_all);
num_plot   = 1000;

t0 = 0;

rel_frame_all      = cell(num_orbits, 1);
inertial_frame_all = cell(num_orbits, 1);

for n = 1:num_orbits
    r0_sc = pos_all{n}{1};
    v0_sc = pos_all{n}{2};
    tof   = pos_all{n}{3};

    t = linspace(0, tof, num_plot);

    rel_seg      = zeros(num_plot, 3);
    inertial_seg = zeros(num_plot, 3);

    for k = 1:num_plot
        tk = t(k) + t0;

        r_sc = get_planet_pos(r0_sc, v0_sc, mu, t(k));
        r_e  = get_planet_pos(r0,   v0,   mu, tk);
        r_m  = get_planet_pos(r0_2, v0_2, mu, tk);

        inertial_seg(k, :) = r_sc(:)';

        r_em = r_m - r_e;
        d_em = norm(r_em);

        xhat = r_em / d_em;
        yhat = [-xhat(2); xhat(1); 0];
        zhat = [0; 0; 1];

        r_rel = r_sc - r_e;

        rel_seg(k, 1) = dot(r_rel, xhat) / d_em;
        rel_seg(k, 2) = dot(r_rel, yhat) / d_em;
        rel_seg(k, 3) = dot(r_rel, zhat) / d_em;
    end

    inertial_frame_all{n} = inertial_seg;
    rel_frame_all{n}      = rel_seg;

    t0 = t0 + tof;
end

return
%% Plot pulsating Earth/Mars-fixed frame
figure(3)
clf
hold on

colors = turbo(num_orbits);

for n = 1:num_orbits
    rel_seg = rel_frame_all{n};
    plot3(rel_seg(:,1), rel_seg(:,2), rel_seg(:,3), ...
        'LineWidth', 2, ...
        'Color', colors(n,:), ...
        'DisplayName', sprintf('Trajectory %d', n), ...
        'Clipping', 'off');
end

plot3(0, 0, 0, '.b', 'MarkerSize', 50, 'DisplayName', 'Earth')
plot3(1, 0, 0, '.r', 'MarkerSize', 50, 'DisplayName', 'Mars')

xlabel('x / |r_{EM}|')
ylabel('y / |r_{EM}|')
zlabel('z / |r_{EM}|')
title('Earth/Mars-Fixed Pulsating Frame')
legend
axis equal

lim = 3;
xlim([-1, 1] * lim)
ylim([-1, 1] * lim)
zlim([-1, 1] * lim)

view(3)
grid on

%% Plot inertial frame
figure(4)
clf
hold on

for n = 1:num_orbits
    inertial_seg = inertial_frame_all{n};
    plot3(inertial_seg(:,1), inertial_seg(:,2), inertial_seg(:,3), ...
        'LineWidth', 2, ...
        'Color', colors(n,:), ...
        'DisplayName', sprintf('Trajectory %d', n), ...
        'Clipping', 'off');
end

plot_orbit(r0,   v0,   mu, 'DisplayName', 'Earth Orbit', 'LineStyle', '--', 'LineWidth', 1.5);
plot_orbit(r0_2, v0_2, mu, 'DisplayName', 'Mars Orbit',  'LineStyle', '--', 'LineWidth', 1.5);

plot3(r0(1),   r0(2),   r0(3),   '.b', 'MarkerSize', 40, 'DisplayName', 'Earth Start')
plot3(r0_2(1), r0_2(2), r0_2(3), '.r', 'MarkerSize', 40, 'DisplayName', 'Mars Start')

xlabel('X')
ylabel('Y')
zlabel('Z')
title('Inertial Frame')
legend
axis equal
view(3)
grid on

function draw_sphere(radius, origin, x_lim, y_lim, z_lim, varargin)
num_points = 200;
azimuth   = linspace(0, 2*pi, num_points);
elevation = linspace(-pi/2, pi/2, num_points);

[az, el] = meshgrid(azimuth, elevation);

[x, y, z] = sph2cart(az, el, radius);
x(x < x_lim(1)) = x_lim(1); x(x > x_lim(2)) = x_lim(2);
y(y < y_lim(1)) = y_lim(1); y(y > y_lim(2)) = y_lim(2);
z(z < z_lim(1)) = z_lim(1); z(z > z_lim(2)) = z_lim(2);

x = x + origin(1);
y = y + origin(2);
z = z + origin(3);

surf(x, y, z, varargin{:}); hold on
end
