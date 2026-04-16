clear; clc; close all
load_params

%% Run Algorithm

p_targ = 1; % Number of synodic periods for the cycler
h_targ = 3; % Number of half-years
s_targ = 1; % Number of identical generic returns
i_targ = 5; % Number of revolutions made by generic return

fast = true;
long = false;

solution_str = '2.5.1.+0';
% solution_str = '4.3.1.-5';

str_parts = split(solution_str, '.');

p_targ = str2double(str_parts(1));
h_targ = str2double(str_parts(2));
s_targ = str2double(str_parts(3));

fast = contains(str_parts(4), '-');

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

plot_intersection(r0, lim/20, "0")
plot_intersection(r_intersect_vec, lim/20, "1")
plot_intersection(r_end, lim/20, "end")
plot_intersection(r0_2, lim/20, "Mars Start")

view(2)

title(sprintf("Solutions for %s", solution_str));

axis equal

% return
figure(2)
C = solution.rotm;

min_z = inf;
for n=1:length(pos_all)
    local = solution.v_local{n};
    z = local(3);
    min_z = min(min_z, z);
end

v_inf_norm = norm(solution.v_inf_minus{1}{1});
draw_sphere(v_inf_norm, [0, 0, 0], [-inf,inf], [-inf,inf], [min_z,inf], 'FaceAlpha', 0.1, 'FaceColor', 'blue', 'LineStyle','none', 'DisplayName', 'Primary Velocity Sphere')

hold on;

num_orbits = length(pos_all);
colors = turbo(num_orbits);
for n = 1:num_orbits
    r0_sc = pos_all{n}{1} / AU;
    v0_sc = pos_all{n}{2};
    tof   = pos_all{n}{3};
    C = solution.rotm;

    v_inf = solution.v_inf_minus{1}{n};

    v = solution.v_local{n};

    fprintf("Segment %d\n" + ...
        "\tr0:    (%0.2f, %0.2f, %0.2f) AU; mag %0.2f\n" + ...
        "\tv0:    (%0.2f, %0.2f, %0.2f) km/s; mag %0.2f\n"+ ...
        "\tv_inf: (%0.2f, %0.2f, %0.2f) km/s; mag %0.2f\n"+ ...% "\tv_mag: %0.2f km/s"+ ...
        "\tv_inf_FRAME: (%0.2f, %0.2f, %0.2f) km/s; mag %0.2f\n"+ ...% "\tv_mag: %0.2f km/s"+ ...
        "\n", n, r0_sc(1), r0_sc(2), r0_sc(3), norm(r0_sc), ...
        v0_sc(1), v0_sc(2), v0_sc(3), norm(v0_sc), ...
        v_inf(1), v_inf(2), v_inf(3), norm(v_inf), ...
        v(1), v(2), v(3), norm(v));

    figure(2)
    quiver3(0, 0, 0, v(1), v(2), v(3), 'LineWidth', 3, 'Color', colors(n, :), 'DisplayName', sprintf('V inf %d', n), 'AutoScale', 'off');
    hold on
end

figure(2)
legend
axis equal
xlabel("X")
ylabel("Y")
zlabel("Z (Earth Velocity)")
title("Local Frame")

return
figure(2);
clf

num_orbits = length(pos_all);

t0 = 0;

rel_frame_all = [0, 0, 0];
inertial_frame_all = [0, 0, 0];
figure(2)
clf

num_orbits = length(pos_all);
t0 = 0;
num_plot = 1000;

for n = 1:num_orbits
    r0_sc = pos_all{n}{1};
    v0_sc = pos_all{n}{2};
    tof   = pos_all{n}{3};

    t = linspace(0, tof, num_plot);

    x_r = zeros(size(t));
    y_r = zeros(size(t));
    z_r = zeros(size(t));

    p_sc_all = zeros(length(t), 3);

    for k = 1:length(t)
        tk = t(k) + t0;

        r_sc = get_planet_pos(r0_sc, v0_sc, mu, t(k));
        r_e  = get_planet_pos(r0,   v0,   mu, tk);
        r_m  = get_planet_pos(r0_2, v0_2, mu, tk);

        p_sc_all(k, :) = r_sc;

        x_sc = r_sc(1); y_sc = r_sc(2); z_sc = r_sc(3);
        x_e  = r_e(1);  y_e  = r_e(2);
        x_m  = r_m(1);  y_m  = r_m(2);

        r_e  = [x_e; y_e; 0];
        r_m  = [x_m; y_m; 0];
        r_sc = [x_sc; y_sc; z_sc];

        r_em = r_m - r_e;
        d_em = norm(r_em);

        xhat = r_em / d_em;
        yhat = [-xhat(2); xhat(1); 0];

        r_rel = r_sc - r_e;

        x_r(k) = dot(r_rel, xhat) / d_em;
        y_r(k) = dot(r_rel, yhat) / d_em;
        z_r(k) = r_rel(3) / d_em;
    end
    
    inertial_frame_all = vertcat(inertial_frame_all, p_sc_all);

    rel_frame_all = vertcat(rel_frame_all, [x_r', y_r', z_r']);

    figure(3)
    plot3(x_r, y_r, z_r, '-k', 'Clipping','off', 'HandleVisibility', 'off'); hold on

    t0 = t0 + tof;
end


rel_frame_all = rel_frame_all(2:end, :);
inertial_frame_all = inertial_frame_all(2:end, :);


figure(2);
clf
plot3(rel_frame_all(:, 1), rel_frame_all(:, 2), rel_frame_all(:, 3), '.k', 'Clipping','off', 'DisplayName', 'Trajectory'); hold on
plot3(0, 0, 0, '.b', 'MarkerSize', 50, 'DisplayName', 'Earth')
plot3(1, 0, 0, '.r', 'MarkerSize', 50, 'DisplayName', 'Mars')
legend
axis equal

lim = 3;
xlim([-1,1]*lim)
ylim([-1,1]*lim)
zlim([-1,1]*lim)

figure(3);
clf
plot3(inertial_frame_all(:, 1), inertial_frame_all(:, 2), inertial_frame_all(:, 3), '.k', 'Clipping','off', 'DisplayName', 'Trajectory'); hold on

legend
axis equal

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
