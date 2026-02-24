clear; clc; close all
[mu_sun, mu_earth, mu_mars,... 
    r_earth_sun, r_mars_sun,... 
    T_earth, T_mars, S, ~] = get_earth_mars_params();

dv_tot = 10; % [km]

num_angles = 100;
thetas = linspace(0, pi, num_angles);

v_earth_sun = sqrt(mu_sun / r_earth_sun);

v_earth_sun_vec = [0, v_earth_sun, 0];
r_earth_sun_vec = [r_earth_sun, 0, 0];

as = zeros(num_angles, 1);

ts = zeros(num_angles, 1);

thetas_orbit = 0:0.01:2*pi;

zlim = 1200;
numz = 1000;

for idx = 1:num_angles
    theta = thetas(idx);

    v_vec = dv_tot * [cos(theta), sin(theta), 0] + v_earth_sun_vec;
    v = norm(v_vec);

    r_vec = r_earth_sun_vec;
    r = norm(r_vec);

    h_vec = cross(r_vec, v_vec);
    h = norm(h_vec);
    eps = v^2/2 - mu_sun / r_earth_sun;

    a = -mu_sun/(2*eps);

    if a < r_mars_sun
        continue
    end

    e_vec = (cross(v_vec, h_vec) / mu_sun) - (r_vec / r);
    e = norm(e_vec);
    theta0 = atan2(dot(cross(e_vec, r_vec), h_vec)/h, dot(e_vec, r_vec));  % radians    

    theta_encounter = acos((h^2/(mu_sun * r_mars_sun) - 1)/e);
    tof_encounter = calc_tof(a, mu_sun, e, theta_encounter - theta0);

    r_encounter = h^2/mu_sun * 1/(1 + e*cos(theta_encounter));

    tof_1 = S - tof_encounter;

    R11 = r_encounter * [cos(theta_encounter - theta0); sin(theta_encounter - theta0); 0];
    R21 = r_earth_sun * [cos(2*pi*S/T_earth); sin(2*pi*S/T_earth); 0];

    [v_D1, v_A1] = lambert(R11, R21, tof_1, mu_sun, zlim, numz);
    v_D1 = cell2mat(v_D1); v_A1 = cell2mat(v_A1);

    %% Plotting

    thetas_encounter = theta0:0.01:theta_encounter;

    ts(idx) = theta_encounter;
    as(idx) = a;

    r_orbit = h^2/mu_sun * 1./(1 + e*cos(thetas_encounter));

    clf
    plot(r_earth_sun * cos(thetas_orbit), r_earth_sun * sin(thetas_orbit), '-b', 'LineWidth', 1.5); hold on
    plot(r_mars_sun * cos(thetas_orbit), r_mars_sun * sin(thetas_orbit), '-r', 'LineWidth', 1.5);
    
    x_orbit = r_orbit .* cos(thetas_encounter - theta0);
    y_orbit = r_orbit .* sin(thetas_encounter - theta0);

    num_orbits = width(v_D1);
    for k = 1:num_orbits
        plot_lambert_arc_time(R11, v_D1(:, k), mu_sun, tof_1, 'k');
    end

    plot(x_orbit, y_orbit, '-g', 'LineWidth', 1.5);
    xlim([-1, 1] * r_mars_sun * 2)
    ylim([-1, 1] * r_mars_sun * 2)

    axis square

    drawnow
    pause(0.1)
end