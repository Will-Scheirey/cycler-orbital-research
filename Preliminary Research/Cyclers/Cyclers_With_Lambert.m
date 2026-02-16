clear; clc; close all;

mu_sun = 132712000000;                    % [km^3 s^-2]

r_earth_sun = 149.6e6;                    % [km]
r_mars_sun = 227.9e6;                     % [km]

T_earth = 2*pi * sqrt(r_earth_sun^3 / mu_sun);
T_mars = 2*pi * sqrt(r_mars_sun^3 / mu_sun);

T_ratio = T_mars / T_earth;

T_transfer = T_ratio * 3600*24*365;

S = 1 / (1/T_earth - 1/T_mars);

ns = 1:8;

thetas = linspace(0, 2*pi, 100);

max_revs = 5;

for i = 1:length(ns)
    subplot(length(ns), max_revs, 1+(i-1)*max_revs)
    ylabel("n = " + i)
    for j = 1:max_revs
        if i == 1
            subplot(length(ns), max_revs, 1+(i-1)*max_revs + j-1)
            title("R = " + (j-1))
        end

        figure(1)
        subplot(length(ns), max_revs, (i-1)*max_revs + j)
        hold on
        plot(r_earth_sun*cos(thetas), r_earth_sun*sin(thetas), 'b', 'LineWidth', 2)
        plot(r_mars_sun*cos(thetas), r_mars_sun*sin(thetas), 'r', 'LineWidth', 2)

        % plot(0, 0, '.y', 'MarkerSize', 80)

        set(gca, 'XTickLabel', {});
        set(gca, 'YTickLabel', {});

        hYLabel = get(gca, 'YLabel');
        set(hYLabel, 'Rotation', 0);
        axis equal
    end
end

for i = 1:length(ns)

    n = ns(i);

    tof = n*S;

    R1 = r_earth_sun * [1; 0; 0];
    omega_earth = 2*pi / T_earth;   % rad/s
    theta = omega_earth * tof;      % radians
    R2 = r_earth_sun * [cos(2*pi*tof/T_earth); sin(2*pi*tof/T_earth); 0];

    [v_D_list, v_A_list] = lambert(R1, R2, tof, mu_sun, 1200, 1000);

    sqr = ceil(sqrt(length(v_D_list)));
    sqr1 = ceil(length(v_D_list) / sqr);

    for z=1:width(v_D_list)
        figure(1)

        v_D = v_D_list{z};

        T = orbital_period(R1, v_D_list{z}, mu_sun);
        revs = tof / T;
        if imag(tof) ~= 0 || revs > max_revs
            continue
        end

        subplot(length(ns), max_revs, 1+(i-1)*max_revs + floor(revs))

        plot_orbit_from_state(R1, v_D, mu_sun, '-k', 1); hold on

        axis equal

        axis('padded');
    end
end

sgtitle("Earth Cycler Orbits from Lambert's Problem")
fontsize(16,"points")


function T = orbital_period(r_vec, v_vec, mu)
r = norm(r_vec);
v = norm(v_vec);
epsilon = v^2/2 - mu/r;
a = -mu/(2*epsilon);
T = 2*pi*sqrt(a^3/mu);
end