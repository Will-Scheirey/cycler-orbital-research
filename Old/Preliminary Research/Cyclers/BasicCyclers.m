clear; clc; close all;
mu_sun = 132712000000;                    % [km^3 s^-2]

r_earth_sun = 149.6e6;                    % [km]
r_mars_sun = 227.9e6;                     % [km]

T_earth = 2*pi * sqrt(r_earth_sun^3 / mu_sun);
T_mars = 2*pi * sqrt(r_mars_sun^3 / mu_sun);

T_ratio = T_mars / T_earth;

T_transfer = T_ratio * 3600*24*365;

synodic_period = 1 / (1/T_earth - 1/T_mars);

angle_turn = mod(T_ratio * 360, 360);
if angle_turn > 180
    angle_turn = angle_turn - 360;
end

[a_transfer, e_transfer] = solve_transfer_geometry(r_earth_sun, T_transfer, mu_sun, 0);

eps_transfer = -mu_sun/(2*a_transfer);
h_transfer = sqrt(-1/(2*eps_transfer)*mu_sun^2*(1-e_transfer^2));

v0_transfer = mu_sun/h_transfer *[
    e_transfer*sin(0)
    1 + e_transfer*cos(0)
]

period = T_earth * (ceil(T_ratio));

num_steps = 300;
ts = linspace(0, synodic_period, num_steps);

thetas = linspace(0, 2*pi, 100);
r_transfer = a_transfer * (1 - e_transfer^2) ./ (1 + e_transfer * cos(thetas));

r_transfer(1)

% return

figure(1)

plot(r_earth_sun * cos(thetas), r_earth_sun * sin(thetas), '-b', 'LineWidth', 1.5); hold on
plot(r_mars_sun * cos(thetas), r_mars_sun * sin(thetas), '-r', 'LineWidth', 1.5);

for z = 1:abs(angle_turn):360
    plot(r_transfer .* cos(thetas), r_transfer .* sin(thetas), '-g', 'LineWidth', 1.5);
end

legend("Earth", "Mars", "Spacecraft")

axis equal
% axis square


precession_angle = 0;
t0 = 0;
last_encounter = 0;

outputVideo = VideoWriter('myVideo.mp4'); % Specify filename and format
outputVideo.FrameRate = num_steps / 10;
open(outputVideo);

for i = 1:1:length(ts)
    t = ts(i);
    
    angle = @(t, T) wrapToPi(t/T * 2*pi);

    angle_earth = angle(t, T_earth);
    angle_mars = angle(t, T_mars);
    angle_transfer = true_anomaly(t - t0, e_transfer, T_transfer);

    diff_turn = abs(angle_earth - abs(deg2rad(angle_turn)));
    diff_transfer = abs(angle_transfer - angle_earth);
    %{
    if mod(t, synodic_period) < 1e6
        if t - last_encounter > T_earth
        angle_transfer = true_anomaly(t - t0, e_transfer, T_transfer);
        precession_angle = precession_angle - angle_transfer;
        last_encounter = t;
        end
    end
    %}
    r_transfer = a_transfer * (1 - e_transfer^2) ./ (1 + e_transfer * cos(angle_transfer + precession_angle));

    clf
    % set(gca,'Units','pixels','Position', [20 20 100 150])

    P_earth = [r_earth_sun * cos(angle_earth), r_earth_sun * sin(angle_earth)];
    P_mars = [r_mars_sun * cos(angle_mars), r_mars_sun * sin(angle_mars)];

    plot(P_earth(1), P_earth(2), '.', 'Color', [125, 201, 240]/255, 'MarkerSize', 55); hold on;
    plot(P_mars(1), P_mars(2), '.', 'Color',     [231, 173, 138]/255, 'MarkerSize', 40)
    plot(r_transfer * cos(angle_transfer), r_transfer * sin(angle_transfer), '.', 'Color',         [202, 116, 199]/255, 'MarkerSize', 35)

    plot(0, 0, '.', 'Color', [252, 243, 81]/255, 'MarkerSize', 50)

    % line([P_earth(1), P_mars(1)], [P_earth(2), P_mars(2)], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5)
    % line([P_earth(1), 0], [P_earth(2), 0], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5)

    legend("Earth", "Mars", "Spacecraft")
    xlim([-1,1] * a_transfer * 1.5);
    ylim([-1,1] * a_transfer * 1.5);
    axis square
    drawnow;

    frame = getframe(gcf); % captures the current figure (gcf)
    writeVideo(outputVideo, frame);
end
close(outputVideo);

function theta = true_anomaly(t, e, period)
    n = 2*pi/period;
    
    % Calculate mean anomaly at time t (M)
    M = wrapTo2Pi(n*t); % Use wrapTo2Pi for angles in [0, 2pi]
    
    % Solve Kepler's Equation for Eccentric Anomaly (E) using fsolve or a custom function
    % Here, we'll use the formula derived from the tangent half-angle relation
    E = fsolve(@(x) x - e*sin(x) - M, pi, optimset('Display','off'));
    
    % Calculate True Anomaly (nu)
    theta = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(E/2));
    
    % Ensure angle is in the range [0, 2pi]
    theta = wrapTo2Pi(theta);
end

function [a, e] = solve_transfer_geometry(r0, T, mu, theta0)
       
    a = ((T/(2*pi))^2 * mu)^(1/3);

    syms e_sym;
    eqn = r0 == a * (1 - e_sym^2) / (1 + e_sym * cosd(theta0));
    e_sol = double(solve(eqn, e_sym));
    
    e = e_sol(find(e_sol > 0, 1));
end