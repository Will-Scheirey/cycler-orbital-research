clear; clc; close all;

mu_sun = 132712000000;                    % [km^3 s^-2]
mu_earth = 398600;

% Force the orbits of Earth and Mars to line up nicely

T_desired_earth = 3600*24*365;
T_desired_mars  = 15/8 * T_desired_earth;

r_earth_sun = ((T_desired_earth/(2*pi))^2 * mu_sun)^(1/3);
r_mars_sun =  ((T_desired_mars/ (2*pi))^2 * mu_sun)^(1/3);

%{
r_earth_sun = 149.6e6;                    % [km]
r_mars_sun = 227.9e6;                     % [km]
%}

T_earth = 2*pi * sqrt(r_earth_sun^3 / mu_sun);
T_mars = 2*pi * sqrt(r_mars_sun^3 / mu_sun);

T_ratio = T_mars / T_earth;

T_transfer = T_ratio * 3600*24*365;

S = 1 / (1/T_earth - 1/T_mars);

angle_turn = mod(T_ratio * 2*pi, 2*pi);
if angle_turn > pi
    angle_turn = angle_turn - 2*pi;
end

thetas = linspace(0, 2*pi, 100);

figure(1)
hold on
plot(r_earth_sun*cos(thetas), r_earth_sun*sin(thetas), 'b', 'LineWidth', 2)
plot(r_mars_sun*cos(thetas), r_mars_sun*sin(thetas), 'r', 'LineWidth', 2)

ns = 1:7;

for i = 1:length(ns)
    n = ns(i);

    tof = n*S;

    R1 = r_earth_sun * [1; 0; 0];
    omega_earth = 2*pi / T_earth;   % rad/s
    theta = omega_earth * tof;      % radians
    R2 = r_earth_sun * [cos(2*pi*tof/T_earth); sin(2*pi*tof/T_earth); 0];

    [v_D_list, v_A_list] = lambert(R1, R2, tof, mu_sun, 1200, 1000);

    for z = 1:numel(v_D_list)
        V1 = v_D_list{z};
        good = eval_transfer(R1, V1, r_mars_sun, mu_sun, mu_earth, 6378, angle_turn);

        if good
            plot_orbit_from_state(R1, V1, mu_sun, '-k', 0.2); hold on
        end
    end
end

axis equal
axis padded
title(sprintf("Cycler Orbits for n=%d to %d", ns(1), ns(end)))

function good = eval_transfer(R1, V1, apo_min, mu, mu0, rp_min, angle_turn)
    r1 = norm(R1);
    v1 = norm(V1);

    h_vec = cross(R1, V1);
    h = norm(h_vec);
    
    e_vec = cross(V1, h_vec)/mu - R1/r1;
    e = norm(e_vec);

    a = h^2/mu * 1/(1-e^2);
    
    r_a = abs(h^2/mu * 1/(e-1));

    
    if r_a < apo_min
        good = false;
        return
    end
    
    
    gamma = atan(e*r1 / (a*(1 - e^2)) * sin(angle_turn));

    V = sqrt(mu/r1 * (2 - r1/a));
    VO = sqrt(mu / r1);

    DV = abs(V*sin(gamma));

    V_inf = sqrt(V^2 + VO^2 - 2*V*VO*cos(gamma));

    r_p = mu0 * (2/(DV*V_inf) - 1/V_inf^2);

    good = r_p > rp_min;  
end