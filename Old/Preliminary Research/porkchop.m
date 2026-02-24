clear; clc; close all;

mu_sun = 132712000000;                    % [km^3 s^-2]

mu_earth = 398600;                        % [km^3 s^-2]
r_earth = 6378;                           % [km]
r_earth_sun = 149.6e6;                    % [km]
v_earth_sun = sqrt(mu_sun / r_earth_sun);  % [km s^-1]
i_earth_sun = 7.25; % Inclination of Earth Orbit [deg]

mu_mars = 42828;                          % [km^3 s^-2] 
r_mars = 3396;                            % [km]
r_mars_sun = 227.9e6;                     % [km]
v_mars_sun = sqrt(mu_sun / r_mars_sun);  % [km s^-1]
i_mars_sun = 1.85; % Inclination of Mars Orbit [deg]

period_earth = 2*pi*sqrt(r_earth_sun^3 / mu_sun); 
period_mars = 2*pi*sqrt(r_mars_sun^3 / mu_sun);

synodic_period = 1 / (1/period_earth - 1/period_mars);

num = 100;
tof_all = linspace(period_mars/10, period_mars, num);

start_all = linspace(0, 1, num) * synodic_period*2;

c3 = NaN(num, num);

for i=1:num
    tof = tof_all(i);
    for n = 1:num
        start = start_all(n);
        disp("i: " + i + "; n: " + n)

        n_earth = 2*pi/period_earth;
        n_mars = 2*pi/period_mars;

        theta_earth_i = n_earth*start;
        theta_mars_f = n_mars*(start+tof);

        r_vec_i = r_earth_sun * [cos(theta_earth_i)*sind(i_earth_sun), sin(theta_earth_i)*sind(i_earth_sun), cosd(i_earth_sun)];
        r_vec_f = r_mars_sun * [cos(theta_mars_f)*sind(i_mars_sun), sin(theta_mars_f)*sind(i_mars_sun), cosd(i_mars_sun)];

        [a, p, v_D, v_A] = lambert_solver(r_vec_i, r_vec_f, tof, mu_sun, "1A");
        
        v_earth_vec = v_earth_sun * [-sin(theta_earth_i), cos(theta_earth_i), 0];
        dv_1 = norm(v_D - v_earth_vec);

        v_mars_vec = v_mars_sun * [-sin(theta_mars_f), cos(theta_mars_f), 0];
        dv_2 = norm(v_A - v_mars_vec);

        dv_A = dv_1 + dv_2;

        [a, p, v_D, v_A] = lambert_solver(r_vec_i, r_vec_f, tof, mu_sun, "1B");
        
        dv_1 = norm(v_D - v_earth_vec);

        dv_2 = norm(v_A - v_mars_vec);

        dv_B = dv_1 + dv_2;

        [a, p, v_D, v_A] = lambert_solver(r_vec_i, r_vec_f, tof, mu_sun, "2A");
        
        dv_1 = norm(v_D - v_earth_vec);

        dv_2 = norm(v_A - v_mars_vec);

        dv_A1 = dv_1 + dv_2;

        [a, p, v_D, v_A] = lambert_solver(r_vec_i, r_vec_f, tof, mu_sun, "2B");
        
        dv_1 = norm(v_D - v_earth_vec);

        dv_2 = norm(v_A - v_mars_vec);

        dv_B1 = dv_1 + dv_2;

        [a, p, v_D, v_A] = lambert_solver(r_vec_i, r_vec_f, tof, mu_sun, "1H");
        
        dv_1 = norm(v_D - v_earth_vec);

        dv_2 = norm(v_A - v_mars_vec);

        dv_H1 = dv_1 + dv_2;

        [a, p, v_D, v_A] = lambert_solver(r_vec_i, r_vec_f, tof, mu_sun, "2H");
        
        dv_1 = norm(v_D - v_earth_vec);

        dv_2 = norm(v_A - v_mars_vec);

        dv_H2 = dv_1 + dv_2;

        dv = min([dv_A, dv_B, dv_A1, dv_B1, dv_H1, dv_H2]);
    
        if dv < 15
            c3(i, n) = dv;
        end
    end
end

%% Plot

figure(1)
clf
contourf(tof_all / 3600 / 24, start_all, c3.')
xlabel("TOF (day)")
ylabel("Start (s)")
zlabel("C3")