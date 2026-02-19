function [mu_sun, mu_earth, mu_mars, r_earth_sun, r_mars_sun, T_earth, T_mars, S, angle_turn] = get_earth_mars_params()
mu_sun = 132712000000;                    % [km^3 s^-2]
mu_earth = 398600;
mu_mars  = 42828;

% Force the orbits of Earth and Mars to line up nicely

T_desired_earth = 3600*24*365;
T_desired_mars  = 15/8 * T_desired_earth;

r_earth_sun = ((T_desired_earth/(2*pi))^2 * mu_sun)^(1/3);
r_mars_sun =  ((T_desired_mars/ (2*pi))^2 * mu_sun)^(1/3);

T_earth = 2*pi * sqrt(r_earth_sun^3 / mu_sun);
T_mars = 2*pi * sqrt(r_mars_sun^3 / mu_sun);

S = 1 / (1/T_earth - 1/T_mars);

angle_turn = mod(T_mars/T_earth * 2*pi, 2*pi);
if angle_turn > pi
    angle_turn = angle_turn - 2*pi;
end
end