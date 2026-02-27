function [mu_sun, r_mercury, r_venus, r_earth, r_mars, r_jupiter, r_saturn] = get_circular_planets()
mu_sun = 132712000000;                    % [km^3 s^-2]

params = solar_orbit_params();

r_mercury = params.Mercury.a / 1e3; % [km]
r_venus   = params.Venus.a / 1e3;   % [km]
r_earth   = params.EarthMoonBarycenter.a / 1e3;   % [km]
r_mars    = params.Mars.a / 1e3;    % [km]
r_jupiter = params.Jupiter.a / 1e3; % [km]
r_saturn  = params.Saturn.a / 1e3;  % [km]

end