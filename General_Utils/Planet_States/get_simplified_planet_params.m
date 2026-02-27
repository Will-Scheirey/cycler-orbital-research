function data_out = get_simplified_planet_params()

mu_mercury    = 22030;                    % [km^3 s^-2] 
r_mercury     = 2440;                     % [km]
r_mercury_primary = 57.91e6;                  % [km]
mercury = struct('mu', mu_mercury, 'r', r_mercury, 'r_primary', r_mercury_primary);

mu_venus = 324900;                        % [km^3 s^-2] 
r_venus = 6051.8;                         % [km]
r_venus_primary = 108.2e6;                    % [km]
venus = struct('mu', mu_venus, 'r', r_venus, 'r_primary', r_venus_primary);

mu_earth = 398600;                        % [km^3 s^-2]
r_earth = 6378;                           % [km]
r_earth_primary = 149.6e6;                    % [km]
earth = struct('mu', mu_earth, 'r', r_earth, 'r_primary', r_earth_primary);

mu_mars = 42828;                          % [km^3 s^-2] 
r_mars = 3396;                            % [km]
r_mars_primary = 227.9e6;                     % [km]
mars = struct('mu', mu_mars, 'r', r_mars, 'r_primary', r_mars_primary);

mu_jupiter = 126686000;                   % [km^3 s^-2] 
r_jupiter = 71490;                        % [km]
r_jupiter_primary = 778.6e6;                  % [km]
jupiter = struct('mu', mu_jupiter, 'r', r_jupiter, 'r_primary', r_jupiter_primary);

mu_saturn = 37931000;                     % [km^3 s^-2] 
r_saturn = 60270;                         % [km]
r_saturn_primary = 1.433e9;                   % [km]
saturn = struct('mu', mu_saturn, 'r', r_saturn, 'r_primary', r_saturn_primary);

mu_primary = 132712000000;                    % [km^3 s^-2]

data_out = struct( ...
    'mu_primary', mu_primary, ...
    'mercury', mercury, ...
    'venus', venus, ...
    'earth', earth, ...
    'mars', mars, ...
    'jupiter', jupiter, ...
    'saturn', saturn);


end
