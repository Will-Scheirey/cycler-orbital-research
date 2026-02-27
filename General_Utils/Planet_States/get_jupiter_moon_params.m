function data_out = get_jupiter_moon_params()
    % Parameters from: https://www.princeton.edu/~willman/planetary_systems/Sol/Jupiter/
    % AND: https://www.nasa.gov/wp-content/uploads/2009/12/moons_of_jupiter_lithograph.pdf

    G = 6.674e-20;

    mu_jupiter = 1.89813e27 * G;
    
    mu_io = 8.93e22 * G;
    r_io = 1821.6;
    r_io_jupiter = 4.216e5;
    io = struct('mu', mu_io, 'r', r_io, 'r_primary', r_io_jupiter);

    mu_ganymede = 1.482e23 * G;
    r_ganymede = 2631;
    r_ganymede_jupiter = 1.070e6;
    ganymede = struct('mu', mu_ganymede, 'r', r_ganymede, 'r_primary', r_ganymede_jupiter);

    % mu_callisto = 1.076e23 * G;
    % r_callisto = 2410;
    % r_callisto_jupiter = 1.883e6;

    mu_europa = 4.80e22 * G;
    r_europa = 1560.8;
    r_europa_jupiter = 6.709e5;
    europa = struct('mu', mu_europa, 'r', r_europa, 'r_primary', r_europa_jupiter);

    data_out = struct( ...
        'mu_primary', mu_jupiter, ...
        'io', io, ...
        'ganymede', ganymede, ...
        'europa', europa);
end