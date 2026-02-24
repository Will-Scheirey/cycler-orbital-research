% mu = 132712000000;

[mu, mu_b1, mu_b2, r_b1, r_b2, T_1, T_2, tau, angle_turn] = get_ideal_earth_mars_params();
angle_turn = -angle_turn;

sec2day = 1 / (3600 * 24);
sec2year = sec2day / 365;
year2sec = 1 / sec2year;

% AU = 149597871; % [km]
AU = r_b1;

% r_b1 = AU; % Earth
v_b1 = sqrt(mu / r_b1); % Velocity of Earth
% mu_b1 = 3.986004418e5;
rp_min = 6378;

% r_b2 = 2.2737e8; % Mars
v_b2 = sqrt(mu / r_b2); % Velocity of Mars

% [tau, angle_turn] = calc_synodic(mu, r_b1, r_b2);
tau = tau * sec2year;

AR_min = 0.9;
TR_min = 0.85;
p_max = 4;

h_max = 5*p_max;
s_max = 3*p_max;

r0 = [r_b1; 0; 0];
r0_2 = [r_b2; 0; 0];

mean_motion = sqrt(mu / r_b1^3);

v0 = v_b1 * [0; 1; 0];
v0_2 = v_b2 * [0; 1; 0];