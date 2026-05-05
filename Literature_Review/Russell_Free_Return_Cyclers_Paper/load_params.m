% [mu, mu_b1, mu_b2, r_b1, r_b2, min_rp] = load_jupiter_moons_russel();
[mu, mu_b1, mu_b2, r_b1, r_b2, min_rp] = load_earth_mars_russell();

year2sec = 2*pi*sqrt(r_b1^3 / mu);
sec2year = 1 / year2sec;

% AU = 149597871; % [km]
AU = r_b1;

% r_b1 = AU; % Earth
v_b1 = sqrt(mu / r_b1); % Velocity of Earth
% mu_b1 = 3.986004418e5;
rp_min = min_rp;

% r_b2 = 2.2737e8; % Mars
v_b2 = sqrt(mu / r_b2); % Velocity of Mars

[tau, angle_turn] = calc_synodic(mu, r_b1, r_b2);
angle_turn = -angle_turn;
tau = tau * sec2year;

AR_min = 0.9;
TR_min = 0.85;

p_min = 1;
p_max = 15;

h_min = 0;
s_min = 1;

h_max = 5*p_max;
s_max = 5*p_max;

r0 = [r_b1; 0; 0];
r0_2 = [r_b2; 0; 0];

mean_motion = sqrt(mu / r_b1^3);
mean_motion_2 = sqrt(mu / r_b2^3);

v0 = v_b1 * [0; 1; 0];
v0_2 = v_b2 * [0; 1; 0];