function [mu, mu_b1, mu_b2, r_b1, r_b2, min_rp] = load_jupiter_moons_russel()
jupiter_moons = get_jupiter_moon_params();

mu = jupiter_moons.mu_primary;

mu_b1 = jupiter_moons.io.mu;
mu_b2 = jupiter_moons.europa.mu;

r_b1 = jupiter_moons.io.r_primary;
r_b2 = jupiter_moons.europa.r_primary;

min_rp = jupiter_moons.io.r;


end