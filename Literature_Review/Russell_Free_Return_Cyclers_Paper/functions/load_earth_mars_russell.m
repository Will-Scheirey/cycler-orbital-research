function [mu, mu_b1, mu_b2, r_b1, r_b2, min_rp] = load_earth_mars_russell()

[mu, mu_b1, mu_b2, r_b1, r_b2] = get_ideal_earth_mars_params();
min_rp = 6378;

%{
mu_b1_temp = mu_b1;
r_b1_temp  = r_b1;

mu_b1 = mu_b2;
r_b1  = r_b2;

mu_b2 = mu_b1_temp;
r_b2  = r_b1_temp;
%}
end

