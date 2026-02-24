clear;
[consts, mu_altaira, planet_data] = load_problem_data();
clc

num_planets = height(planet_data);

periods = zeros(1, num_planets);
rows = cell(11, 5);

max_year = 200;

for i=1:num_planets
    period = orbital_period(mu_altaira, planet_data.a_km(i));
    period_year = period / consts.year_s;

    max_flybys = min(floor(max_year / (period_year / 3)) + 1, 13);

    rows{i, 1} = planet_data.Name{i};
    rows{i, 2} = period_year;
    rows{i, 3} = max_flybys;
    rows{i, 4} = max_flybys * planet_data.Weight__(i);
    rows{i, 5} = planet_data.i_deg(i);
end

FormattedTable.Display({'Planet', 'Orbital Period (year)', 'Max Flybys', 'Max Solo Points', 'Inclination (deg)'}, rows)

function T = orbital_period(mu, a)
    T = 2*pi * sqrt(a^3 / mu);
end