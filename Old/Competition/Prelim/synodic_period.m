clear;
[consts, mu_altaira, planet_data] = load_problem_data();
clc;

num_planets = height(planet_data);

periods = zeros(1, num_planets);
rows = cell(num_planets, num_planets + 1);
vars = cell(1, num_planets + 1);

vars{1} = 'PLANET';

for i=1:num_planets
    p1 = orbital_period(mu_altaira, planet_data.a_km(i));
    rows{i, 1} = planet_data.Name{i};
    vars{i+1} = planet_data.Name{i};

    for j = 1:num_planets
        p2 = orbital_period(mu_altaira, planet_data.a_km(j));

        S = 1 / (abs(1/p1 - 1/p2));
        rows{i, j + 1} = S / consts.year_s;

        if isinf(S)
            rows{i, j + 1} = '---';
        end
    end
end

FormattedTable.Display(vars, rows)


function T = orbital_period(mu, a)
    T = 2*pi * sqrt(a^3 / mu);
end