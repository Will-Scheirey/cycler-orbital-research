clear; clc; close all

disp("--- Russell Table 2 Test ---")
for hj = 0:8
    [fj, dt_years] = russell_table2_fcn(hj);

    if isempty(dt_years)
        year_str = '--';
    else
        year_str = join(string(dt_years), ', ');
    end

    fprintf("hj: %d | fj: %d | dt: %s\n", hj, fj, year_str);

end