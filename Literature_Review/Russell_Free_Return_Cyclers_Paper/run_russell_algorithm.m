% Based on: https://doi.org/10.2514/1.1011

clear; clc; close all
load_params

% Build all (p,h,s) combinations first

num_cases = (p_max - p_min + 1) * (h_max - h_min + 1) * (s_max - s_min + 1);
case_list = zeros(num_cases, 3);

idx = 1;

fprintf('Generating cases...\n')
for p = p_min:p_max
    for h = h_min:h_max
        for s = s_min:s_max
            case_list(idx, :) = [p, h, s];
            idx = idx + 1;
        end
    end
end

solutions = cell(num_cases, 1);
fprintf('Running cases (%d)...\n', num_cases)

for idx = 1:num_cases
    p = case_list(idx, 1);
    h = case_list(idx, 2);
    s = case_list(idx, 3);

    [N_max, tof, rev_all] = russell_algorithm( ...
    tau, p, h, s, ...
    mean_motion, r_b1, v_b1, r_b2, ...
    r0, v0, mu_b1, mu, ...
    rp_min, TR_min, AR_min, year2sec);

    if isempty(rev_all)
        solutions{idx} = [];
    else
        solutions{idx} = struct( ...
        'p', p, ...
        'h', h, ...
        's', s, ...
        'N_max', N_max, ...
        'tof', tof, ...
        'rev_data', rev_all);
    end

    fprintf('Finished %d/%d\n', idx, num_cases)
end

clc

solutions = solutions(~cellfun(@isempty, solutions));
%% Processing

if isempty(solutions)
    disp("No solutions!")
    return
end

feasible_solutions = get_feasible_solutions(solutions, "retrograde_only", false);

if isempty(feasible_solutions)
    disp("No feasible solutions!")
    return
end

%% PRINT
clc
% print_solutions_latex(feasible_solutions)
diary('cycler_table.txt'); 
[num_prograde, num_retrograde] = print_solutions_latex(feasible_solutions); 
diary off