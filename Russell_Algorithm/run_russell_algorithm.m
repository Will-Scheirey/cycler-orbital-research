clear; clc; close all

load_params
all_sol_idx = 1;

for p = 1:p_max
    for h = 0:h_max
        for s = 1:s_max
            
            [N_max, tof, rev_all] = russell_algorithm(tau, p, h, s, mean_motion, r_b1, v_b1, r_b2, r0, v0, mu_b1, mu, rp_min, TR_min, AR_min, year2sec);

            if isempty(rev_all)
                continue
            end

            data = struct('p', p, 'h', h, 's', s, 'N_max', N_max, 'tof', tof, 'rev_data', rev_all);
            solutions{all_sol_idx} = data;
            all_sol_idx = all_sol_idx + 1;
        end
    end
end

%% Processing

feasible_solutions = get_feasible_solutions(solutions);
print_solutions(feasible_solutions)