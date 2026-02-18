function feasible_solutions = get_feasible_solutions(solutions, allow_all)

if nargin < 2
    allow_all = false;
end

feasible_solutions = {};
feasible_idx = 1;

for i = 1:length(solutions)

    solution = solutions{i};
    solution_data = solution.solution_data;

    solution_names = fieldnames(solution_data);

    for n = 1:length(solution_names)
        this_solution = solution_data.(solution_names{n});

        fast_solution = this_solution.fast;

        if isfield(fast_solution, "feasible")
            if fast_solution.feasible(1) || allow_all
                fast_solution.p = solution.p;
                fast_solution.h = solution.h;
                fast_solution.s = solution.s;
                fast_solution.tof = solution.tof;
                fast_solution.direction = '-';
                feasible_solutions{feasible_idx} = fast_solution;
                feasible_idx = feasible_idx + 1;
            end
        end

        slow_solution = this_solution.slow;

        if isfield(slow_solution, "feasible")
            if slow_solution.feasible(1) || allow_all
                slow_solution.p = solution.p;
                slow_solution.h = solution.h;
                slow_solution.s = solution.s;
                slow_solution.tof = solution.tof;
                slow_solution.direction = '+';
                feasible_solutions{feasible_idx} = slow_solution;
                feasible_idx = feasible_idx + 1;
            end
        end
    end
end
end