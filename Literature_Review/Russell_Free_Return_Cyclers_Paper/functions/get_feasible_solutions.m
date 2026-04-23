function feasible_solutions = get_feasible_solutions(solutions, allow_all)

if nargin < 2
    allow_all = false;
end

feasible_solutions = {};
feasible_idx = 1;

for i = 1:length(solutions)

    solution = solutions{i};
    rev_data = solution.rev_data;

    if isempty(rev_data)
        continue;
    end

    rev_num = fieldnames(rev_data);

    for n = 1:length(rev_num)
        this_rev = rev_data.(rev_num{n});

        rev_solutions = fieldnames(this_rev);

        for z = 1:length(rev_solutions)
            this_solution = this_rev.(rev_solutions{z});

            if ~isfield(this_solution, "feasible")
                continue
            end
            if ~this_solution.feasible(1) && ~allow_all
                continue
            end
            this_solution.p = solution.p;
            this_solution.h = solution.h;
            this_solution.s = solution.s;
            this_solution.tof = solution.tof;
            if this_solution.fast
                this_solution.direction = '-';
            else
                this_solution.direction = '+';
            end
            if this_solution.long
                this_solution.orbit_direction = '-';
            else
                this_solution.orbit_direction = '+';
            end
            feasible_solutions{feasible_idx} = this_solution;
            feasible_idx = feasible_idx + 1;
        end
    end
end
end