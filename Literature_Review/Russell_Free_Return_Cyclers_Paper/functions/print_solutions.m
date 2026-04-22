function print_solutions(solutions)
num_feasible = length(solutions);

if num_feasible == 1
    sol = solutions;
    fprintf('Feasible solution for %d.%d.%d.%s%d', sol.p, sol.h, sol.s, sol.direction, sol.i);
    fprintf('\t AR: %0.2f; TR: %0.2f; Earth V_inf: %0.1f\n', sol.AR, sol.TR, sol.v_inf);
    return
end

names = cell(1, num_feasible);
ARs = cell(1, num_feasible);
TRs = cell(1, num_feasible);
V_infs = cell(1, num_feasible);

rows = cell(num_feasible, 4);

max_num_turn_angles = 1;

for i = 1:num_feasible
    sol = solutions{i};
    max_num_turn_angles = max(max_num_turn_angles, length(horzcat(sol.deltas{:})));
end

variable_names = cell(4 + max_num_turn_angles, 1);
variable_names{1} = 'Cycler';
variable_names{2} = 'AR';
variable_names{3} = 'TR';
variable_names{4} = 'V_inf';

for n = 1:max_num_turn_angles
    variable_names{4 + n} = sprintf('TA %d', n);
end

for i = 1:num_feasible
    sol = solutions{i};
    if isempty(sol)
        continue
    end
    rows{i, 1} = sprintf('%s%d.%d.%d.%s%d', sol.orbit_direction, sol.p, sol.h, sol.s, sol.direction, sol.i);
    rows{i, 2}   = sol.AR;
    rows{i, 3}   = sol.TR;
    rows{i, 4} = sol.v_inf;
    deltas = horzcat(sol.deltas{:});

    for n = 1:max_num_turn_angles
        if n <= length(deltas)
            delta = round(rad2deg(deltas(n)));
        else
            delta = '--';
        end
        rows{i, 4+n} = delta;
    end

end

FormattedTable.Display(variable_names, rows);


end