function [num_prograde, num_retrograde] = print_solutions_latex(solutions)

num_feasible = length(solutions);

if num_feasible == 1 && ~iscell(solutions)
    solutions = {solutions};
    num_feasible = 1;
end

max_num_turn_angles = 1;

for i = 1:num_feasible
    sol = solutions{i};
    if isempty(sol)
        continue
    end
    max_num_turn_angles = max(max_num_turn_angles, length(horzcat(sol.deltas{:})));
end

max_num_turn_angles = min(max_num_turn_angles, 10);

headers = cell(1, 4 + max_num_turn_angles);
headers{1} = 'Cycler';
headers{2} = 'AR';
headers{3} = 'TR';
headers{4} = '$V_{\infty}$';

for n = 1:max_num_turn_angles
    headers{4+n} = sprintf('TA %d', n);
end

fprintf('\\begin{table}[h!]\n');
fprintf('\\centering\n');
fprintf('\\begin{tabular}{%s}\n', repmat('c ', 1, length(headers)));
fprintf('\\hline\n');

for j = 1:length(headers)
    if j < length(headers)
        fprintf('%s & ', headers{j});
    else
        fprintf('%s \\\\\n', headers{j});
    end
end

fprintf('\\hline\n');

num_prograde = 0;
num_retrograde = 1;

for i = 1:num_feasible
    sol = solutions{i};

    if isempty(sol)
        continue
    end

    if strcmp(sol.orbit_direction, '-')
        num_retrograde = num_retrograde + 1;
    else
        num_prograde = num_prograde + 1;
    end


    cycler = sprintf('%s%d.%d.%d.%s%d', ...
        sol.orbit_direction, sol.p, sol.h, sol.s, sol.direction, sol.i);

    deltas = horzcat(sol.deltas{:});

    row = cell(1, length(headers));
    row{1} = cycler;
    row{2} = sprintf('%.3f', sol.AR);
    row{3} = sprintf('%.3f', sol.TR);
    row{4} = sprintf('%.3f', sol.v_inf);

    for n = 1:max_num_turn_angles
        if n <= length(deltas)
            row{4+n} = sprintf('%d', round(rad2deg(deltas(n))));
        else
            row{4+n} = '--';
        end
    end

    for j = 1:length(row)
        if j < length(row)
            fprintf('%s & ', row{j});
        else
            fprintf('%s \\\\\n', row{j});
        end
    end
end

fprintf('\\hline\n');
fprintf('\\end{tabular}\n');
fprintf('\\caption{Cycler Data Table}\n');
fprintf('\\end{table}\n');

end