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

    variable_names = {'Cycler', 'AR', 'TR', 'V_inf'};

    for i = 1:num_feasible
        sol = solutions{i};
        if isempty(sol)
            continue
        end
        rows{i, 1} = sprintf('%d.%d.%d.%s%d', sol.p, sol.h, sol.s, sol.direction, sol.i);
        rows{i, 2}   = sol.AR;
        rows{i, 3}   = sol.TR;
        rows{i, 4} = sol.v_inf;
    end

    FormattedTable.Display(variable_names, rows);


end