function print_solutions(solutions)
    num_feasible = length(solutions);
    
    for i = 1:num_feasible
        sol = solutions{i};
        fprintf("Feasible solution for %d.%d.%d.%s%d", sol.p, sol.h, sol.s, sol.direction, sol.i);
        fprintf("\t AR: %0.2f; TR: %0.2f; Earth V_inf: %0.1f\n", sol.AR, sol.TR, sol.v_inf);
    end
end