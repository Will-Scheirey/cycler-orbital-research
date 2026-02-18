function sol = get_solution(p, h, s, i, direction, solutions)
    sol = [];

    num_sol = numel(solutions);

    for k = 1:num_sol
        sol_temp = solutions{k};
        if sol_temp.p == p && sol_temp.h == h && sol_temp.s == s && sol_temp.i == i && sol_temp.direction == direction
            sol = sol_temp;
            break
        end
    end
end