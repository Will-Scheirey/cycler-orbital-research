function sols = get_solution(p, h, s, i, solutions, fast, long)
    sols = {};
    sol_idx = 1;

    num_sol = numel(solutions);

    for k = 1:num_sol
        sol_temp = solutions{k};
        if sol_temp.p == p && sol_temp.h == h && sol_temp.s == s && sol_temp.i == i
            if nargin >= 6
                if sol_temp.fast ~= fast
                    continue
                end
            end
            if nargin >= 7
                if sol_temp.long ~= long
                    continue
                end
            end
            sols{sol_idx} = sol_temp;
            sol_idx = sol_idx + 1;
        end
    end
end