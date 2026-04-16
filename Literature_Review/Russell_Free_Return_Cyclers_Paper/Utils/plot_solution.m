function [orbits, pos_all] = plot_solution(solution, mu, vd_generic, r_b1, v_b1, theta_generic, x_offset)

orbits = get_orbits(solution, mu, vd_generic, r_b1, v_b1, theta_generic);
num_orbits = numel(orbits);
colors = turbo(num_orbits);
lw = 3;
%{
    plot_orbit(r0, vd_generic, mu, 'DisplayName', 'Initial Generic', 'LineWidth', 1.5, 'Color', color);
    quiver(0, 0, r0(1), r0(2), 0, 'LineWidth', 1.5, 'Color', color, 'HandleVisibility', 'off');
    quiver(r0(1), r0(2), vd_generic(1), vd_generic(2), 5e6, 'LineWidth', 1.5, 'Color', color, 'HandleVisibility', 'off');
%}

pos_all = cell(num_orbits, 1);

for n = 1:num_orbits
    orbit = orbits{n};

    r0  = orbit{1};
    vd  = orbit{2};
    tof = orbit{3};
    ve  = orbit{4};
    
    figure(1)
    plot_orbit_time(r0, vd, mu, [0, tof], 'DisplayName', sprintf('Trajectory %d', n), 'LineWidth', lw, 'Color', colors(n, :), 'Clipping', 'off');
    
    quiver3(r0(1), r0(2), r0(3), vd(1), vd(2), vd(3), 5e6, 'HandleVisibility', 'off', 'Color', [colors(n, :), 0.2], 'LineWidth', lw);
    quiver3(r0(1), r0(2), r0(3), ve(1), ve(2), ve(3), 5e6, 'HandleVisibility', 'off', 'Color', 'b', 'LineWidth', lw);

    if ~isempty(solution.theta_earth_all)
        theta_earth = solution.theta_earth_all(n);
        x = r_b1 * cos(theta_earth);
        y = r_b1 * sin(theta_earth);

        plot3(x, y, 0, 'x', 'Color', colors(n, :), 'MarkerSize', 30, 'LineWidth', 5, 'DisplayName', sprintf('Earth Start %d', n))
    end

    pos_all{n} = {r0, vd, tof};
end

r0_locs = {orbits{2}{1}};
r0_nums = {{'1'}};
r0_idx = 2;

for n = 3:num_orbits
    r0  = orbits{n}{1};
    match_idx = 0;
    for z = 1:r0_idx-1
        r0_comp = r0_locs{z};
        diff = norm(r0_comp - r0);
        if diff < 1e4
            match_idx = z;
            break
        end
    end

    if match_idx
        r0_nums{match_idx} = vertcat(r0_nums{match_idx}, sprintf('%d', n));
    else
        r0_locs = vertcat(r0_locs, r0);
        r0_nums{r0_idx} = {sprintf('%d', n)};
        r0_idx = r0_idx + 1;
    end
end

num_intersect = length(r0_nums);
for n = 1:num_intersect
    r0  = r0_locs{n};
    thestr = strjoin(r0_nums{n}, ', ');
    plot_intersection(r0, x_offset, sprintf('%s', thestr));
end

legend;
axis equal
end