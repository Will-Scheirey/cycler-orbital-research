function [orbits, pos_all] = plot_solution(solution, mu, vd_generic, r_b1, v_b1, theta_generic, x_offset)

orbits = get_orbits(solution, mu, vd_generic, r_b1, v_b1, theta_generic);
num_orbits = numel(orbits);

colors = turbo(num_orbits);
lw = 3;

pos_all = cell(num_orbits, 1);

for n = 1:num_orbits
    orbit = orbits{n};

    r0  = orbit{1};
    vd  = orbit{2};
    tof = orbit{3};
    ve  = orbit{4};

    figure(1)
    plot_orbit_time(r0, vd, mu, [0, tof], ...
        'DisplayName', sprintf('Trajectory %d', n), ...
        'LineWidth', lw, ...
        'Color', colors(n, :), ...
        'Clipping', 'off');
    hold on

    if n == 1
        quiver3(r0(1), r0(2), r0(3), vd(1), vd(2), vd(3), 5e6, ...
            'HandleVisibility', 'off', ...
            'Color', [colors(n, :), 0.2], ...
            'LineWidth', lw);
    end
    %{
    quiver3(r0(1), r0(2), r0(3), ve(1), ve(2), ve(3), 5e6, ...
        'HandleVisibility', 'off', ...
        'Color', 'b', ...
        'LineWidth', lw);
    %}
    %{
    % Mark the Earth position where this orbit starts
    plot3(r0(1), r0(2), r0(3), 'x', ...
        'Color', colors(n, :), ...
        'MarkerSize', 18, ...
        'LineWidth', 3, ...
        'HandleVisibility', 'off');
    %}

    pos_all{n} = {r0, vd, tof};
end

% Group identical/sufficiently-close start locations so they get one label
r0_locs = {};
r0_nums = {};
r0_idx = 1;

for n = 2:num_orbits
    r0 = orbits{n}{1};

    match_idx = 0;
    for z = 1:r0_idx-1
        r0_comp = r0_locs{z};
        if norm(r0_comp - r0) < 1e4
            match_idx = z;
            break
        end
    end

    if match_idx ~= 0
        r0_nums{match_idx} = vertcat(r0_nums{match_idx}, sprintf('%d', n));
    else
        r0_locs{r0_idx} = r0;
        r0_nums{r0_idx} = {sprintf('%d', n)};
        r0_idx = r0_idx + 1;
    end
end

num_intersect = numel(r0_nums);
for n = 1:num_intersect
    r0 = r0_locs{n};
    thestr = strjoin(r0_nums{n}, ', ');
    plot_intersection(r0, x_offset, thestr);
end

legend
axis equal
set(gca, 'SortMethod', 'childorder')
end