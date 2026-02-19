function p = plot_orbit(r0, v0, mu, display_name)
[~, ~, e, p, w, ~] = calc_orbit(r0, v0, mu);

thetas = linspace(0, 2*pi, 1000);

r = p ./ (1 + e*cos(thetas - w));

r_vec = r .* [cos(thetas); sin(thetas)];

if nargin < 4
    display_name = "Orbit";
end

p = plot(r_vec(1, :), r_vec(2, :), 'LineWidth', 1.5, 'DisplayName', display_name);
end