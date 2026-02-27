clc; clear;

AU = 1.496e+8;
TU = 5022642.89376;
mu = 132712000000;

r1 = 1;
r2 = 1;

num_N = 9;

num_a = 1000;
a = linspace(0.89, 5, num_a);

figure(1)
clf

theta = @(N) 2*N*pi + 4*pi/7;

lw = 0.8;

for N=0:num_N

    tof_eqn_fast = lambert_lagrange_eqn(r1, r2, theta(N), 1, N, true);
    tof_eqn_slow = lambert_lagrange_eqn(r1, r2, theta(N), 1, N, false);

    plot_fast = real(tof_eqn_fast(a));
    plot_slow = real(tof_eqn_slow(a));

    if N == 0
        semilogx(a, plot_fast, '-k', 'DisplayName', 'Fast Transfers (\alpha=\alpha_0)', 'LineWidth', lw); hold on
        semilogx(a, plot_slow, '--k', 'DisplayName', 'Slow Transfers (\alpha=2\pi-\alpha_0)', 'LineWidth', lw);
    else
        semilogx(a, plot_fast, '-k', 'HandleVisibility', 'off', 'LineWidth', lw);
        semilogx(a, plot_slow, '--k', 'HandleVisibility', 'off', 'LineWidth', lw);
    end

    text(a(1) - 2e-2, plot_fast(1), sprintf("N=%d", N), 'HorizontalAlignment', 'right', 'FontSize', 12)
end

xlabel("Semi-major Axis (AU)")
ylabel("Time of Flight (TU)")

ylim([0, 55])
xlim([0.7, 4.2])

legend('Location', 'best')

function tof_eqn = lambert_lagrange_eqn(r1, r2, theta, mu, N, fast)

    c  = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(theta));
    S  = 0.5*(c + r1 + r2);

    alpha0 = @(a) 2 * asin(sqrt(S ./ (2*a)));
    beta0  = @(a) 2 * asin(sqrt((S - c) ./ (2*a)));

    if fast
        alpha = @(a) alpha0(a);
    else
        alpha = @(a) 2*pi - alpha0(a);
    end

    if theta - 2*pi*N < pi
        beta = @(a) beta0(a);
    else
        beta = @(a) -beta0(a);
    end

    tof_eqn = @(a) sqrt(a.^3 / mu) .* (2*N*pi + alpha(a) - beta(a) - sin(alpha(a)) + sin(beta(a)));
end