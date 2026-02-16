function plot_lambert_arc_time(R1, V1, mu, tof, color)
    % Two-body ODE
    f = @(t,x) [ x(4:6);
                -mu * x(1:3) / norm(x(1:3))^3 ];
    x0 = [R1; V1];

    % Sample uniformly in time (avoids anomaly decisions)
    tgrid = linspace(0, tof, 400);
    opts  = odeset('RelTol',1e-12, 'AbsTol',1e-12);
    [~, X] = ode113(f, tgrid, x0, opts);

    plot3(X(:,1), X(:,2), X(:,3), color, 'LineWidth', 2);
end