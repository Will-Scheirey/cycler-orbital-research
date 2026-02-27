function t = calc_tof(a, mu, e, dtheta)
    T = (2 * pi) * sqrt(a^3/mu);

    % Find E, Me, and t
    E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(dtheta/2));
    Me = E - e * sin(E);
    t = T / (2 * pi) * Me;
end