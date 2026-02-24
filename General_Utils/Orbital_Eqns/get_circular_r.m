function r_vec = get_circular_r(mu, r, t)
    n = mean_motion(r, mu);

    theta = n * t;

    r_vec = r * [cos(theta); sin(theta); 0];
end