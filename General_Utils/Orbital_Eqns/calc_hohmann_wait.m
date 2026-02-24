function t_wait = calc_hohmann_wait(r1, r2, mu)
    n1 = mean_motion(r1, mu);
    n2 = mean_motion(r2, mu);

    [~, ~, ~, tof_h] = calc_hohmann(mu, r1, r2);

    phaseReq = mod(pi - n2*tof_h, 2*pi);

    dn = n2 - n1;
    
    if dn > 0
        t_wait = phaseReq / dn;
    else
        t_wait = (2*pi - phaseReq) / (-dn);
    end
end