function [dv1, dv2, a_t, tof] = calc_hohmann(mu, r1,r2)

    a_t = (r1 + r2) / 2;
    dv1 = sqrt(mu/r1) * (    sqrt((2*r2) / (r1 + r2)) - 1);
    dv2 = sqrt(mu/r2) * (1 - sqrt((2*r1) / (r1 + r2))    );
    tof = pi * sqrt(a_t^3 / mu);
end