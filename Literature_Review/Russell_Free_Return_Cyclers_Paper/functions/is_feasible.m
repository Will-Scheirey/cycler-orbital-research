function [AR, TR, feasible, delta_max] = is_feasible(a, rp_min, mu, r_b2, v_inf_mag, deltas, TR_min, AR_min, ra)

AR = ra / r_b2;

e = 1 + (rp_min * v_inf_mag^2) / mu;
delta_max = 2*asin(1/e);

max_delta = -inf;
for i = 1:length(deltas)
    max_delta = max(max_delta, max(deltas{i}));
end

TR = delta_max / max_delta;

feasible = TR > TR_min && AR > AR_min;
end