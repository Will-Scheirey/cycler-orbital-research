function theta = theta_from_r(p, e, w, r)
    theta = acos((p / r - 1) / e) + w;
end