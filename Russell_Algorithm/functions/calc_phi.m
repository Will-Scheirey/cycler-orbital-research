function [phi_fr, phi_gr] = calc_phi(vd_vec, va_vec, ve_vec_dep, ve_vec_arr)
ve_dep = norm(ve_vec_dep);
ve_arr = norm(ve_vec_arr);

v_inf_vec   = vd_vec - ve_vec_dep;   % dep v-infinity
v_inf_1_vec = va_vec - ve_vec_arr;   % arr v-infinity

v_inf   = norm(v_inf_vec);
v_inf_1 = norm(v_inf_1_vec);

q = v_inf / (2*ve_dep);
if ~isfinite(q) || q > 1
    phi_fr = NaN; phi_gr = NaN;
    return
end
phi_fr = -asin(q);

u = (v_inf_1_vec' * ve_vec_arr) / (v_inf_1 * ve_arr);
u = max(-1,min(1,u));
phi_gr = asin(u);
end