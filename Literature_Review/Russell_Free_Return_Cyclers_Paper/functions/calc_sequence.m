function [v_inf_minus_all, deltas_all, hi, dt_years_all, C, v_inf_minus_local_all, theta_earth_all] = calc_sequence(vd, va, h, s, v_e_dep, v_e_arr, theta_generic)

v_inf_minus_all        = [];
v_inf_minus_local_all  = [];
deltas_all             = [];
hi                     = [];
dt_years_all           = [];
C                      = [];
theta_earth_all        = [];

[phi_fr, phi_gr] = calc_phi(vd, va, v_e_dep, v_e_arr);

if ~isfinite(phi_fr) || ~isfinite(phi_gr)
    return
end

v_inf_1_minus = va - v_e_arr;
v_inf_mag = norm(v_inf_1_minus);
C = calc_C(v_inf_1_minus, v_e_arr);

if v_inf_mag == 0
    return
end

hi = calc_hi(phi_gr, phi_fr, h, s);
num_hi = length(hi);

v_inf_minus_all       = cell(num_hi, 1);
v_inf_minus_local_all = cell(num_hi, 1);
deltas_all            = cell(num_hi, 1);
dt_years_all          = cell(num_hi, 1);
theta_earth_all       = cell(num_hi, 1);

theta_earth = 0;

for j = 1:num_hi

    hj = hi(j);

    [fj, dt_years] = russell_table2_fcn(hj);
    [delta_c, delta_a, delta_minimax, delta_min, lambda_a, lambda_b, lambda] = calc_deltas(phi_gr, phi_fr, fj);

    v_inf_minus_j       = cell(fj + 1, 1);
    v_inf_minus_local_j = cell(fj + 1, 1);
    theta_earth_j       = zeros(fj + 1, 1);

    if fj == 1
        deltas = delta_c;

        [vx1, vy1, vz1] = sph2cart(0, phi_gr, v_inf_mag);
        vd0 = [vx1; vy1; vz1];

        theta_earth_j(1) = theta_earth;
        theta_earth = theta_earth + theta_generic;
        theta_earth_j(2) = theta_earth;

        v_inf_minus_local_j{1} = vd0;
        v_inf_minus_j{1} = C * vd0;

        vd_restart_local = [-vd0(1); vd0(2); vd0(3)];
        v_inf_minus_local_j{2} = vd_restart_local;
        v_inf_minus_j{2} = C * vd_restart_local;

        v_inf_minus_all{j}       = v_inf_minus_j;
        v_inf_minus_local_all{j} = v_inf_minus_local_j;
        deltas_all{j}            = deltas;
        dt_years_all{j}          = dt_years;
        theta_earth_all{j}       = theta_earth_j;
        continue
    end

    if fj == 2
        [vx1, vy1, vz1] = sph2cart(0,    phi_gr, v_inf_mag);
        [vx2, vy2, vz2] = sph2cart(pi/2, phi_fr, v_inf_mag);

        vd1 = [vx1; vy1; vz1];
        vd2 = [vx2; vy2; vz2];
        vd3 = -vd2;

        theta_earth_j(1) = theta_earth;
        theta_earth = theta_earth + theta_generic;
        theta_earth_j(2) = theta_earth;
        theta_earth = theta_earth + dt_years(1) * 2*pi;
        theta_earth_j(3) = theta_earth;

        v_inf_minus_local_j{1} = vd1;
        v_inf_minus_local_j{2} = vd2;
        v_inf_minus_local_j{3} = vd3;

        v_inf_minus_j{1} = C * vd1;
        v_inf_minus_j{2} = C * vd2;
        v_inf_minus_j{3} = C * vd3;

        deltas = [delta_minimax, delta_minimax];

        v_inf_minus_all{j}       = v_inf_minus_j;
        v_inf_minus_local_all{j} = v_inf_minus_local_j;
        deltas_all{j}            = deltas;
        dt_years_all{j}          = dt_years;
        theta_earth_all{j}       = theta_earth_j;
        continue
    end

    deltas = zeros(fj - 1, 1);

    for k = 1:fj
        theta_earth_j(k) = theta_earth;

        if k == 1
            lat = phi_gr;
            lon = 0;
        else
            lat = phi_fr;
            if delta_min > delta_a
                lon = lambda_a * (k - 2);
            else
                lon = lambda + lambda_b * (k - 2);
            end
        end

        [vx, vy, vz] = sph2cart(-lon, lat, v_inf_mag);
        vd_local = [vx; vy; vz];

        if k == 1
            vd0 = vd_local;
            theta_earth = theta_earth + theta_generic;
        else
            theta_earth = theta_earth + dt_years(k - 1) * 2*pi;
        end

        v_inf_minus_local_j{k} = vd_local;
        v_inf_minus_j{k} = C * vd_local;
    end

    theta_earth_j(fj + 1) = theta_earth;

    vd_restart_local = [-vd0(1); vd0(2); vd0(3)];
    v_inf_minus_local_j{fj + 1} = vd_restart_local;
    v_inf_minus_j{fj + 1} = C * vd_restart_local;

    v_inf_minus_all{j}       = v_inf_minus_j;
    v_inf_minus_local_all{j} = v_inf_minus_local_j;
    deltas_all{j}            = deltas;
    dt_years_all{j}          = dt_years;
    theta_earth_all{j}       = theta_earth_j;
end

end

function C = calc_C(v_local, v_e)
z_hat = v_e / norm(v_e);

x_vec = v_local - dot(v_local, z_hat)*z_hat;
x_hat = x_vec / norm(x_vec);

y_hat = cross(z_hat, x_hat);
y_hat = y_hat / norm(y_hat);

C = [x_hat, y_hat, z_hat];          % local->global
end

function hi = calc_hi(phi_gr, phi_fr, h, s)

hj_temp = floor(h / s);

[fj, ~] = russell_table2_fcn(hj_temp);

[delta_c, ~, delta_minimax, ~, ~, ~, ~] = calc_deltas(phi_gr, phi_fr, fj);

hi = zeros(s, 1);

if delta_c > delta_minimax
    hi(1:s-1) = floor(h/s);
    hi(s)     = floor(h/s) + mod(h, s);
else
    hi(1)   = h;
    hi(2:s) = 0;
end

end

function [delta_c, delta_a, delta_minimax, delta_min, lambda_a, lambda_b, lambda] = calc_deltas(phi_gr, phi_fr, fj)

lambda_a  = pi / (fj - 2);
lambda_b  = 0;
lambda    = 0;

delta_a   = 0;
delta_min = 0;
delta_c   = pi - 2 * abs(phi_gr);

if fj == 1
    delta_minimax = delta_c;
elseif fj == 2
    delta_minimax = acos(sin(phi_gr) * sin(phi_fr));
elseif fj > 2
    delta_min = acos(cos(phi_fr)*cos(phi_gr) + sin(phi_fr)*sin(phi_gr));
    delta_a   = acos(cos(phi_fr)^2 * cos(lambda_a) + sin(phi_fr)^2);

    if delta_min >= delta_a
        delta_minimax = delta_min;
    else
        lb  = @(l) (pi - 2*l) / (fj - 2);
        fcn = @(l) cos(phi_fr)^2*cos(l)*cos(lb(l) + l) - cos(phi_fr)*cos(l)*cos(phi_gr) + ...
            cos(phi_fr)^2*sin(l)*sin(lb(l) + l) + sin(phi_fr)^2 - sin(phi_fr)*sin(phi_gr);
        l   = fsolve(fcn, 0, optimset('Display', 'off'));

        delta_minimax = acos(cos(phi_gr)*cos(phi_fr)*cos(l) + sin(phi_gr)*sin(phi_fr));

        lambda   = l;
        lambda_b = lb(l);
    end
end
end