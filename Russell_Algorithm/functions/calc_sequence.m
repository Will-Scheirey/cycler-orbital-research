function [v_inf_minus, deltas, v_inf_local] = calc_sequence(vd, va, h, s, v_e_dep, v_e_arr)

v_inf_minus = [];
v_inf_local = [];
deltas = [];

[phi_fr, phi_gr] = calc_phi(vd, va, v_e_dep, v_e_arr);

if ~isfinite(phi_fr) || ~isfinite(phi_gr)
    return
end

delta_c = pi - 2 * abs(phi_gr);

hj = floor(h/s);

if mod(hj, 2) == 0
    fj = hj/2 + 1;
else
    fj = 2 * floor(hj/4 + 1);
end

lambda_a = pi / (fj - 2);

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
        lb = @(l) (pi - 2*l) / (fj - 2);
        fcn = @(l) cos(phi_fr)^2*cos(l)*cos(lb(l) + l) - cos(phi_fr)*cos(l)*cos(phi_gr) + ...
            cos(phi_fr)^2*sin(l)*sin(lb(l) + l) + sin(phi_fr)^2 - sin(phi_fr)*sin(phi_gr);
        l = fsolve(fcn, 0, optimset('Display', 'off'));

        delta_minimax = acos(cos(phi_gr)*cos(phi_fr)*cos(l) + sin(phi_gr)*sin(phi_fr));
    end
end

h_all = zeros(s, 1);

if delta_c > delta_minimax
    h_all(1:s-1) = floor(h/s);
    h_all(s)   = floor(h/s) + mod(h, s);
else
    h_all(1) = h;
    h_all(2:s) = 0;
end

% Assemble reference frame

v_inf_1_minus = va - v_e_arr;
v_inf_mag = norm(v_inf_1_minus);

if v_inf_mag == 0
    return
end

z_hat = v_e_arr / norm(v_e_arr);

x_vec = v_inf_1_minus - dot(v_inf_1_minus, z_hat)*z_hat;
x_hat = x_vec / norm(x_vec);

y_hat = cross(z_hat, x_hat);
y_hat = y_hat / norm(y_hat);

C = [x_hat, y_hat, z_hat];          % local->global

if fj == 1
    deltas = delta_c;

    [vx1, vy1, vz1] = sph2cart(0, phi_gr, v_inf_mag);
    v_inf_minus = {[vx1; vy1; vz1]};
    v_global = v_inf_minus{1};
    % v_global = C * [vx1;vy1;vz1];

    v_inf_local = {v_global, v_global};
    v_inf_local{2}(1) = -v_inf_local{2}(1);
    return
end

if fj == 2
    [vx1, vy1, vz1] = sph2cart(-0,      phi_gr, v_inf_mag);
    [vx2, vy2, vz2] = sph2cart(-pi/2,   phi_gr, v_inf_mag);

    v_inf_minus = cell(2,1);
    v_inf_minus{1} = C * [vx1; vy1; vz1];
    v_inf_minus{2} = C * [vx2; vy2; vz2];

    deltas = [delta_minimax, delta_minimax];
    return
end

v_inf_local = cell(fj+1, 1);

v_inf_minus = cell(fj, 1);
deltas = zeros(fj - 1, 1);

for k=1:fj
    if k == 1
        lat = phi_gr;
        lon = 0;
    else
        lat = phi_fr;
        if delta_min > delta_a
            lon = lambda_a * (k - 2);
        else
            lon = l + lb(l) * (k-2);
        end
    end

    [vx, vy, vz] = sph2cart(-lon, lat, v_inf_mag);

    v_global = C * [vx;vy;vz];

    v_inf_local{k} = [vx; vy; vz];
    v_inf_minus{k} = v_global;

    if k >= 2
        v_prev = v_inf_minus{k-1};
        v_curr = v_inf_minus{k};

        deltas(k-1) = acos(dot(v_prev, v_curr) / (norm(v_prev)*norm(v_curr)));
    end
end

v_inf_local{fj+1} = v_inf_local{1};
v_inf_local{fj+1}(1) = -v_inf_local{fj+1}(1);

end