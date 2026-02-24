function [r_out, h_out, e_out, p_out, w_out, nu_out, theta_intersect, rp_sol, info] = ...
    calc_grav_assist_direct(r0x, r0y, v0x, v0y, mu_sun, mu_b2, r_b2, ...
                            theta_E_target, sign_turn, which_intersection, opts)
%CALC_GRAV_ASSIST_DIRECT
%   Computes a gravity assist at body-2 circular orbit radius r_b2 and SOLVES rp
%   such that the OUTGOING heliocentric orbit intersects r = r0x (Earth orbit radius)
%   at inertial polar angle theta_E_target.
%
%   This avoids any derived scalar term that can admit spurious roots.
%
% Inputs:
%   r0x, r0y, v0x, v0y   spacecraft heliocentric state at start (km, km/s)
%   mu_sun               Sun GM (km^3/s^2)
%   mu_b2                Body-2 GM (km^3/s^2) (Jupiter/Mars etc.)
%   r_b2                 Body-2 orbital radius (km) (assumed circular, coplanar)
%   theta_E_target       inertial polar angle where you want r = r0x on outbound orbit
%   sign_turn            +1 or -1 for flyby bend direction (choose one)
%   which_intersection   1 => use theta_int = w_sc + nu_int, 2 => w_sc - nu_int
%   opts                 struct:
%       rp_min, rp_max   bounds to search (km)
%       n_grid           log grid samples for bracketing
%       verbose          true/false
%
% Outputs:
%   ... usual outgoing orbit elements plus:
%   rp_sol               solved periapsis radius (km)
%   info                 diagnostics (status, bracket, samples)

    if nargin < 11 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'rp_min'),   opts.rp_min = 1e2; end
    if ~isfield(opts,'rp_max'),   opts.rp_max = 1e7; end
    if ~isfield(opts,'n_grid'),   opts.n_grid = 120; end
    if ~isfield(opts,'verbose'),  opts.verbose = false; end

    % ---- Spacecraft initial orbit about the Sun ----
    [~, h_sc, e_sc, p_sc, w_sc, ~] = calc_orbit([r0x; r0y; 0], [v0x; v0y; 0], mu_sun);

    % ---- Find intersection with body-2 circular orbit (r = r_b2) ----
    c = (p_sc/r_b2 - 1)/e_sc;
    if ~isfinite(c) || abs(c) > 1
        % no intersection
        [r_out,h_out,e_out,p_out,w_out,nu_out,theta_intersect,rp_sol,info] = deal([]);
        return
    end

    nu_int = acos(max(-1,min(1,c)));

    phi_a = mod(w_sc + nu_int, 2*pi);
    phi_b = mod(w_sc - nu_int, 2*pi);

    if which_intersection == 1
        theta_intersect = phi_a;
    else
        theta_intersect = phi_b;
    end

    % true anomaly of SC at encounter
    nu_sc_enc = wrapToPi(theta_intersect - w_sc);

    % SC heliocentric RT at encounter
    vr_sc = mu_sun/h_sc * e_sc * sin(nu_sc_enc);
    vt_sc = mu_sun/h_sc * (1 + e_sc*cos(nu_sc_enc));

    % Body-2 heliocentric RT at encounter (circular orbit)
    vr_b2 = 0;
    vt_b2 = sqrt(mu_sun / r_b2);

    % Planet-relative v-infinity in RT frame
    v_inf_in = [vr_sc - vr_b2; vt_sc - vt_b2];
    v_inf    = norm(v_inf_in);

    % ---- Solve rp by DIRECT intersection residual ----
    [rp_sol, info] = solve_rp_by_direct_intersection( ...
        mu_sun, mu_b2, r0x, r_b2, theta_intersect, theta_E_target, ...
        v_inf_in, vr_b2, vt_b2, sign_turn, opts);

    if ~isfinite(rp_sol)
        [r_out,h_out,e_out,p_out,w_out,nu_out] = deal([]);
        return
    end

    % ---- With rp_sol, compute outbound orbit elements for output ----
    [r_vec_out, v_vec_out] = outgoing_state_from_rp( ...
        mu_sun, mu_b2, r_b2, theta_intersect, v_inf_in, vr_b2, vt_b2, rp_sol, sign_turn);

    [r_out, h_out, e_out, p_out, w_out, nu_out] = calc_orbit(r_vec_out(:), v_vec_out(:), mu_sun);

end


function [rp_sol, info] = solve_rp_by_direct_intersection( ...
    mu_sun, mu_b2, r_earth, r_b2, theta_enc, theta_E_target, ...
    v_inf_in, vr_b2, vt_b2, sign_turn, opts)
% Solves f(rp)=0 where f(rp) = r_out(theta_E_target) - r_earth

    v_inf = norm(v_inf_in);

    % Residual scaling: make O(1)
    % r residual scaled by r_earth is very well-conditioned.
    function y = f_of_rp(rp)
        if ~isfinite(rp) || rp <= 0
            y = NaN; return
        end

        % hyperbola e and turn angle
        e_hyp = 1 + rp * v_inf^2 / mu_b2;
        if e_hyp <= 1 || ~isfinite(e_hyp)
            y = NaN; return
        end

        delta = sign_turn * 2 * asin(1/e_hyp);

        % rotate v_inf_in by delta in RT plane
        R = [cos(delta) -sin(delta); sin(delta) cos(delta)];
        v_inf_out = R * v_inf_in;

        % heliocentric RT out
        vr_out = vr_b2 + v_inf_out(1);
        vt_out = vt_b2 + v_inf_out(2);

        % build heliocentric inertial state at encounter
        rx = r_b2*cos(theta_enc);
        ry = r_b2*sin(theta_enc);

        vx = vr_out*cos(theta_enc) - vt_out*sin(theta_enc);
        vy = vr_out*sin(theta_enc) + vt_out*cos(theta_enc);

        [~, ~, e_out, p_out, w_out, ~] = calc_orbit([rx; ry; 0], [vx; vy; 0], mu_sun);

        % If outbound is parabolic-ish or invalid, reject
        if ~isfinite(e_out) || ~isfinite(p_out)
            y = NaN; return
        end

        nuE = wrapToPi(theta_E_target - w_out);
        r_at_thetaE = p_out / (1 + e_out*cos(nuE));

        if ~isfinite(r_at_thetaE)
            y = NaN; return
        end

        y = (r_at_thetaE - r_earth) / r_earth;  % scaled residual
    end

    % ---- Bracket search in logspace ----
    r_grid = logspace(log10(opts.rp_min), log10(opts.rp_max), opts.n_grid);
    f_grid = arrayfun(@f_of_rp, r_grid);

    good = isfinite(f_grid);
    rg = r_grid(good);
    fg = f_grid(good);

    info = struct();
    info.r_grid = r_grid;
    info.f_grid = f_grid;
    info.good_mask = good;
    info.status = "no_sign_change";
    info.bracket = [NaN NaN];
    info.best_rp = NaN;
    info.best_abs_resid = NaN;

    rp_sol = NaN;

    if isempty(fg) || numel(fg) < 2
        info.status = "no_valid_samples";
        if opts.verbose
            fprintf("[solve_rp] no valid samples in [%.3g, %.3g] km\n", opts.rp_min, opts.rp_max);
        end
        return
    end

    [bestAbs, kbest] = min(abs(fg));
    info.best_rp = rg(kbest);
    info.best_abs_resid = bestAbs;

    idx = find(sign(fg(1:end-1)) ~= sign(fg(2:end)), 1, 'first');
    if isempty(idx)
        if opts.verbose
            fprintf("[solve_rp] no sign change in [%.3g, %.3g] km. best |resid|=%.3e at rp=%.3g\n", ...
                opts.rp_min, opts.rp_max, bestAbs, info.best_rp);
        end
        return
    end

    a = rg(idx);
    b = rg(idx+1);
    info.bracket = [a b];

    % ---- Solve ----
    try
        rp_sol = fzero(@f_of_rp, [a b]);
        info.status = "ok";
        info.resid_at_sol = f_of_rp(rp_sol);
        if opts.verbose
            fprintf("[solve_rp] rp=%.6g km, resid=%.3e, bracket=[%.3g, %.3g]\n", rp_sol, info.resid_at_sol, a, b);
        end
    catch ME
        info.status = "fzero_failed";
        info.error = ME.message;
        rp_sol = NaN;
        if opts.verbose
            fprintf("[solve_rp] fzero failed: %s\n", ME.message);
        end
    end
end


function [r_vec_out, v_vec_out] = outgoing_state_from_rp( ...
    mu_sun, mu_b2, r_b2, theta_enc, v_inf_in, vr_b2, vt_b2, rp, sign_turn)
% Returns heliocentric inertial state at encounter after applying flyby with rp

    v_inf = norm(v_inf_in);
    e_hyp = 1 + rp * v_inf^2 / mu_b2;
    delta = sign_turn * 2 * asin(1/e_hyp);

    R = [cos(delta) -sin(delta); sin(delta) cos(delta)];
    v_inf_out = R * v_inf_in;

    vr_out = vr_b2 + v_inf_out(1);
    vt_out = vt_b2 + v_inf_out(2);

    rx = r_b2*cos(theta_enc);
    ry = r_b2*sin(theta_enc);

    vx = vr_out*cos(theta_enc) - vt_out*sin(theta_enc);
    vy = vr_out*sin(theta_enc) + vt_out*cos(theta_enc);

    r_vec_out = [rx; ry; 0];
    v_vec_out = [vx; vy; 0];
end


function [r, h, e, p, w, nu] = calc_orbit(r_vec, v_vec, mu)
%CALC_ORBIT Basic planar orbit elements from inertial state (2D in XY-plane)
    r_vec = r_vec(:); v_vec = v_vec(:);

    h_vec = cross(r_vec, v_vec);
    e_vec = cross(v_vec, h_vec)/mu - r_vec/norm(r_vec);

    h = norm(h_vec);
    e = norm(e_vec);
    p = h^2/mu;

    % Argument of periapsis in inertial XY-plane (longitude of periapsis)
    if e < 1e-12
        w = 0; % circular fallback
    else
        w = atan2(e_vec(2), e_vec(1));
    end

    % inertial polar angle of position
    theta = atan2(r_vec(2), r_vec(1));
    nu = wrapToPi(theta - w);

    r = p/(1 + e*cos(nu));
end


function ang = wrapToPi(ang)
% Wrap angle to (-pi, pi]
    ang = mod(ang + pi, 2*pi) - pi;
end