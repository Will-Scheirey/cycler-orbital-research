function [rp_sol, info] = solve_rp_numeric(mu_1, mu_2, r1, theta_r, ...
                                          v_inf_in_r, v_inf_in_t, theta, ...
                                          Vr_b2, Vt_b2, sign_turn, opts)
%SOLVE_RP_NUMERIC  Numerically solve for periapsis radius rp (no symbolic).
%
%   Solves term(rp)=0 using fzero with automatic log-space bracketing.
%   Includes robust residual scaling to improve numeric conditioning.
%
%   Outputs:
%     rp_sol : solution in km (NaN if no root found)
%     info  : struct with diagnostics, including scaling used.

    if nargin < 10 || isempty(sign_turn), sign_turn = +1; end
    if nargin < 11, opts = struct(); end

    % ---------- Defaults ----------
    if ~isfield(opts,'rp_min'),   opts.rp_min = 1e2;   end   % km
    if ~isfield(opts,'rp_max'),   opts.rp_max = 1e7;   end   % km
    if ~isfield(opts,'n_grid'),   opts.n_grid = 90;    end
    if ~isfield(opts,'penalty'),  opts.penalty = 1e30; end
    if ~isfield(opts,'verbose'),  opts.verbose = false; end

    % Scaling options (new)
    % "auto" chooses a characteristic magnitude based on r1 and speeds.
    if ~isfield(opts,'scale_mode'),   opts.scale_mode = "auto"; end
    if ~isfield(opts,'scale_value'),  opts.scale_value = 1;      end % used if scale_mode="fixed"
    if ~isfield(opts,'scale_floor'),  opts.scale_floor = 1;      end % avoid divide-by-zero

    rp_min = opts.rp_min;
    rp_max = opts.rp_max;

    % Precompute
    v_inf = hypot(v_inf_in_r, v_inf_in_t);

    % ---------- Choose a characteristic scale ----------
    % The raw term has units of (km^3/s^2)^2 = km^6/s^4-ish because of r^2*v^4 terms.
    % We just want an O(1)-ish residual; any positive scale works.
    %
    % A decent scale is ~ (r1^2 * v_char^4) where v_char is some heliocentric speed scale.
    v_char = max([abs(Vr_b2) abs(Vt_b2) v_inf 1]);  % km/s
    scale_auto = (max(r1,1)^2) * (v_char^4);        % km^6/s^4-like
    scale_auto = max(scale_auto, opts.scale_floor);

    if isstring(opts.scale_mode) || ischar(opts.scale_mode)
        mode = string(opts.scale_mode);
    else
        mode = "auto";
    end

    switch mode
        case "auto"
            term_scale = scale_auto;
        case "mu"
            % Sometimes scaling by mu^2 is fine if your equation is mu-dominated
            term_scale = max(mu_1^2, opts.scale_floor);
        case "fixed"
            term_scale = max(opts.scale_value, opts.scale_floor);
        otherwise
            term_scale = scale_auto;
    end

    % ---------- Residual function ----------
    function y = term_of_rp(rp)
        % Hard invalid checks
        if ~isfinite(rp) || rp <= 0
            y = opts.penalty; return
        end

        % beta in (0,1) for trig identity form
        beta = 1 / (1 + (rp * v_inf^2)/mu_2);
        if ~(isfinite(beta) && beta > 0 && beta < 1)
            y = opts.penalty; return
        end

        % Signed rotation parameters
        c = (1 - 2*beta^2);
        s = sign_turn * (2*beta*sqrt(1 - beta^2));

        % Post-flyby heliocentric RT components
        v_r = Vr_b2 + v_inf_in_r*c - v_inf_in_t*s;
        v_t = Vt_b2 + v_inf_in_r*s + v_inf_in_t*c;

        A = v_t*sin(theta) - v_r*cos(theta);

        % Raw residual
        raw = mu_1*sin(theta)^2 - r1^2*A^2 - r1*A*v_r*cos(theta_r);

        % Scale for conditioning (DOES NOT change the root location)
        y = raw / term_scale;

        if ~isfinite(y)
            y = opts.penalty;
        end
    end

    % ---------- Grid scan for bracketing ----------
    r_grid = logspace(log10(rp_min), log10(rp_max), opts.n_grid);
    f_grid = arrayfun(@term_of_rp, r_grid);

    % Mark "valid" (not penalty-ish)
    good = isfinite(f_grid) & abs(f_grid) < 0.1*opts.penalty;
    r_good = r_grid(good);
    f_good = f_grid(good);

    info = struct();
    info.r_grid = r_grid;
    info.f_grid = f_grid;
    info.good_mask = good;
    info.term_scale = term_scale;
    info.scale_mode = mode;
    info.v_char = v_char;

    rp_sol = NaN;
    info.status = "no_root";
    info.bracket = [NaN NaN];
    info.best_rp = NaN;
    info.best_abs_term = NaN;

    if ~isempty(f_good)
        [bestAbs, kbest] = min(abs(f_good));
        info.best_rp = r_good(kbest);
        info.best_abs_term = bestAbs;
    end

    % Find first sign change in good region
    idx = [];
    if numel(f_good) >= 2
        idx = find(sign(f_good(1:end-1)) ~= sign(f_good(2:end)), 1, 'first');
    end

    if isempty(idx)
        info.status = "no_sign_change";
        if opts.verbose
            fprintf("[solve_rp_numeric] No sign change in [%.3g, %.3g] km.\n", rp_min, rp_max);
            fprintf("  scale_mode=%s, term_scale=%.3e, v_char=%.3g\n", mode, term_scale, v_char);
            fprintf("  Best |scaled term|=%.3e at rp=%.3g km\n", info.best_abs_term, info.best_rp);
        end
        return
    end

    rp_lo = r_good(idx);
    rp_hi = r_good(idx+1);
    info.bracket = [rp_lo rp_hi];

    % ---------- Solve with fzero ----------
    try
        rp_sol = fzero(@term_of_rp, [rp_lo rp_hi]);
        info.status = "ok";
        info.term_at_solution_scaled = term_of_rp(rp_sol);
        info.term_at_solution_raw = info.term_at_solution_scaled * term_scale;

        if opts.verbose
            fprintf("[solve_rp_numeric] rp=%.6g km\n", rp_sol);
            fprintf("  scaled term=%.3e (raw=%.3e)\n", info.term_at_solution_scaled, info.term_at_solution_raw);
            fprintf("  bracket=[%.3g, %.3g] km\n", rp_lo, rp_hi);
            fprintf("  scale_mode=%s, term_scale=%.3e, v_char=%.3g\n", mode, term_scale, v_char);
        end
    catch ME
        info.status = "fzero_failed";
        info.error = ME.message;
        rp_sol = NaN;

        if opts.verbose
            fprintf("[solve_rp_numeric] fzero failed: %s\n", ME.message);
        end
    end
end