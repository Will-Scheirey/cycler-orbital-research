function [dv_min_kms, best] = min_dv_lambert_universal(r1_km, v1p_kms, r2_km, v2p_kms, tof_s, mu_km, z_lim, num_z)
% Returns minimum dv (km/s) over all Lambert solutions returned by lambert()

    best = struct('V1', [], 'V2', [], 'dv1', NaN, 'dv2', NaN, 'dv', Inf);

    [V1_list, V2_list] = lambert(r1_km, r2_km, tof_s, mu_km, z_lim, num_z);

    if isempty(V1_list)
        dv_min_kms = NaN;
        return;
    end

    dv_min_kms = Inf;

    for k = 1:numel(V1_list)
        V1 = V1_list{k};
        V2 = V2_list{k};

        dv1 = norm(V1 - v1p_kms);
        dv2 = norm(V2 - v2p_kms);
        dv  = dv1 + dv2;

        if dv < dv_min_kms
            dv_min_kms = dv;
            best.V1 = V1; best.V2 = V2;
            best.dv1 = dv1; best.dv2 = dv2; best.dv = dv;
        end
    end

    if ~isfinite(dv_min_kms)
        dv_min_kms = NaN;
    end
end