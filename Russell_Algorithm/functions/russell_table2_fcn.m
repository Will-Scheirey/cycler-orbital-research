function [fj, dt_years] = russell_table2_fcn(hj)

if mod(hj,2) == 0
    % even case
    fj = hj/2 + 1;
    dt_years = ones(fj-1,1);

else
    % odd case
    fj = 2*floor(hj/4 + 1);

    dt_years = ones(fj-1,1);

    % middle interval index
    k_mid = fj/2;

    dt_years(k_mid) = 0.5 * mod(hj,4);
end
end