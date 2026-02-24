function JD = datetime_to_jd(t)
%DATETIME_TO_JD Convert MATLAB datetime to Julian Date (UTC-based).
%
% Notes:
% - This treats the input as UTC if no TimeZone is set.
% - It does NOT convert UTC->TT/TDB. The difference is ~minutes/seconds scale.
%   For visualization / rough positions this is typically fine.

    if isempty(t.TimeZone)
        t.TimeZone = 'UTC';
    else
        t = datetime(t, 'TimeZone', 'UTC');
    end

    % Julian Date for Unix epoch 1970-01-01 00:00:00 UTC is 2440587.5
    JD = 2440587.5 + posixtime(t) / 86400;
end