function out = planet_state_from_datetime(params, planetName, t)
% Convenience wrapper: datetime -> JD -> state
    JD = datetime_to_jd(t);
    out = planet_state_from_jd(params, planetName, JD);
    out.JD = JD;
    out.datetime_utc = datetime(t, 'TimeZone', 'UTC');
end