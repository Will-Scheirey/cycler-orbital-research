function [S, angle_turn] = calc_synodic(mu, r_b1, r_b2)
    T1 = 2*pi * sqrt(r_b1 ^ 3 / mu);
    T2 = 2*pi * sqrt(r_b2 ^ 3 / mu);
    
    T_ratio = T2 / T1;
    
    S = 1 / abs(1/T1 - 1/T2);
    
    angle_turn = mod(T_ratio * 2*pi, 2*pi);
    if angle_turn > pi
        angle_turn = angle_turn - 2*pi;
    end
    
    angle_turn = -angle_turn;
end