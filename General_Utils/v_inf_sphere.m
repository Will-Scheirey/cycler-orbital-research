clear; clc; close all

AU = 1.496e+8;
mu = 132712000000;

ve = sqrt(mu / AU);

v_inf = 10;

ve_vec = [0, -ve, 0];
v_inf_vec = v_inf * [1, 1, 1] / sqrt(3);

draw_v_inf_sphere(ve_vec, v_inf_vec)
