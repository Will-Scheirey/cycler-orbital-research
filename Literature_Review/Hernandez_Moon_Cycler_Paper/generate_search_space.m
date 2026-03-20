function out = generate_search_space(T_syn, syn_per_max, T_min)
out = zeros(syn_per_max, 1);

for syn_per = 1:syn_per_max
    out(syn_per) = floor(T_syn * syn_per / T_min);
end
end