function families = generate_families(num_flybys, num_bodies)

    num_families = num_bodies^num_flybys;
    
    families = cell(num_families, num_flybys);

    for n = 1:num_flybys
        for b=1:num_bodies
            families{n, }
        end
    end
end