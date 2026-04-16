function plot_intersection(r, offset_x, text_str)
plot(r(1), r(2), 'k.', 'MarkerSize', 25, 'LineWidth', 2, 'HandleVisibility', 'off')
text(r(1) + offset_x, r(2), text_str)
end