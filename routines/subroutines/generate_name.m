function name = generate_name(c, m, n)
	name = sprintf('%08d_%03d_%03d', round(1000.0 * c), m, n);
end
