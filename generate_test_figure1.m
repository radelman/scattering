function generate_test_figure1(k, description)
	a = 1.0;
	theta0 = pi;
	xi1 = 1.5;
	w = 5.0;
	nels = 200;
	[ ...
	x, y, z ...
	] = create_slice([-w; 0.0; w], [0.0; 0.0; -2.0 * w], [2.0 * w; 0.0; 0.0], nels, nels);
	x = reshape(x, 1, nels * nels);
	y = reshape(y, 1, nels * nels);
	z = reshape(z, 1, nels * nels);
	v_in = plane_wave_in(k, sin(theta0), 0.0, cos(theta0), x, y, z);
	[ ...
	v_scat, max_abs_change ...
	] = pro_plane_wave_scat_hard(k, a, theta0, 'saved', xi1, x, y, z);
	v = v_in + v_scat;
	bx = dot(pro_to_cart(a, [0.0; xi1; 0.0]), [1.0; 0.0; 0.0]);
	by = dot(pro_to_cart(a, [1.0; xi1; 0.0]), [0.0; 0.0; 1.0]);
	x = reshape(x, nels, nels);
	y = reshape(x, nels, nels);
	z = reshape(x, nels, nels);
	v_in = reshape(v_in, nels, nels);
	v_scat = reshape(v_scat, nels, nels);
	v = reshape(v, nels, nels);
	cool_plot([-w, w], [-w, w], 'x (m)', 'z (m)', 'Field Strength (Non-Dimensional)', v, 0.0, [-2.0, 2.0], jet(256), 0.0, 0.0, bx, by);
	title(description);
end
