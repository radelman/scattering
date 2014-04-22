function generate_figure6(k, description)
	a = 1.0;
	x0 = 3.0;
	z0 = 3.0;
	xi1 = 0.0;
	w = 5.0;
	nels = 200;
	[ ...
	x, y, z ...
	] = create_slice([-w; 0.0; w], [0.0; 0.0; -2.0 * w], [2.0 * w; 0.0; 0.0], nels, nels);
	x = reshape(x, 1, nels * nels);
	y = reshape(y, 1, nels * nels);
	z = reshape(z, 1, nels * nels);
	v_in = point_source_in(k, x0, 0.0, z0, x, y, z);
	cart = [x0; 0.0; z0];
	obl = cart_to_obl(a, cart);
	eta0 = obl(1);
	xi0 = obl(2);
	[ ...
	v_scat, max_abs_change ...
	] = obl_point_source_scat_hard(k, a, eta0, xi0, 'saved', xi1, x, y, z);
	v = v_in + v_scat;
	bx = dot(obl_to_cart(a, [0.0; xi1; 0.0]), [1.0; 0.0; 0.0]);
	by = dot(obl_to_cart(a, [1.0; xi1; 0.0]), [0.0; 0.0; 1.0]);
	s = abs(point_source_in(k, x0, 0.0, z0, 0.0, 0.0, 0.0));
	v_in = v_in / s;
	v_scat = v_scat / s;
	v = v / s;
	x = reshape(x, nels, nels);
	y = reshape(x, nels, nels);
	z = reshape(x, nels, nels);
	v_in = reshape(v_in, nels, nels);
	v_scat = reshape(v_scat, nels, nels);
	v = reshape(v, nels, nels);
	cool_plot([-w, w], [-w, w], 'x (m)', 'z (m)', 'Field Strength (Non-Dimensional)', v, 0.0, [-2.0, 2.0], jet(256), 0.0, 0.0, bx, by);
	title(description);
end
