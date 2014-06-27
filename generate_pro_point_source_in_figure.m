function generate_pro_point_source_in_figure(k, a, x0, z0, w, nels)
	[ ...
	x, y, z ...
	] = create_slice([-w; 0.0; w], [0.0; 0.0; -2.0 * w], [2.0 * w; 0.0; 0.0], nels, nels);
	x = reshape(x, 1, nels * nels);
	y = reshape(y, 1, nels * nels);
	z = reshape(z, 1, nels * nels);
	[ ...
	v_in, grad_in_cart ...
	] = point_source_in(k, x0, 0.0, z0, x, y, z);
	cart = [x0; 0.0; z0];
	pro = cart_to_pro(a, cart);
	eta0 = pro(1);
	xi0 = pro(2);
	[ ...
	v_scat, grad_scat_cart, max_abs_change_scat ...
	] = pro_point_source_in(k, a, eta0, xi0, 'saved', x, y, z);
	v_in = reshape(v_in, nels, nels);
	v_scat = reshape(v_scat, nels, nels);
	
	cool_plot([-w, w], [-w, w], 'x (m)', 'z (m)', '', log10(abs(v_scat - v_in)), 0.0, [-6.0, 0.0], cubehelix(250), 0.0, 0.0);
end
