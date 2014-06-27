function pro = cart_to_pro(a, cart)
	x = cart(1, :);
	y = cart(2, :);
	z = cart(3, :);
	b = -(a ^ 2);
	c = (a ^ 2) + (x .^ 2) + (y .^ 2) + (z .^ 2);
	d = -(z .^ 2);
	xi = sqrt((-c - sqrt((c .^ 2) - 4.0 * b * d)) / (2.0 * b));
	eta = z ./ (a * xi);
	phi = atan2(y, x);
	pro = [eta; xi; phi];
end
