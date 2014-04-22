function pro = cart_to_pro(a, cart)
	x = cart(1, :);
	y = cart(2, :);
	z = cart(3, :);
	
	d = a ^ 2;
	e = -(x .^ 2) - (y .^ 2) - (z .^ 2) - (a ^ 2);
	g = z .^ 2;
	xi = sqrt((-e + sqrt((e .^ 2) - 4.0 * d * g)) / (2.0 * d));
	eta = z ./ (a * xi);
	phi = atan2(y, x);
	
	pro = [eta; xi; phi];
end
