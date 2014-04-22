function obl = cart_to_obl(a, cart)
	x = cart(1, :);
	y = cart(2, :);
	z = cart(3, :);
	
	d = a ^ 2;
	e = -(x .^ 2) - (y .^ 2) - (z .^ 2) + (a ^ 2);
	g = -(z .^ 2);
	xi = sqrt((-e + sqrt((e .^ 2) - 4.0 * d * g)) / (2.0 * d));
	eta = sqrt(max(zeros(1, length(x)), 1.0 - ((x .^ 2) + (y .^ 2)) ./ ((a ^ 2) * ((xi .^ 2) + 1.0))));
	eta(z < 0.0) = -eta(z < 0.0);
	phi = atan2(y, x);
	
	obl = [eta; xi; phi];
end
