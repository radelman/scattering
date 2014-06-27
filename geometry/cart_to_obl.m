function obl = cart_to_obl(a, cart)
	x = cart(1, :);
	y = cart(2, :);
	z = cart(3, :);
	b = -(a ^ 2);
	c = (x .^ 2) + (y .^ 2) + (z .^ 2) - (a ^ 2);
	d = z .^ 2;
	xi = sqrt((-c - sqrt((c .^ 2) - 4.0 * b * d)) / (2.0 * b));
	eta = sqrt(max(zeros(1, length(x)), 1.0 - ((x .^ 2) + (y .^ 2)) ./ ((a ^ 2) * ((xi .^ 2) + 1.0))));
	eta(z < 0.0) = -eta(z < 0.0);
	phi = atan2(y, x);
	obl = [eta; xi; phi];
end
