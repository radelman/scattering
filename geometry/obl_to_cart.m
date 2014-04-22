function cart = obl_to_cart(a, obl)
	eta = obl(1, :);
	xi = obl(2, :);
	phi = obl(3, :);
	
	x = a * sqrt((1.0 - (eta .^ 2)) .* ((xi .^ 2) + 1.0)) .* cos(phi);
	y = a * sqrt((1.0 - (eta .^ 2)) .* ((xi .^ 2) + 1.0)) .* sin(phi);
	z = a * eta .* xi;
	
	cart = [x; y; z];
end
