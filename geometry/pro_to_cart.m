function cart = pro_to_cart(a, pro)
	eta = pro(1, :);
	xi = pro(2, :);
	phi = pro(3, :);
	
	x = a * sqrt((1.0 - (eta .^ 2)) .* ((xi .^ 2) - 1.0)) .* cos(phi);
	y = a * sqrt((1.0 - (eta .^ 2)) .* ((xi .^ 2) - 1.0)) .* sin(phi);
	z = a * eta .* xi;
	
	cart = [x; y; z];
end
