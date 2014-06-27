function grad_obl = grad_cart_to_obl(a, cart, grad_cart)
	obl = cart_to_obl(a, cart);
	grad_obl = zeros(3, size(cart, 2));
	for i = 1 : size(cart, 2)
		D = obl_calculate_D(a, obl(1 : 3, i));
		if (abs(obl(1, i)) < 1.0)
			grad_obl(1 : 3, i) = D * grad_cart(1 : 3, i);
		else
			grad_obl(1 : 3, i) = [nan; D(2, 3) * grad_cart(3, i); nan];
		end
	end
end
