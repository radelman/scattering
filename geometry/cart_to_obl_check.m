a = 10.0 * rand();
for i = 1 : 10000
	cart = -10.0 + 20.0 * rand(3, 1);
	obl = cart_to_obl(a, cart);
	if (~isreal(obl) || obl(1) < -1.0 || obl(1) > 1.0 || obl(2) < 0.0)
		cart
		obl
		error('complex or out of range');
	end
	cart_back = obl_to_cart(a, obl);
	norm(cart_back - cart)
end
