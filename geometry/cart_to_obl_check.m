for i = 1 : 10000
	cart = -1.0 + 2.0 * rand(3, 1);
	obl = cart_to_obl(2.0, cart);
	cart_back = obl_to_cart(2.0, obl);
	norm(cart_back - cart)
end
