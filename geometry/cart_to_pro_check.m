a = 10.0 * rand();
for i = 1 : 10000
	cart = -10.0 + 20.0 * rand(3, 1);
	pro = cart_to_pro(a, cart);
	if (~isreal(pro) || pro(1) < -1.0 || pro(1) > 1.0 || pro(2) < 1.0)
		cart
		pro
		error('complex or out of range');
	end
	cart_back = pro_to_cart(a, pro);
	norm(cart_back - cart)
end
