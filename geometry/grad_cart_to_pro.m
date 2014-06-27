function grad_pro = grad_cart_to_pro(a, cart, grad_cart)
	pro = cart_to_pro(a, cart);
	grad_pro = zeros(3, size(cart, 2));
	for i = 1 : size(cart, 2)
		D = pro_calculate_D(a, pro(1 : 3, i));
		if (abs(pro(1, i)) < 1.0)
			grad_pro(1 : 3, i) = D * grad_cart(1 : 3, i);
		else
			grad_pro(1 : 3, i) = [nan; D(2, 3) * grad_cart(3, i); nan];
		end
	end
end
