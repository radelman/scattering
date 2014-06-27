function grad_cart = grad_pro_to_cart(a, pro, grad_pro)
	grad_cart = zeros(3, size(pro, 2));
	for i = 1 : size(pro, 2)
		D = pro_calculate_D(a, pro(1 : 3, i));
		if (abs(pro(1, i)) < 1.0)
			grad_cart(1 : 3, i) = D \ grad_pro(1 : 3, i);
		else
			grad_cart(1 : 3, i) = [nan; nan; grad_pro(2, i) / D(2, 3)];
		end
	end
end
