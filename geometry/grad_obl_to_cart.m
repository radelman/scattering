function grad_cart = grad_obl_to_cart(a, obl, grad_obl)
	grad_cart = zeros(3, size(obl, 2));
	for i = 1 : size(obl, 2)
		D = obl_calculate_D(a, obl(1 : 3, i));
		if (abs(obl(1, i)) < 1.0)
			grad_cart(1 : 3, i) = D \ grad_obl(1 : 3, i);
		else
			grad_cart(1 : 3, i) = [nan; nan; grad_obl(2, i) / D(2, 3)];
		end
	end
end
