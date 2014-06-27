function [x, y, z] = create_slice(corner, a, b, a_nels, b_nels)
	x = zeros(a_nels, b_nels);
	y = zeros(a_nels, b_nels);
	z = zeros(a_nels, b_nels);
	for i = 1 : a_nels
		for j = 1 : b_nels
			p = corner + ((i - 1) / (a_nels - 1)) * a + ((j - 1) / (b_nels - 1)) * b;
			x(i, j) = p(1);
			y(i, j) = p(2);
			z(i, j) = p(3);
		end
	end
end
