function generate_figure3()
	figure();
	hold('on');
	XI = 0.0 : 1.0 : 4.0;
	for xi = XI
		eta = linspace(-1.0, 1.0, 1001);
		obl = [eta; xi * ones(1, length(eta)); 0.0 * ones(1, length(eta))];
		cart = obl_to_cart(1.0, obl);
		plot(cart(1, :), cart(3, :), 'linewidth', 1.0, 'color', [0.0, 0.0, 0.0]);
		obl = [eta; xi * ones(1, length(eta)); pi * ones(1, length(eta))];
		cart = obl_to_cart(1.0, obl);
		plot(cart(1, :), cart(3, :), 'linewidth', 1.0, 'color', [0.0, 0.0, 0.0]);
	end
	ETA = -1.0 : 0.25 : 1.0;
	for eta = ETA
		xi = linspace(0.0, 10.0, 1001);
		obl = [eta * ones(1, length(xi)); xi; 0.0 * ones(1, length(xi))];
		cart = obl_to_cart(1.0, obl);
		plot(cart(1, :), cart(3, :), 'linewidth', 1.0, 'color', [0.0, 0.0, 0.0]);
		obl = [eta * ones(1, length(xi)); xi; pi * ones(1, length(xi))];
		cart = obl_to_cart(1.0, obl);
		plot(cart(1, :), cart(3, :), 'linewidth', 1.0, 'color', [0.0, 0.0, 0.0]);
	end
	axis('equal');
	xlim([-4.5, 4.5]);
	ylim([-4.5, 4.5]);
	axis('off');
	set(gcf, 'position', [100, 100, 800, 600]);
	set(gcf, 'color', [1.0, 1.0, 1.0]);
end
