function generate_figure1()
	figure();
	hold('on');
	XI = 1.0 : 1.0 : 4.0;
	for xi = XI
		eta = linspace(-1.0, 1.0, 1001);
		pro = [eta; xi * ones(1, length(eta)); 0.0 * ones(1, length(eta))];
		cart = pro_to_cart(1.0, pro);
		plot(cart(1, :), cart(3, :), 'linewidth', 1.0, 'color', [0.0, 0.0, 0.0]);
		pro = [eta; xi * ones(1, length(eta)); pi * ones(1, length(eta))];
		cart = pro_to_cart(1.0, pro);
		plot(cart(1, :), cart(3, :), 'linewidth', 1.0, 'color', [0.0, 0.0, 0.0]);
	end
	ETA = -1.0 : 0.25 : 1.0;
	for eta = ETA
		xi = linspace(1.0, 10.0, 1001);
		pro = [eta * ones(1, length(xi)); xi; 0.0 * ones(1, length(xi))];
		cart = pro_to_cart(1.0, pro);
		plot(cart(1, :), cart(3, :), 'linewidth', 1.0, 'color', [0.0, 0.0, 0.0]);
		pro = [eta * ones(1, length(xi)); xi; pi * ones(1, length(xi))];
		cart = pro_to_cart(1.0, pro);
		plot(cart(1, :), cart(3, :), 'linewidth', 1.0, 'color', [0.0, 0.0, 0.0]);
	end
	axis('equal');
	xlim([-4.5, 4.5]);
	ylim([-4.5, 4.5]);
	axis('off');
% 	text(4.55, -0.05, 'x, y', 'fontsize', 20);
% 	text(-0.1, 4.65, 'z', 'fontsize', 20);
% 	text(-0.5, -4.8, '\eta = -1', 'fontsize', 20);
% 	text(3.2, -4.8, '\eta = -0.75', 'fontsize', 20);
% 	text(4.55, -2.65, '\eta = -0.5', 'fontsize', 20);
% 	text(4.55, -1.2, '\eta = -0.25', 'fontsize', 20);
% 	text(4.55, 1.2, '\eta = 0.25', 'fontsize', 20);
% 	text(4.55, 2.65, '\eta = 0.5', 'fontsize', 20);
% 	text(3.2, 4.65, '\eta = 0.75', 'fontsize', 20);
% 	text(0.15, 2.23, '\xi = 2', 'fontsize', 20);
% 	text(0.25, 3.23, '\xi = 3', 'fontsize', 20);
% 	text(0.35, 4.23, '\xi = 4', 'fontsize', 20);
	set(gcf, 'position', [100, 100, 800, 600]);
	set(gcf, 'color', [1.0, 1.0, 1.0]);
end
