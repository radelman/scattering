%
% Generate a figure showing the prolate spheroidal coordinate system.
%
function generate_figure1a()
	w = 2.0;
	nels = 100;
	x = linspace(-w, w, nels);
	y = linspace(-w, w, nels);
	z = linspace(-w, w, nels);
	[ ...
	x, y, z ...
	] = meshgrid(x, y, z);
	x = reshape(x, 1, nels * nels * nels);
	y = reshape(y, 1, nels * nels * nels);
	z = reshape(z, 1, nels * nels * nels);
	cart = [x; y; z];
	pro = cart_to_pro(1.0, cart);
	eta = pro(1, :);
	xi = pro(2, :);
	phi = pro(3, :);
	x = reshape(x, nels, nels, nels);
	y = reshape(y, nels, nels, nels);
	z = reshape(z, nels, nels, nels);
	eta = reshape(eta, nels, nels, nels);
	xi = reshape(xi, nels, nels, nels);
	phi = reshape(phi, nels, nels, nels);
	fv = isosurface(x, y, z, eta, -0.5);
	trisurf(fv.faces, fv.vertices(:, 1), fv.vertices(:, 2), fv.vertices(:, 3), 'edgecolor', 'none', 'facecolor', [1.0, 0.0, 0.0]);
	hold('on');
	fv = isosurface(x, y, z, eta, 0.5);
	trisurf(fv.faces, fv.vertices(:, 1), fv.vertices(:, 2), fv.vertices(:, 3), 'edgecolor', 'none', 'facecolor', [1.0, 0.0, 0.0]);
	fv = isosurface(x, y, z, xi, 1.5);
	trisurf(fv.faces, fv.vertices(:, 1), fv.vertices(:, 2), fv.vertices(:, 3), 'edgecolor', 'none', 'facecolor', [0.0, 0.7, 0.0]);
	x = linspace(0.0, w, nels);
	y = linspace(-w, w, nels);
	z = linspace(-w, w, nels);
	[ ...
	x, y, z ...
	] = meshgrid(x, y, z);
	x = reshape(x, 1, nels * nels * nels);
	y = reshape(y, 1, nels * nels * nels);
	z = reshape(z, 1, nels * nels * nels);
	cart = [x; y; z];
	pro = cart_to_pro(1.0, cart);
	eta = pro(1, :);
	xi = pro(2, :);
	phi = pro(3, :);
	x = reshape(x, nels, nels, nels);
	y = reshape(y, nels, nels, nels);
	z = reshape(z, nels, nels, nels);
	eta = reshape(eta, nels, nels, nels);
	xi = reshape(xi, nels, nels, nels);
	phi = reshape(phi, nels, nels, nels);
	fv = isosurface(x, y, z, phi, 0.0);
	trisurf(fv.faces, fv.vertices(:, 1), fv.vertices(:, 2), fv.vertices(:, 3), 'edgecolor', 'none', 'facecolor', [0.0, 0.0, 1.0]);
	alpha(0.5);
	axis('equal');
	xlim([-w, w]);
	ylim([-w, w]);
	zlim([-w, w]);
	axis('off');
	view(50.0, 10.0);
	set(gcf, 'position', [100, 100, 800, 600]);
	set(gcf, 'color', [1.0, 1.0, 1.0]);
	export_fig(gcf, 'images/figure1a.png');
	close();
end
