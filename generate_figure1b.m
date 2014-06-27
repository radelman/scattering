%
% Generate a figure showing the oblate spheroidal coordinate system.
%
function generate_figure1b()
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
	obl = cart_to_obl(1.0, cart);
	eta = obl(1, :);
	xi = obl(2, :);
	phi = obl(3, :);
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
	fv = isosurface(x, y, z, xi, 0.5);
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
	obl = cart_to_obl(1.0, cart);
	eta = obl(1, :);
	xi = obl(2, :);
	phi = obl(3, :);
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
	export_fig(gcf, 'images/figure1b.png');
	close();
end
