%
% Copyright (c) 2014, Ross Adelman, Nail A. Gumerov, and Ramani Duraiswami
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%

function generate_pro_plane_wave_scat_figure(k, a, theta0, xi1, alpha, w, nels)
	eta = linspace(-1.0, 1.0, 1000);
	xi = xi1 * ones(1, length(eta));
	phi = zeros(1, length(eta));
	pro = [eta; xi; phi];
	cart = pro_to_cart(a, pro);
	x = cart(1, :);
	y = cart(2, :);
	z = cart(3, :);
	[ ...
	v_in, grad_in_cart ...
	] = plane_wave_in(k, sin(theta0), 0.0, cos(theta0), x, y, z);
	grad_in_pro = grad_cart_to_pro(a, [x; y; z], grad_in_cart);
	gradxi_in = grad_in_pro(2, :);
	if (ischar(alpha))
		if (strcmp(alpha, 'soft'))
			[ ...
			v_scat, grad_scat_cart, max_abs_change_scat ...
			] = pro_plane_wave_scat_soft(k, a, theta0, 'saved', xi1, x, y, z);
		else
			[ ...
			v_scat, grad_scat_cart, max_abs_change_scat ...
			] = pro_plane_wave_scat_hard(k, a, theta0, 'saved', xi1, x, y, z);
		end
	else
		[ ...
		v_scat, grad_scat_cart, max_abs_change_scat ...
		] = pro_plane_wave_scat_robin(k, a, theta0, 'saved', xi1, alpha, x, y, z);
	end
	grad_scat_pro = grad_cart_to_pro(a, [x; y; z], grad_scat_cart);
	gradxi_scat = grad_scat_pro(2, :);
	v = v_in + v_scat;
	gradxi = gradxi_in + gradxi_scat;
	
	figure();
	semilogy(eta, abs(v), 'k');
	hold('on');
	semilogy(eta, abs(gradxi), 'r');
	if (~ischar(alpha))
		semilogy(eta, abs(v + alpha * gradxi), 'b');
	end
	
	[ ...
	x, y, z ...
	] = create_slice([-w; 0.0; w], [0.0; 0.0; -2.0 * w], [2.0 * w; 0.0; 0.0], nels, nels);
	x = reshape(x, 1, nels * nels);
	y = reshape(y, 1, nels * nels);
	z = reshape(z, 1, nels * nels);
	[ ...
	v_in, grad_in_cart ...
	] = plane_wave_in(k, sin(theta0), 0.0, cos(theta0), x, y, z);
	grad_in_pro = grad_cart_to_pro(a, [x; y; z], grad_in_cart);
	gradxi_in = grad_in_pro(2, :);
	if (ischar(alpha))
		if (strcmp(alpha, 'soft'))
			[ ...
			v_scat, grad_scat_cart, max_abs_change_scat ...
			] = pro_plane_wave_scat_soft(k, a, theta0, 'saved', xi1, x, y, z);
		else
			[ ...
			v_scat, grad_scat_cart, max_abs_change_scat ...
			] = pro_plane_wave_scat_hard(k, a, theta0, 'saved', xi1, x, y, z);
		end
	else
		[ ...
		v_scat, grad_scat_cart, max_abs_change_scat ...
		] = pro_plane_wave_scat_robin(k, a, theta0, 'saved', xi1, alpha, x, y, z);
	end
	grad_scat_pro = grad_cart_to_pro(a, [x; y; z], grad_scat_cart);
	gradxi_scat = grad_scat_pro(2, :);
	v = v_in + v_scat;
	gradxi = gradxi_in + gradxi_scat;
	v = reshape(v, nels, nels);
	gradxi = reshape(gradxi, nels, nels);
	ax = dot(pro_to_cart(a, [0.0; xi1; 0.0]), [1.0; 0.0; 0.0]);
	ay = dot(pro_to_cart(a, [1.0; xi1; 0.0]), [0.0; 0.0; 1.0]);
	
	cool_plot([-w, w], [-w, w], 'x (m)', 'z (m)', '', v, 0.0, [-2.0, 2.0], cubehelix(250), ax, ay);
	cool_plot([-w, w], [-w, w], 'x (m)', 'z (m)', '|V|', abs(v), 0.0, [0.0, 2.0], cubehelix(250), ax, ay);
	cool_plot([-w, w], [-w, w], 'x (m)', 'z (m)', '', log10(abs(v)), 0.0, [-6.0, 0.0], cubehelix(250), ax, ay);
	cool_plot([-w, w], [-w, w], 'x (m)', 'z (m)', '', log10(abs(gradxi)), 0.0, [-6.0, 0.0], cubehelix(250), ax, ay);
	if (~ischar(alpha))
		cool_plot([-w, w], [-w, w], 'x (m)', 'z (m)', '', log10(abs(v + alpha * gradxi)), 0.0, [-6.0, 0.0], cubehelix(250), ax, ay);
	end
end
