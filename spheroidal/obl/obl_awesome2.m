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

function obl_awesome2(path, c, m, N, n_dr, n_dr_neg, n_c2k, n_B2r, description, dr_caxis, c2k_caxis, B2r_caxis)
	log_abs_dr = [];
	log_abs_c2k = [];
	log_abs_B2r = [];
	for n = m : m + N - 1
		everything = obl_open_everything(path, c, m, n);
		log_abs_dr = [log_abs_dr; everything.coefficients.log_abs_dr_neg(n_dr_neg : -1 : 1), everything.coefficients.log_abs_dr(1 : n_dr)];
		log_abs_c2k = [log_abs_c2k; everything.coefficients.log_abs_c2k(1 : n_c2k)];
		log_abs_B2r = [log_abs_B2r; everything.coefficients.log_abs_B2r(1 : n_B2r)];
	end
	
	figure();
	line([0, N - 1], [-n_dr_neg, n_dr - 1]);
	hold('on');
	imagesc([0, N - 1], [-n_dr_neg, n_dr - 1], log_abs_dr.' / log(10.0));
	xlim([0, N - 1]);
	set(gca, 'xtick', round(linspace(0, N - 1, 4)));
	ylim([-n_dr_neg, n_dr - 1]);
	draw_border();
	set(gca, 'fontsize', 20);
	xlabel('n - m');
	ylabel('r');
	title(sprintf('log_{10}(|d_r^{mn}(-ic)|)%s', description));
	colormap(cubehelix(250));
	cb = colorbar();
	if (~isempty(dr_caxis))
		caxis([min(dr_caxis), max(dr_caxis)]);
		set(cb, 'ylim', [min(dr_caxis), max(dr_caxis)]);
		set(cb, 'ytick', dr_caxis);
	end
	set(cb, 'ticklength', [1.0e-9, 1.0e-9]);
	set(cb, 'layer', 'top');
	set(cb, 'fontsize', 20);
	set(gcf, 'position', [100, 100, 800, 600]);
	set(gcf, 'color', [1.0, 1.0, 1.0]);
	
	figure();
	line([0, N - 1], [0, n_c2k - 1]);
	hold('on');
	imagesc([0, N - 1], [0, n_c2k - 1], log_abs_c2k.' / log(10.0));
	xlim([0, N - 1]);
	set(gca, 'xtick', round(linspace(0, N - 1, 4)));
	ylim([0, n_c2k - 1]);
	draw_border();
	set(gca, 'fontsize', 20);
	xlabel('n - m');
	ylabel('k');
	title(sprintf('log_{10}(|c_{2k}^{mn}(-ic)|)%s', description));
	colormap(cubehelix(250));
	cb = colorbar();
	if (~isempty(c2k_caxis))
		caxis([min(c2k_caxis), max(c2k_caxis)]);
		set(cb, 'ylim', [min(c2k_caxis), max(c2k_caxis)]);
		set(cb, 'ytick', c2k_caxis);
	end
	set(cb, 'ticklength', [1.0e-9, 1.0e-9]);
	set(cb, 'layer', 'top');
	set(cb, 'fontsize', 20);
	set(gcf, 'position', [100, 100, 800, 600]);
	set(gcf, 'color', [1.0, 1.0, 1.0]);
	
	figure();
	line([0, N - 1], [0, n_B2r - 1]);
	hold('on');
	imagesc([0, N - 1], [0, n_B2r - 1], log_abs_B2r.' / log(10.0));
	xlim([0, N - 1]);
	set(gca, 'xtick', round(linspace(0, N - 1, 4)));
	ylim([0, n_B2r - 1]);
	draw_border();
	set(gca, 'fontsize', 20);
	xlabel('n - m');
	ylabel('k');
	title(sprintf('log_{10}(|B_{2r}^{mn}(-ic)|)%s', description));
	colormap(cubehelix(250));
	cb = colorbar();
	if (~isempty(B2r_caxis))
		caxis([min(B2r_caxis), max(B2r_caxis)]);
		set(cb, 'ylim', [min(B2r_caxis), max(B2r_caxis)]);
		set(cb, 'ytick', B2r_caxis);
	end
	set(cb, 'ticklength', [1.0e-9, 1.0e-9]);
	set(cb, 'layer', 'top');
	set(cb, 'fontsize', 20);
	set(gcf, 'position', [100, 100, 800, 600]);
	set(gcf, 'color', [1.0, 1.0, 1.0]);
end
