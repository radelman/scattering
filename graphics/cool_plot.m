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

function cool_plot(xr, yr, xlb, ylb, zlb, phis, angles, cr, cm, ax, ay)
	figure();
	line(xr, yr);
	hold('on');
	axis('equal');
	xlim(xr);
	ylim(yr);
	draw_border();
	set(gca, 'fontsize', 30);
	xlabel(xlb);
	ylabel(ylb);
	set(gca, 'fontsize', 20);
	set(gcf, 'position', [100, 100, 800, 600]);
	set(gcf, 'color', [1.0, 1.0, 1.0]);
	imh = 0;
	for theta = angles
		if (imh == 0)
			imh = imagesc(xr, yr(2 : -1 : 1), real(phis * exp(-1i * theta)));
			if (ax > 0.0 && ay > 0.0)
				rectangle('position', [-ax, -ay, 2.0 * ax, 2.0 * ay], 'curvature', [1.0, 1.0], 'edgecolor', [0.001, 0.001, 0.001], 'facecolor', [0.001, 0.001, 0.001]);
			end
			colormap(cm);
			cb = colorbar();
			caxis(cr);
			set(cb, 'ylim', cr);
			set(cb, 'ticklength', [1.0e-9, 1.0e-9]);
			set(cb, 'layer', 'top');
			set(cb, 'fontsize', 30);
			set(get(cb, 'ylabel'), 'fontsize', 30, 'string', zlb);
			set(cb, 'fontsize', 20);
		else
			set(imh, 'cdata', real(phis * exp(-1i * theta)));
		end
		drawnow();
		pause(0.1);
	end
end
