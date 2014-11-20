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

function obl_awesome3(path, c, m, N, max_xi, description, axis_width, axis_height)
	figure();
	subplot(2, 3, 1);
	hold('on');
	nice_ylim1 = [inf, -inf];
	subplot(2, 3, 4);
	hold('on');
	nice_ylim4 = [inf, -inf];
	subplot(2, 3, 2);
	hold('on');
	nice_ylim2 = [inf, -inf];
	subplot(2, 3, 5);
	hold('on');
	nice_ylim5 = [inf, -inf];
	subplot(2, 3, 3);
	hold('on');
	nice_ylim3 = [inf, -inf];
	subplot(2, 3, 6);
	hold('on');
	nice_ylim6 = [inf, -inf];
	colors = jet(N);
	
	for n = m : m + N - 1
		everything = obl_open_everything(path, c, m, n);
		
		subplot(2, 3, 1);
		plot(everything.S1.eta, everything.S1.S1, 'color', colors(n - m + 1, 1 : 3));
		nice_ylim = calculate_nice_ylim(everything.S1.S1);
		nice_ylim1(1) = min([nice_ylim1(1), nice_ylim(1)]);
		nice_ylim1(2) = max([nice_ylim1(2), nice_ylim(2)]);
		
		subplot(2, 3, 4);
		plot(everything.S1.eta, everything.S1.S1p, 'color', colors(n - m + 1, 1 : 3));
		nice_ylim = calculate_nice_ylim(everything.S1.S1p);
		nice_ylim4(1) = min([nice_ylim4(1), nice_ylim(1)]);
		nice_ylim4(2) = max([nice_ylim4(2), nice_ylim(2)]);
		
		subplot(2, 3, 2);
		plot(everything.R.xi, everything.R.R1, 'color', colors(n - m + 1, 1 : 3));
		nice_ylim = calculate_nice_ylim(everything.R.R1);
		nice_ylim2(1) = min([nice_ylim2(1), nice_ylim(1)]);
		nice_ylim2(2) = max([nice_ylim2(2), nice_ylim(2)]);
		
		subplot(2, 3, 5);
		plot(everything.R.xi, everything.R.R1p, 'color', colors(n - m + 1, 1 : 3));
		nice_ylim = calculate_nice_ylim(everything.R.R1p);
		nice_ylim5(1) = min([nice_ylim5(1), nice_ylim(1)]);
		nice_ylim5(2) = max([nice_ylim5(2), nice_ylim(2)]);
		
		subplot(2, 3, 3);
		plot(everything.R.xi, everything.R.R2, 'color', colors(n - m + 1, 1 : 3));
		nice_ylim = calculate_nice_ylim(everything.R.R2);
		nice_ylim3(1) = min([nice_ylim3(1), nice_ylim(1)]);
		nice_ylim3(2) = max([nice_ylim3(2), nice_ylim(2)]);
		
		subplot(2, 3, 6);
		plot(everything.R.xi, everything.R.R2p, 'color', colors(n - m + 1, 1 : 3));
		nice_ylim = calculate_nice_ylim(everything.R.R2p);
		nice_ylim6(1) = min([nice_ylim6(1), nice_ylim(1)]);
		nice_ylim6(2) = max([nice_ylim6(2), nice_ylim(2)]);
	end
	
	subplot(2, 3, 1);
	xlim([-1.0, 1.0]);
	set(gca, 'xtick', [-1.0, -0.5, 0.0, 0.5, 1.0]);
	ylim(nice_ylim1 + 0.13 * [-1.0, 1.0] * (nice_ylim1(2) - nice_ylim1(1)));
	draw_border();
	set(gca, 'fontsize', 22);
	xlabel('\eta')
	title('S_{mn}^{(1)}(c, \eta) / {N_{mn}(c)}^{1/2}');
	set(gca, 'fontsize', 16);
	
	subplot(2, 3, 4);
	xlim([-1.0, 1.0]);
	set(gca, 'xtick', [-1.0, -0.5, 0.0, 0.5, 1.0]);
	ylim(nice_ylim4 + 0.13 * [-1.0, 1.0] * (nice_ylim4(2) - nice_ylim4(1)));
	draw_border();
	set(gca, 'fontsize', 22);
	xlabel('\eta')
	title('S_{mn}^{(1)\prime}(c, \eta) / {N_{mn}(c)}^{1/2}');
	set(gca, 'fontsize', 16);
	
	subplot(2, 3, 2);
	xlim([0.0, max_xi]);
	set(gca, 'xtick', [0.0, 0.25 * max_xi, max_xi / 2.0, 0.75 * max_xi, max_xi]);
	ylim(nice_ylim2 + 0.13 * [-1.0, 1.0] * (nice_ylim2(2) - nice_ylim2(1)));
	draw_border();
	set(gca, 'fontsize', 22);
	xlabel('\xi')
	title('R_{mn}^{(1)}(c, \xi)');
	set(gca, 'fontsize', 16);
	
	subplot(2, 3, 5);
	xlim([0.0, max_xi]);
	set(gca, 'xtick', [0.0, 0.25 * max_xi, max_xi / 2.0, 0.75 * max_xi, max_xi]);
	ylim(nice_ylim5 + 0.13 * [-1.0, 1.0] * (nice_ylim5(2) - nice_ylim5(1)));
	draw_border();
	set(gca, 'fontsize', 22);
	xlabel('\xi')
	title('R_{mn}^{(1)\prime}(c, \xi)');
	set(gca, 'fontsize', 16);
	
	subplot(2, 3, 3);
	xlim([0.0, max_xi]);
	set(gca, 'xtick', [0.0, 0.25 * max_xi, max_xi / 2.0, 0.75 * max_xi, max_xi]);
	ylim(nice_ylim3 + 0.13 * [-1.0, 1.0] * (nice_ylim3(2) - nice_ylim3(1)));
	draw_border();
	set(gca, 'fontsize', 22);
	xlabel('\xi')
	title('R_{mn}^{(2)}(c, \xi)');
	set(gca, 'fontsize', 16);
	
	subplot(2, 3, 6);
	xlim([0.0, max_xi]);
	set(gca, 'xtick', [0.0, 0.25 * max_xi, max_xi / 2.0, 0.75 * max_xi, max_xi]);
	ylim(nice_ylim6 + 0.13 * [-1.0, 1.0] * (nice_ylim6(2) - nice_ylim6(1)));
	draw_border();
	set(gca, 'fontsize', 22);
	xlabel('\xi')
	title('R_{mn}^{(2)\prime}(c, \xi)');
	set(gca, 'fontsize', 16);
	
	position = [100, 100, 3 * axis_width + 40, 2 * axis_height + 70];
	set(gcf, 'position', position);
	spaceplots(gcf, [10, 10, 50, 10] ./ position([3, 3, 4, 4]), [10, 10] ./ position([3, 4]));
	ha = axes('parent', gcf, 'units', 'normalized', 'position', [0.0, 0.0, 1.0, 1.0], 'visible', 'off');
	uistack(ha, 'bottom');
	text('parent', ha, 'horizontalalignment', 'center', 'verticalalignment', 'bottom', 'units', 'pixels', 'position', [position(3) / 2, position(4) - 40], 'string', description, 'fontsize', 22);
	set(gcf, 'color', [1.0, 1.0, 1.0]);
end
