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

function obl_check_everything(everything, description)
	figure();
	plot(everything.S1.eta, everything.S1.S1_log_abs_difference / log(10.0), 'color', [0.0, 0.0, 0.0], 'linewidth', 2.0);
	hold('on');
	plot(everything.S1.eta, everything.S1.S1p_log_abs_difference / log(10.0), 'color', [1.0, 0.0, 0.0], 'linewidth', 2.0);
	xlim([-1.0, 1.0]);
	set(gca, 'xtick', [-1.0, -0.5, 0.0, 0.5, 1.0]);
	draw_border();
	set(gca, 'fontsize', 20);
	xlabel('\eta');
	ylabel('log_{10}(|S1\_1 - S1\_2|)');
	title(sprintf('Difference Between Methods S1\\_1 and S1\\_2%s', description));
	legendflex_pre2014b({'S_{mn}^{(1)}(-ic, \eta)', 'S_{mn}^{(1)\prime}(-ic, \eta)'}, 'ref', gca, 'anchor', [5, 5], 'buffer', [-30, 30]);
	set(gcf, 'position', [100, 100, 800, 600]);
	set(gcf, 'color', [1.0, 1.0, 1.0]);
	
	figure();
	plot(everything.R.xi, everything.R.R1_log_abs_difference / log(10.0), 'color', [0.0, 0.0, 0.0], 'linewidth', 2.0);
	hold('on');
	plot(everything.R.xi, everything.R.R1p_log_abs_difference / log(10.0), 'color', [1.0, 0.0, 0.0], 'linewidth', 2.0);
	xlim([0.0, max(everything.R.xi)]);
	set(gca, 'xtick', [0.0, 0.25 * max(everything.R.xi), max(everything.R.xi) / 2.0, 0.75 * max(everything.R.xi), max(everything.R.xi)]);
	draw_border();
	set(gca, 'fontsize', 20);
	xlabel('\xi');
	ylabel('log_{10}(|R1\_1 - R1\_2|)');
	title(sprintf('Difference Between Methods R1\\_1 and R1\\_2%s', description));
	legendflex_pre2014b({'R_{mn}^{(1)}(-ic, i\xi)', 'R_{mn}^{(1)\prime}(-ic, i\xi)'}, 'ref', gca, 'anchor', [1, 1], 'buffer', [30, -30]);
	set(gcf, 'position', [100, 100, 800, 600]);
	set(gcf, 'color', [1.0, 1.0, 1.0]);
	
	figure();
	plot(everything.R.xi, everything.R.R2_log_abs_difference_1_2 / log(10.0), 'color', [0.0, 0.0, 0.0], 'linewidth', 2.0);
	hold('on');
	plot(everything.R.xi, everything.R.R2_log_abs_difference_1_3 / log(10.0), 'color', [1.0, 0.0, 0.0], 'linewidth', 2.0);
	plot(everything.R.xi, everything.R.R2_log_abs_difference_2_3 / log(10.0), 'color', [0.0, 0.0, 1.0], 'linewidth', 2.0);
	xlim([0.0, max(everything.R.xi)]);
	set(gca, 'xtick', [0.0, 0.25 * max(everything.R.xi), max(everything.R.xi) / 2.0, 0.75 * max(everything.R.xi), max(everything.R.xi)]);
	draw_border();
	set(gca, 'fontsize', 20);
	xlabel('\xi');
	title(sprintf('Difference Between Methods for R_{mn}^{(2)}(-ic, i\\xi)%s', description));
	legendflex_pre2014b({'log_{10}(|R1\_1 - R1\_2|)', 'log_{10}(|R1\_1 - R1\_3|)', 'log_{10}(|R1\_2 - R1\_3|)'}, 'ref', gca, 'anchor', [3, 3], 'buffer', [-30, -30]);
	set(gcf, 'position', [100, 100, 800, 600]);
	set(gcf, 'color', [1.0, 1.0, 1.0]);
	
	figure();
	plot(everything.R.xi, everything.R.R2p_log_abs_difference_1_2 / log(10.0), 'color', [0.0, 0.0, 0.0], 'linewidth', 2.0);
	hold('on');
	plot(everything.R.xi, everything.R.R2p_log_abs_difference_1_3 / log(10.0), 'color', [1.0, 0.0, 0.0], 'linewidth', 2.0);
	plot(everything.R.xi, everything.R.R2p_log_abs_difference_2_3 / log(10.0), 'color', [0.0, 0.0, 1.0], 'linewidth', 2.0);
	xlim([0.0, max(everything.R.xi)]);
	set(gca, 'xtick', [0.0, 0.25 * max(everything.R.xi), max(everything.R.xi) / 2.0, 0.75 * max(everything.R.xi), max(everything.R.xi)]);
	draw_border();
	set(gca, 'fontsize', 20);
	xlabel('\xi');
	title(sprintf('Difference Between Methods for R_{mn}^{(2)\\prime}(-ic, i\\xi)%s', description));
	legendflex_pre2014b({'log_{10}(|R1\_1 - R1\_2|)', 'log_{10}(|R1\_1 - R1\_3|)', 'log_{10}(|R1\_2 - R1\_3|)'}, 'ref', gca, 'anchor', [3, 3], 'buffer', [-30, -30]);
	set(gcf, 'position', [100, 100, 800, 600]);
	set(gcf, 'color', [1.0, 1.0, 1.0]);
	
	figure();
	plot(everything.R.xi, everything.R.W_1_1_log_abs_error / log(10.0) - everything.R.log_W / log(10.0), 'color', [0.0, 0.0, 0.0], 'linewidth', 2.0);
	hold('on');
	plot(everything.R.xi, everything.R.W_1_2_log_abs_error / log(10.0) - everything.R.log_W / log(10.0), 'color', [1.0, 0.0, 0.0], 'linewidth', 2.0);
	plot(everything.R.xi, everything.R.W_1_31_log_abs_error / log(10.0) - everything.R.log_W / log(10.0), 'color', [0.0, 0.7, 0.0], 'linewidth', 2.0);
	plot(everything.R.xi, everything.R.W_2_1_log_abs_error / log(10.0) - everything.R.log_W / log(10.0), 'color', [0.0, 0.0, 1.0], 'linewidth', 2.0);
	plot(everything.R.xi, everything.R.W_2_2_log_abs_error / log(10.0) - everything.R.log_W / log(10.0), 'color', [1.0, 0.0, 1.0], 'linewidth', 2.0);
	plot(everything.R.xi, everything.R.W_2_32_log_abs_error / log(10.0) - everything.R.log_W / log(10.0), 'color', [0.0, 1.0, 1.0], 'linewidth', 2.0);
	xlim([0.0, max(everything.R.xi)]);
	set(gca, 'xtick', [0.0, 0.25 * max(everything.R.xi), max(everything.R.xi) / 2.0, 0.75 * max(everything.R.xi), max(everything.R.xi)]);
	draw_border();
	set(gca, 'fontsize', 20);
	xlabel('\xi');
	ylabel('log_{10}(|Relative Error|)');
	title(sprintf('Relative Error of Wronskian Using Different Methods%s', description));
	legendflex_pre2014b({'R1\_1, R2\_1', 'R1\_1, R2\_2', 'R1\_1, R2\_31', 'R1\_2, R2\_1', 'R1\_2, R2\_2', 'R1\_2, R2\_32'}, 'ref', gca, 'anchor', [3, 3], 'buffer', [-30, -30]);
	set(gcf, 'position', [100, 100, 800, 600]);
	set(gcf, 'color', [1.0, 1.0, 1.0]);
	
	figure();
	plot(everything.S1.eta, everything.S1.S1, 'color', [0.0, 0.0, 0.0], 'linewidth', 2.0);
	hold('on');
	xlim([-1.0, 1.0]);
	set(gca, 'xtick', [-1.0, -0.5, 0.0, 0.5, 1.0]);
	nice_ylim = calculate_nice_ylim(everything.S1.S1);
	ylim(nice_ylim + 0.1 * [-1.0, 1.0] * (nice_ylim(2) - nice_ylim(1)));
	draw_border();
	set(gca, 'fontsize', 20);
	xlabel('\eta');
	title(sprintf('Oblate Angle Functions%s', description));
	legendflex_pre2014b({'S_{mn}^{(1)}(-ic, \eta) / {N_{mn}(-ic)}^{1/2}'}, 'ref', gca, 'anchor', [3, 3], 'buffer', [-30, -30]);
	set(gcf, 'position', [100, 100, 800, 600]);
	set(gcf, 'color', [1.0, 1.0, 1.0]);
	
	figure();
	plot(everything.S1.eta, everything.S1.S1p, 'color', [0.0, 0.0, 0.0], 'linewidth', 2.0);
	hold('on');
	xlim([-1.0, 1.0]);
	set(gca, 'xtick', [-1.0, -0.5, 0.0, 0.5, 1.0]);
	nice_ylim = calculate_nice_ylim(everything.S1.S1p);
	ylim(nice_ylim + 0.1 * [-1.0, 1.0] * (nice_ylim(2) - nice_ylim(1)));
	draw_border();
	set(gca, 'fontsize', 20);
	xlabel('\eta');
	title(sprintf('Derivative of Oblate Angle Functions%s', description));
	legendflex_pre2014b({'S_{mn}^{(1)\prime}(-ic, \eta) / {N_{mn}(-ic)}^{1/2}'}, 'ref', gca, 'anchor', [3, 3], 'buffer', [-30, -30]);
	set(gcf, 'position', [100, 100, 800, 600]);
	set(gcf, 'color', [1.0, 1.0, 1.0]);
	
	figure();
	plot(everything.R.xi, everything.R.R1, 'color', [0.0, 0.0, 0.0], 'linewidth', 2.0);
	hold('on');
	plot(everything.R.xi, everything.R.R2, 'color', [1.0, 0.0, 0.0], 'linewidth', 2.0);
	xlim([0.0, max(everything.R.xi)]);
	set(gca, 'xtick', [0.0, 0.25 * max(everything.R.xi), max(everything.R.xi) / 2.0, 0.75 * max(everything.R.xi), max(everything.R.xi)]);
	nice_ylim1 = calculate_nice_ylim(everything.R.R1);
	nice_ylim2 = calculate_nice_ylim(everything.R.R2);
	nice_ylim = [min([nice_ylim1(1), nice_ylim2(1)]), max([nice_ylim1(2), nice_ylim2(2)])];
	ylim(nice_ylim + 0.1 * [-1.0, 1.0] * (nice_ylim(2) - nice_ylim(1)));
	draw_border();
	set(gca, 'fontsize', 20);
	xlabel('\xi');
	title(sprintf('Oblate Radial Functions%s', description));
	legendflex_pre2014b({'R_{mn}^{(1)}(-ic, i\xi)', 'R_{mn}^{(2)}(-ic, i\xi)'}, 'ref', gca, 'anchor', [3, 3], 'buffer', [-30, -30]);
	set(gcf, 'position', [100, 100, 800, 600]);
	set(gcf, 'color', [1.0, 1.0, 1.0]);
	
	figure();
	plot(everything.R.xi, everything.R.R1p, 'color', [0.0, 0.0, 0.0], 'linewidth', 2.0);
	hold('on');
	plot(everything.R.xi, everything.R.R2p, 'color', [1.0, 0.0, 0.0], 'linewidth', 2.0);
	xlim([0.0, max(everything.R.xi)]);
	set(gca, 'xtick', [0.0, 0.25 * max(everything.R.xi), max(everything.R.xi) / 2.0, 0.75 * max(everything.R.xi), max(everything.R.xi)]);
	nice_ylim1 = calculate_nice_ylim(everything.R.R1p);
	nice_ylim2 = calculate_nice_ylim(everything.R.R2p);
	nice_ylim = [min([nice_ylim1(1), nice_ylim2(1)]), max([nice_ylim1(2), nice_ylim2(2)])];
	ylim(nice_ylim + 0.1 * [-1.0, 1.0] * (nice_ylim(2) - nice_ylim(1)));
	draw_border();
	set(gca, 'fontsize', 20);
	xlabel('\xi');
	title(sprintf('Derivative of Oblate Radial Functions%s', description));
	legendflex_pre2014b({'R_{mn}^{(1)\prime}(-ic, i\xi)', 'R_{mn}^{(2)\prime}(-ic, i\xi)'}, 'ref', gca, 'anchor', [3, 3], 'buffer', [-30, -30]);
	set(gcf, 'position', [100, 100, 800, 600]);
	set(gcf, 'color', [1.0, 1.0, 1.0]);
	
	figure();
	plot(everything.R.xi, everything.R.W_log_abs_error / log(10.0) - everything.R.log_W / log(10.0), 'color', [0.0, 0.0, 0.0], 'linewidth', 2.0);
	hold('on');
	xlim([0.0, max(everything.R.xi)]);
	set(gca, 'xtick', [0.0, 0.25 * max(everything.R.xi), max(everything.R.xi) / 2.0, 0.75 * max(everything.R.xi), max(everything.R.xi)]);
	draw_border();
	set(gca, 'fontsize', 20);
	xlabel('\xi');
	ylabel('log_{10}(|Relative Error|)');
	title(sprintf('Relative Error of Wronskian Using Best Combination of Methods%s', description));
	set(gcf, 'position', [100, 100, 800, 600]);
	set(gcf, 'color', [1.0, 1.0, 1.0]);
end
