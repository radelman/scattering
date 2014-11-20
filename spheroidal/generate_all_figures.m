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

function generate_all_figures()
	pro_awesome1('saved', 10.0, 30, 30, ', c = 10', [0.0, 1000.0, 2000.0, 3000.0], [0.0, 25.0, 50.0, 75.0, 100.0], [0.0, 25.0, 50.0, 75.0, 100.0], [-20.0, -10.0, 0.0, 10.0, 20.0, 30.0]);
	export_fig(1, 'images/pro_lambda_00010000.pdf');
	export_fig(2, 'images/pro_N_00010000.pdf');
	export_fig(3, 'images/pro_k1_00010000.pdf');
	export_fig(4, 'images/pro_k2_00010000.pdf');
	close('all');
	
	pro_awesome3('saved', 10.0, 10, 30, 9.0, 'c = 10, m = 10, n = 10, 11, ..., 39 (Blue to Red)', 350, 325);
	export_fig(1, 'images/pro_awesome3_00010000_010.pdf');
	close('all');
	
	obl_awesome1('saved', 10.0, 30, 30, ', c = 10', [0.0, 1000.0, 2000.0, 3000.0], [0.0, 25.0, 50.0, 75.0, 100.0], [0.0, 25.0, 50.0, 75.0, 100.0], [-20.0, -10.0, 0.0, 10.0, 20.0, 30.0], [0.0, 25.0, 50.0, 75.0, 100.0]);
	export_fig(1, 'images/obl_lambda_00010000.pdf');
	export_fig(2, 'images/obl_N_00010000.pdf');
	export_fig(3, 'images/obl_k1_00010000.pdf');
	export_fig(4, 'images/obl_k2_00010000.pdf');
	export_fig(5, 'images/obl_Q_00010000.pdf');
	close('all');
	
	obl_awesome3('saved', 10.0, 10, 30, 8.0, 'c = 10, m = 10, n = 10, 11, ..., 39 (Blue to Red)', 350, 325);
	export_fig(1, 'images/obl_awesome3_00010000_010.pdf');
	close('all');
	
	everything = pro_open_everything('saved', 10.0, 10, 39);
	pro_check_everything(everything, ' for Prolate Case');
	figure(4);
	ylim([-100.0, 0.0]);
	title('Error of Wronskian for Prolate Case');
	export_fig(4, 'images/pro_different_W_log10_abs_errors_00010000_010_039.pdf');
	close('all');
	
	everything = obl_open_everything('saved', 10.0, 10, 39);
	obl_check_everything(everything, ' for Oblate Case');
	figure(5);
	ylim([-100.0, 0.0]);
	title('Error of Wronskian for Oblate Case');
	export_fig(5, 'images/obl_different_W_log10_abs_errors_00010000_010_039.pdf');
	close('all');
	
	obl_awesome2('saved', 25.0, 49, 50, 51, 50, 101, 101, ', c = 10, m = 10', [-40.0, -30.0, -20.0, -10.0, 0.0], [-50.0, 0.0, 50.0, 100.0], [-50.0, 0.0, 50.0, 100.0]);
	export_fig(3, 'images/obl_B2r_00025000_049.pdf');
	close('all');
end
