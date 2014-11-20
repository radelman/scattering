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

%
% Calculate and save all of the various coefficients needed to calculate the
% prolate spheroidal wave functions.  This function takes care of calling
% pro_sphwv for you.  Because the output from pro_sphwv is in ASCII and, thus,
% takes up a lot of disk space, this function ZIPs up everything at the end.
%
% Arguments:
%     path - the directory in which pro_sphwv is located
%     max_memory - the maximum amount of memory, in MB, that pro_sphwv can use
%                  before automatically terminating
%     precision - the number of bits of precision pro_sphwv should use
%     c - k * a, where k is the wavenumber and 2a is the interfocal distance
%     m - one modal value
%     n - the other modal value, which is equal to m, m + 1, ...
%     n_dr - the minimum number of dr coefficients to calculate
%     log10_dr_min - have pro_sphwv keep calculating more and more dr
%                    coefficients until the log10 of their absolute magnitude
%                    drops below this value
%     n_dr_neg, log10_dr_neg_min - the same as above, but for the dr_neg
%                                  coefficients
%     n_c2k, log10_c2k_min - the same as above, but for the c2k coefficients
% Return Values:
%     None.
%
function pro_calculate_everything(path, max_memory, precision, c, m, n, n_dr, log10_dr_min, n_dr_neg, log10_dr_neg_min, n_c2k, log10_c2k_min)
	fprintf('calculating lambda_approx...\n');
	pro_calculate_lambdamn_approx(c, m, n);
	fprintf('calculating lambda, dr, dr_neg, N, F, k1, k2, and c2k...\n');
	command = sprintf('"%s/pro_sphwv" -max_memory %d -precision %d -verbose y -c %s -m %d -n %d -w everything -n_dr %d -dr_min 1.0e%d -n_dr_neg %d -dr_neg_min 1.0e%d -n_c2k %d -c2k_min 1.0e%d', path, max_memory, precision, nice_number(c, 20), m, n, n_dr, log10_dr_min, n_dr_neg, log10_dr_neg_min, n_c2k, log10_c2k_min);
	fprintf('%s\n', command);
	[ ...
	status, result ...
	] = system(command, '-echo');
	if (status == 1)
		return;
	end
	files = {sprintf('pro_%s_lambda_approx.txt', generate_name(c, m, n)), ...
	         sprintf('pro_%s_lambda.txt', generate_name(c, m, n)), ...
	         sprintf('pro_%s_log_abs_lambda.txt', generate_name(c, m, n)), ...
	         sprintf('pro_%s_dr.txt', generate_name(c, m, n)), ...
	         sprintf('pro_%s_log_abs_dr.txt', generate_name(c, m, n)), ...
	         sprintf('pro_%s_dr_neg.txt', generate_name(c, m, n)), ...
	         sprintf('pro_%s_log_abs_dr_neg.txt', generate_name(c, m, n)), ...
	         sprintf('pro_%s_N.txt', generate_name(c, m, n)), ...
	         sprintf('pro_%s_log_abs_N.txt', generate_name(c, m, n)), ...
	         sprintf('pro_%s_F.txt', generate_name(c, m, n)), ...
	         sprintf('pro_%s_log_abs_F.txt', generate_name(c, m, n)), ...
	         sprintf('pro_%s_k1.txt', generate_name(c, m, n)), ...
	         sprintf('pro_%s_log_abs_k1.txt', generate_name(c, m, n)), ...
	         sprintf('pro_%s_k2.txt', generate_name(c, m, n)), ...
	         sprintf('pro_%s_log_abs_k2.txt', generate_name(c, m, n)), ...
	         sprintf('pro_%s_c2k.txt', generate_name(c, m, n)), ...
	         sprintf('pro_%s_log_abs_c2k.txt', generate_name(c, m, n))};
	zip(sprintf('data/pro_%s.zip', generate_name(c, m, n)), files, 'data');
	delete(sprintf('data/pro_%s_*.txt', generate_name(c, m, n)));
end
