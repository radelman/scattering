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
% Unzip all of the ZIP'd-up prolate spheroidal coefficients and wave functions,
% and organize them into a nice, easily-accessible MATLAB struct.  This is
% where the best combination of prolate spheroidal radial functions is selected
% to minimize the error in the Wronskian.
%
% Arguments:
%     path - the directory in which the prolate spheroidal coefficients and
%            wave functions have been saved and ZIP'd up
%     c - k * a, where k is the wavenumber and 2a is the interfocal distance
%     m - one modal value
%     n - the other modal value, which is equal to m, m + 1, ...
%     name - the same name from pro_calculate_functions
% Return Values:
%     everything - a MATLAB struct with all of the prolate spheroidal
%                  coefficients and wave functions
%
function everything = pro_load_everything(path, c, m, n, name)
	everything = struct();
	everything.c = c;
	everything.m = m;
	everything.n = n;
	try
		everything.coefficients = struct();
		unzip(sprintf('%s/pro_%s.zip', path, generate_name(c, m, n)), path);
		everything.coefficients.lambda = load(sprintf('%s/pro_%s_lambda.txt', path, generate_name(c, m, n)));
		everything.coefficients.log_abs_lambda = load(sprintf('%s/pro_%s_log_abs_lambda.txt', path, generate_name(c, m, n)));
		everything.coefficients.dr = load(sprintf('%s/pro_%s_dr.txt', path, generate_name(c, m, n))).';
		everything.coefficients.log_abs_dr = load(sprintf('%s/pro_%s_log_abs_dr.txt', path, generate_name(c, m, n))).';
		everything.coefficients.dr_neg = load(sprintf('%s/pro_%s_dr_neg.txt', path, generate_name(c, m, n))).';
		everything.coefficients.log_abs_dr_neg = load(sprintf('%s/pro_%s_log_abs_dr_neg.txt', path, generate_name(c, m, n))).';
		everything.coefficients.N = load(sprintf('%s/pro_%s_N.txt', path, generate_name(c, m, n)));
		everything.coefficients.log_abs_N = load(sprintf('%s/pro_%s_log_abs_N.txt', path, generate_name(c, m, n)));
		everything.coefficients.F = load(sprintf('%s/pro_%s_F.txt', path, generate_name(c, m, n)));
		everything.coefficients.log_abs_F = load(sprintf('%s/pro_%s_log_abs_F.txt', path, generate_name(c, m, n)));
		everything.coefficients.k1 = load(sprintf('%s/pro_%s_k1.txt', path, generate_name(c, m, n)));
		everything.coefficients.log_abs_k1 = load(sprintf('%s/pro_%s_log_abs_k1.txt', path, generate_name(c, m, n)));
		everything.coefficients.k2 = load(sprintf('%s/pro_%s_k2.txt', path, generate_name(c, m, n)));
		everything.coefficients.log_abs_k2 = load(sprintf('%s/pro_%s_log_abs_k2.txt', path, generate_name(c, m, n)));
		everything.coefficients.c2k = load(sprintf('%s/pro_%s_c2k.txt', path, generate_name(c, m, n))).';
		everything.coefficients.log_abs_c2k = load(sprintf('%s/pro_%s_log_abs_c2k.txt', path, generate_name(c, m, n))).';
	catch
		fprintf('can''t open one or more of the coefficients...\n');
	end
	try
		everything.S1 = struct();
		unzip(sprintf('%s/pro_%s_%s.zip', path, generate_name(c, m, n), name), path);
		temp = load(sprintf('%s/pro_%s_S1.txt', path, generate_name(c, m, n))).';
		everything.S1.eta = temp(2, :);
		[ ...
		everything.S1.eta, ix ...
		] = sort(everything.S1.eta);
		everything.S1.S1_1 = temp(3, ix);
		everything.S1.S1p_1 = temp(4, ix);
		everything.S1.S1_2 = temp(5, ix);
		everything.S1.S1p_2 = temp(6, ix);
		everything.S1.S1_log_abs_difference = temp(7, ix);
		everything.S1.S1p_log_abs_difference = temp(8, ix);
		everything.S1.S1 = everything.S1.S1_1;
		everything.S1.S1p = everything.S1.S1p_1;
	catch
		fprintf('can''t open S1...\n');
	end
	try
		everything.R = pro_load_everything_R(sprintf('%s/pro_%s_R.txt', path, generate_name(c, m, n)));
	catch
		fprintf('can''t open R...\n');
	end
	try
		everything.R_small = pro_load_everything_R(sprintf('%s/pro_%s_R_small.txt', path, generate_name(c, m, n)));
	catch
		fprintf('can''t open R_small...\n');
	end
	delete(sprintf('%s/pro_%s_*.txt', path, generate_name(c, m, n)));
end

function R = pro_load_everything_R(name)
	R = struct();
	temp = load(name).';
	R.xi = temp(2, :);
	R.log_xi = temp(3, :);
	R.R1_1 = temp(4, :);
	R.R1p_1 = temp(5, :);
	R.R1_2 = temp(6, :);
	R.R1p_2 = temp(7, :);
	R.R1_log_abs_difference = temp(8, :);
	R.R1p_log_abs_difference = temp(9, :);
	R.R2_1 = temp(10, :);
	R.R2p_1 = temp(11, :);
	R.R2_2 = temp(12, :);
	R.R2p_2 = temp(13, :);
	R.R2_log_abs_difference = temp(14, :);
	R.R2p_log_abs_difference = temp(15, :);
	R.W = temp(16, :);
	R.log_W = temp(17, :);
	R.W_1_1_log_abs_error = temp(18, :);
	R.W_1_2_log_abs_error = temp(19, :);
	R.W_2_1_log_abs_error = temp(20, :);
	R.W_2_2_log_abs_error = temp(21, :);
	R.R1 = zeros(1, length(R.xi));
	R.R1p = zeros(1, length(R.xi));
	R.R2 = zeros(1, length(R.xi));
	R.R2p = zeros(1, length(R.xi));
	R.W_log_abs_error = zeros(1, length(R.xi));
	for i = 1 : length(R.xi)
		if (R.xi(i) > 1.0)
			[ ...
			R.W_log_abs_error(i), idx ...
			] = min([R.W_1_1_log_abs_error(i), R.W_1_2_log_abs_error(i), R.W_2_1_log_abs_error(i), R.W_2_2_log_abs_error(i)]);
		else
			idx = 5;
		end
		if (idx == 1)
			R.R1(i) = R.R1_1(i);
			R.R1p(i) = R.R1p_1(i);
			R.R2(i) = R.R2_1(i);
			R.R2p(i) = R.R2p_1(i);
		elseif (idx == 2)
			R.R1(i) = R.R1_1(i);
			R.R1p(i) = R.R1p_1(i);
			R.R2(i) = R.R2_2(i);
			R.R2p(i) = R.R2p_2(i);
		elseif (idx == 3)
			R.R1(i) = R.R1_2(i);
			R.R1p(i) = R.R1p_2(i);
			R.R2(i) = R.R2_1(i);
			R.R2p(i) = R.R2p_1(i);
		elseif (idx == 4)
			R.R1(i) = R.R1_2(i);
			R.R1p(i) = R.R1p_2(i);
			R.R2(i) = R.R2_2(i);
			R.R2p(i) = R.R2p_2(i);
		else
			R.R1(i) = R.R1_1(i);
			R.R1p(i) = R.R1p_1(i);
			R.R2(i) = sign(R.R2_2(i + 1)) * inf;
			R.R2p(i) = sign(R.R2p_2(i + 1)) * inf;
			R.W_log_abs_error(i) = nan;
		end
	end
end
