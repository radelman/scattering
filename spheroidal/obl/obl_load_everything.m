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
% Unzip all of the ZIP'd-up oblate spheroidal coefficients and wave functions,
% and organize them into a nice, easily-accessible MATLAB struct.  This is
% where the best combination of oblate spheroidal radial functions is selected
% to minimize the error in the Wronskian.
%
% Arguments:
%     path - the directory in which the oblate spheroidal coefficients and
%            wave functions have been saved and ZIP'd up
%     c - k * a, where k is the wavenumber and 2a is the interfocal distance
%     m - one modal value
%     n - the other modal value, which is equal to m, m + 1, ...
%     name - the same name from pro_calculate_functions
% Return Values:
%     everything - a MATLAB struct with all of the oblate spheroidal
%                  coefficients and wave functions
%
function everything = obl_load_everything(path, c, m, n, name)
	everything = struct();
	everything.c = c;
	everything.m = m;
	everything.n = n;
	try
		everything.coefficients = struct();
		unzip(sprintf('%s/obl_%s.zip', path, generate_name(c, m, n)), path);
		everything.coefficients.lambda = load(sprintf('%s/obl_%s_lambda.txt', path, generate_name(c, m, n)));
		everything.coefficients.log_abs_lambda = load(sprintf('%s/obl_%s_log_abs_lambda.txt', path, generate_name(c, m, n)));
		everything.coefficients.dr = load(sprintf('%s/obl_%s_dr.txt', path, generate_name(c, m, n))).';
		everything.coefficients.log_abs_dr = load(sprintf('%s/obl_%s_log_abs_dr.txt', path, generate_name(c, m, n))).';
		everything.coefficients.dr_neg = load(sprintf('%s/obl_%s_dr_neg.txt', path, generate_name(c, m, n))).';
		everything.coefficients.log_abs_dr_neg = load(sprintf('%s/obl_%s_log_abs_dr_neg.txt', path, generate_name(c, m, n))).';
		everything.coefficients.N = load(sprintf('%s/obl_%s_N.txt', path, generate_name(c, m, n)));
		everything.coefficients.log_abs_N = load(sprintf('%s/obl_%s_log_abs_N.txt', path, generate_name(c, m, n)));
		everything.coefficients.F = load(sprintf('%s/obl_%s_F.txt', path, generate_name(c, m, n)));
		everything.coefficients.log_abs_F = load(sprintf('%s/obl_%s_log_abs_F.txt', path, generate_name(c, m, n)));
		everything.coefficients.k1 = load(sprintf('%s/obl_%s_k1.txt', path, generate_name(c, m, n)));
		everything.coefficients.log_abs_k1 = load(sprintf('%s/obl_%s_log_abs_k1.txt', path, generate_name(c, m, n)));
		everything.coefficients.k2 = load(sprintf('%s/obl_%s_k2.txt', path, generate_name(c, m, n)));
		everything.coefficients.log_abs_k2 = load(sprintf('%s/obl_%s_log_abs_k2.txt', path, generate_name(c, m, n)));
		everything.coefficients.c2k = load(sprintf('%s/obl_%s_c2k.txt', path, generate_name(c, m, n))).';
		everything.coefficients.log_abs_c2k = load(sprintf('%s/obl_%s_log_abs_c2k.txt', path, generate_name(c, m, n))).';
		everything.coefficients.Q = load(sprintf('%s/obl_%s_Q.txt', path, generate_name(c, m, n)));
		everything.coefficients.log_abs_Q = load(sprintf('%s/obl_%s_log_abs_Q.txt', path, generate_name(c, m, n)));
		everything.coefficients.B2r = load(sprintf('%s/obl_%s_B2r.txt', path, generate_name(c, m, n))).';
		everything.coefficients.log_abs_B2r = load(sprintf('%s/obl_%s_log_abs_B2r.txt', path, generate_name(c, m, n))).';
	catch
		fprintf('can''t open one or more of the coefficients...\n');
	end
	try
		everything.S1 = struct();
		unzip(sprintf('%s/obl_%s_%s.zip', path, generate_name(c, m, n), name), path);
		temp = load(sprintf('%s/obl_%s_S1.txt', path, generate_name(c, m, n))).';
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
		everything.R = obl_load_everything_R(sprintf('%s/obl_%s_R.txt', path, generate_name(c, m, n)));
	catch
		fprintf('can''t open R...\n');
	end
	delete(sprintf('%s/obl_%s_*.txt', path, generate_name(c, m, n)));
end

function R = obl_load_everything_R(name)
	R = struct();
	temp = load(name).';
	R.xi = temp(2, :);
	R.R1_1 = temp(3, :);
	R.R1p_1 = temp(4, :);
	R.R1_2 = temp(5, :);
	R.R1p_2 = temp(6, :);
	R.R1_log_abs_difference = temp(7, :);
	R.R1p_log_abs_difference = temp(8, :);
	R.R2_1 = temp(9, :);
	R.R2p_1 = temp(10, :);
	R.R2_2 = temp(11, :);
	R.R2p_2 = temp(12, :);
	R.R2_31 = temp(13, :);
	R.R2p_31 = temp(14, :);
	R.R2_32 = temp(15, :);
	R.R2p_32 = temp(16, :);
	R.R2_log_abs_difference_1_2 = temp(17, :);
	R.R2p_log_abs_difference_1_2 = temp(18, :);
	R.R2_log_abs_difference_1_31 = temp(19, :);
	R.R2p_log_abs_difference_1_31 = temp(20, :);
	R.R2_log_abs_difference_1_32 = temp(21, :);
	R.R2p_log_abs_difference_1_32 = temp(22, :);
	R.R2_log_abs_difference_2_31 = temp(23, :);
	R.R2p_log_abs_difference_2_31 = temp(24, :);
	R.R2_log_abs_difference_2_32 = temp(25, :);
	R.R2p_log_abs_difference_2_32 = temp(26, :);
	R.R2_log_abs_difference_31_32 = temp(27, :);
	R.R2p_log_abs_difference_31_32 = temp(28, :);
	R.W = temp(29, :);
	R.log_W = temp(30, :);
	R.W_1_1_log_abs_error = temp(31, :);
	R.W_1_2_log_abs_error = temp(32, :);
	R.W_1_31_log_abs_error = temp(33, :);
	R.W_1_32_log_abs_error = temp(34, :);
	R.W_2_1_log_abs_error = temp(35, :);
	R.W_2_2_log_abs_error = temp(36, :);
	R.W_2_31_log_abs_error = temp(37, :);
	R.W_2_32_log_abs_error = temp(38, :);
	R.R1 = zeros(1, length(R.xi));
	R.R1p = zeros(1, length(R.xi));
	R.R2 = zeros(1, length(R.xi));
	R.R2p = zeros(1, length(R.xi));
	R.R2_log_abs_difference_1_3 = zeros(1, length(R.xi));
	R.R2p_log_abs_difference_1_3 = zeros(1, length(R.xi));
	R.R2_log_abs_difference_2_3 = zeros(1, length(R.xi));
	R.R2p_log_abs_difference_2_3 = zeros(1, length(R.xi));
	R.W_log_abs_error = zeros(1, length(R.xi));
	for i = 1 : length(R.xi)
		[ ...
		R.W_log_abs_error(i), idx ...
		] = min([R.W_1_1_log_abs_error(i), R.W_1_2_log_abs_error(i), R.W_1_31_log_abs_error(i), R.W_2_1_log_abs_error(i), R.W_2_2_log_abs_error(i), R.W_2_32_log_abs_error(i)]);
		if (idx == 1)
			R.R1(i) = R.R1_1(i);
			R.R1p(i) = R.R1p_1(i);
			R.R2(i) = R.R2_1(i);
			R.R2p(i) = R.R2p_1(i);
			R.R2_log_abs_difference_1_3(i) = R.R2_log_abs_difference_1_31(i);
			R.R2p_log_abs_difference_1_3(i) = R.R2p_log_abs_difference_1_31(i);
			R.R2_log_abs_difference_2_3(i) = R.R2_log_abs_difference_2_31(i);
			R.R2p_log_abs_difference_2_3(i) = R.R2p_log_abs_difference_2_31(i);
		elseif (idx == 2)
			R.R1(i) = R.R1_1(i);
			R.R1p(i) = R.R1p_1(i);
			R.R2(i) = R.R2_2(i);
			R.R2p(i) = R.R2p_2(i);
			R.R2_log_abs_difference_1_3(i) = R.R2_log_abs_difference_1_31(i);
			R.R2p_log_abs_difference_1_3(i) = R.R2p_log_abs_difference_1_31(i);
			R.R2_log_abs_difference_2_3(i) = R.R2_log_abs_difference_2_31(i);
			R.R2p_log_abs_difference_2_3(i) = R.R2p_log_abs_difference_2_31(i);
		elseif (idx == 3)
			R.R1(i) = R.R1_1(i);
			R.R1p(i) = R.R1p_1(i);
			R.R2(i) = R.R2_31(i);
			R.R2p(i) = R.R2p_31(i);
			R.R2_log_abs_difference_1_3(i) = R.R2_log_abs_difference_1_31(i);
			R.R2p_log_abs_difference_1_3(i) = R.R2p_log_abs_difference_1_31(i);
			R.R2_log_abs_difference_2_3(i) = R.R2_log_abs_difference_2_31(i);
			R.R2p_log_abs_difference_2_3(i) = R.R2p_log_abs_difference_2_31(i);
		elseif (idx == 4)
			R.R1(i) = R.R1_2(i);
			R.R1p(i) = R.R1p_2(i);
			R.R2(i) = R.R2_1(i);
			R.R2p(i) = R.R2p_1(i);
			R.R2_log_abs_difference_1_3(i) = R.R2_log_abs_difference_1_32(i);
			R.R2p_log_abs_difference_1_3(i) = R.R2p_log_abs_difference_1_32(i);
			R.R2_log_abs_difference_2_3(i) = R.R2_log_abs_difference_2_32(i);
			R.R2p_log_abs_difference_2_3(i) = R.R2p_log_abs_difference_2_32(i);
		elseif (idx == 5)
			R.R1(i) = R.R1_2(i);
			R.R1p(i) = R.R1p_2(i);
			R.R2(i) = R.R2_2(i);
			R.R2p(i) = R.R2p_2(i);
			R.R2_log_abs_difference_1_3(i) = R.R2_log_abs_difference_1_32(i);
			R.R2p_log_abs_difference_1_3(i) = R.R2p_log_abs_difference_1_32(i);
			R.R2_log_abs_difference_2_3(i) = R.R2_log_abs_difference_2_32(i);
			R.R2p_log_abs_difference_2_3(i) = R.R2p_log_abs_difference_2_32(i);
		else
			R.R1(i) = R.R1_2(i);
			R.R1p(i) = R.R1p_2(i);
			R.R2(i) = R.R2_32(i);
			R.R2p(i) = R.R2p_32(i);
			R.R2_log_abs_difference_1_3(i) = R.R2_log_abs_difference_1_32(i);
			R.R2p_log_abs_difference_1_3(i) = R.R2p_log_abs_difference_1_32(i);
			R.R2_log_abs_difference_2_3(i) = R.R2_log_abs_difference_2_32(i);
			R.R2p_log_abs_difference_2_3(i) = R.R2p_log_abs_difference_2_32(i);
		end
	end
end
