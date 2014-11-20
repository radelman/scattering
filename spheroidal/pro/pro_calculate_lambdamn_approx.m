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
% Calculate and save an approximate value of the characteristic value.
%
% Arguments:
%     c - k * a, where k is the wavenumber and 2a is the interfocal distance
%     m - one modal value
%     n - the other modal value, which is equal to m, m + 1, ...
% Return Values:
%     None.
%
function pro_calculate_lambdamn_approx(c, m, n)
	N = m + n + 200;
	A = zeros(N, N);
	if (mod(n - m, 2) == 0)
		r = 0;
	else
		r = 1;
	end
	for i = 1 : N
		if (i == 1)
			A(1, 1) = calculate_betar(c, m, r);
			A(1, 2) = calculate_alphar(c, m, r);
		elseif (i == N)
			A(N, N - 1) = calculate_gammar(c, m, r);
			A(N, N) = calculate_betar(c, m, r);
		else
			A(i, i - 1) = calculate_gammar(c, m, r);
			A(i, i) = calculate_betar(c, m, r);
			A(i, i + 1) = calculate_alphar(c, m, r);
		end
		r = r + 2;
	end
	d = eig(A);
	if (mod(n - m, 2) == 0)
		lambda_approx = d((n - m + 2) / 2);
	else
		lambda_approx = d((n - m + 1) / 2);
	end
	fid = fopen(sprintf('data/pro_%s_lambda_approx.txt', generate_name(c, m, n)), 'w');
	fprintf(fid, '%s\n', nice_number(lambda_approx, 20));
	fclose(fid);
end

function alpha = calculate_alphar(c, m, r)
	alpha = (((2 * m + r + 2) * (2 * m + r + 1)) / ((2 * m + 2 * r + 5) * (2 * m + 2 * r + 3))) * (c ^ 2);
end

function beta = calculate_betar(c, m, r)
	beta = (m + r) * (m + r + 1) + ((2 * (m + r) * (m + r + 1) - 2 * (m ^ 2) - 1) / ((2 * m + 2 * r - 1) * (2 * m + 2 * r + 3))) * (c ^ 2);
end

function gamma = calculate_gammar(c, m, r)
	gamma = ((r * (r - 1)) / ((2 * m + 2 * r - 3) * (2 * m + 2 * r - 1))) * (c ^ 2);
end
