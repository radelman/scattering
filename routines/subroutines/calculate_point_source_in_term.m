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
% Calculate a term in the spheroidal wave function expansion for calculating
% the incident field due to a point source.  This code works for both the
% prolate and oblate cases.
%
% Arguments:
%     k - the wavenumber
%     a - the interfocal distance divided by two
%     c - k * a
%     m, n - the term to calculate
%     everything - the struct containing the precomputed spheroidal wave
%                  functions
%     eta, xi, phi - the positions of the evaluation points in spheroidal
%                    coordinates
%     eta0, xi0 - the position of the point source in spheroidal coordinates
% Return Values:
%     v - the term
%     dgrad - the gradient of the term in spheroidal coordinates
%     max_abs_change - the maximum absolute change
%
function [dv, dgrad, max_abs_change] = calculate_point_source_in_term(k, a, c, m, n, everything, eta, xi, phi, eta0, xi0)
	if (m > 0)
		epsilon = 2.0;
	else
		epsilon = 1.0;
	end
	xi_less = min(xi0 * ones(1, length(eta)), xi);
	xi_greater = max(xi0 * ones(1, length(eta)), xi);
	S1 = interp1(everything.S1.eta(everything.S1.S1_idxs), everything.S1.S1(everything.S1.S1_idxs), eta0, 'spline');
	A = ((1i * k) / (2.0 * pi)) * epsilon * S1;
	S1 = interp1(everything.S1.eta(everything.S1.S1_idxs), everything.S1.S1(everything.S1.S1_idxs), eta, 'spline');
	S1p = interp1(everything.S1.eta(everything.S1.S1p_idxs), everything.S1.S1p(everything.S1.S1p_idxs), eta, 'spline');
	R1_less = interp1(everything.R.xi(everything.R.idxs), everything.R.R1(everything.R.idxs), xi_less, 'spline');
	R1p_less = interp1(everything.R.xi(everything.R.idxs), everything.R.R1p(everything.R.idxs), xi_less, 'spline');
	R3_greater = interp1(everything.R.xi(everything.R.idxs), everything.R.R3(everything.R.idxs), xi_greater, 'spline');
	R3p_greater = interp1(everything.R.xi(everything.R.idxs), everything.R.R3p(everything.R.idxs), xi_greater, 'spline');
	Phi = cos(m * phi);
	Phip = -m * sin(m * phi);
	dv = A * S1 .* R1_less .* R3_greater .* Phi;
	dgrad = zeros(3, length(eta)); % TODO.
	max_abs_change = max(abs([dv, dgrad(1, :), dgrad(2, :), dgrad(3, :)]));
end
