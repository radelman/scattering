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
