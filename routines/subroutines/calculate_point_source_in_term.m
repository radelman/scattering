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
%     change - the term
%     max_abs_change - the maximum absolute term over all the evaluation points
%
function [change, max_abs_change] = calculate_point_source_in_term(k, a, c, m, n, everything, eta, xi, phi, eta0, xi0)
	if (m == 0)
		epsilon = 1.0;
	else
		epsilon = 2.0;
	end
	xi_less = min(xi0 * ones(1, length(eta)), xi);
	xi_greater = max(xi0 * ones(1, length(eta)), xi);
	change = ...
	((1i * k) / (2.0 * pi)) * ...
	epsilon * ...
	interp1(everything.S1.eta, everything.S1.S1, eta0, 'spline') * ...
	interp1(everything.S1.eta, everything.S1.S1, eta, 'spline') .* ...
	interp1(everything.R.xi, everything.R.R1, xi_less, 'spline') .* ...
	interp1(everything.R.xi, everything.R.R1 + 1i * everything.R.R2, xi_greater, 'spline') .* ...
	cos(m * phi);
	max_abs_change = max(abs(change));
end
