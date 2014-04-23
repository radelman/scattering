%
% Calculate a term in the spheroidal wave function expansion for calculating
% the incident field due to a plane wave.  This code works for both the prolate
% and oblate cases.
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
%     theta0 - the direction of the plane wave in spherical coordinates
% Return Values:
%     change - the term
%     max_abs_change - the maximum absolute term over all the evaluation points
%
function [change, max_abs_change] = calculate_plane_wave_in_term(k, a, c, m, n, everything, eta, xi, phi, theta0)
	if (m == 0)
		epsilon = 1.0;
	else
		epsilon = 2.0;
	end
	change = ...
	2.0 * ...
	epsilon * ...
	(1i ^ n) * ...
	interp1(everything.S1.eta, everything.S1.S1, cos(theta0), 'spline') * ...
	interp1(everything.S1.eta, everything.S1.S1, eta, 'spline') .* ...
	interp1(everything.R.xi, everything.R.R1, xi, 'spline') .* ...
	cos(m * phi);
	max_abs_change = max(abs(change));
end
