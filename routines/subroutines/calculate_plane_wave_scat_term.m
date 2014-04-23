%
% Calculate a term in the spheroidal wave function expansion for calculating
% the the scattered field from a plane wave striking a spheroid.  This code
% works for both the prolate and oblate cases.
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
%     xi1 - the size of the spheroid
%     imped - whether the spheroid is sound soft or hard: use 'soft' for sound
%             soft and 'hard' or sound hard
% Return Values:
%     change - the term
%     max_abs_change - the maximum absolute term over all the evaluation points
%
function [change, max_abs_change] = calculate_plane_wave_scat_term(k, a, c, m, n, everything, eta, xi, phi, theta0, xi1, imped)
	if (m == 0)
		epsilon = 1.0;
	else
		epsilon = 2.0;
	end
	if (strcmp(imped, 'soft'))
		temp = interp1(everything.R.xi, everything.R.R1, xi1, 'spline') / interp1(everything.R.xi, everything.R.R1 + 1i * everything.R.R2, xi1, 'spline');
	else
		temp = interp1(everything.R.xi, everything.R.R1p, xi1, 'spline') / interp1(everything.R.xi, everything.R.R1p + 1i * everything.R.R2p, xi1, 'spline');
	end
	change = ...
	-2.0 * ...
	epsilon * ...
	(1i ^ n) * ...
	interp1(everything.S1.eta, everything.S1.S1, cos(theta0), 'spline') * ...
	temp * ...
	interp1(everything.S1.eta, everything.S1.S1, eta, 'spline') .* ...
	interp1(everything.R.xi, everything.R.R1 + 1i * everything.R.R2, xi, 'spline') .* ...
	cos(m * phi);
	max_abs_change = max(abs(change(xi >= xi1)));
end
