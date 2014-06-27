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
%     v - the term
%     dgrad - the gradient of the term in spheroidal coordinates
%     max_abs_change - the maximum absolute change
%
function [dv, dgrad, max_abs_change] = calculate_plane_wave_in_term(k, a, c, m, n, everything, eta, xi, phi, theta0)
	if (m > 0)
		epsilon = 2.0;
	else
		epsilon = 1.0;
	end
	S1 = interp1(everything.S1.eta(everything.S1.S1_idxs), everything.S1.S1(everything.S1.S1_idxs), cos(theta0), 'spline');
	A = 2.0 * epsilon * (1i ^ n) * S1;
	S1 = interp1(everything.S1.eta(everything.S1.S1_idxs), everything.S1.S1(everything.S1.S1_idxs), eta, 'spline');
	S1p = interp1(everything.S1.eta(everything.S1.S1p_idxs), everything.S1.S1p(everything.S1.S1p_idxs), eta, 'spline');
	R1 = interp1(everything.R.xi(everything.R.idxs), everything.R.R1(everything.R.idxs), xi, 'spline');
	R1p = interp1(everything.R.xi(everything.R.idxs), everything.R.R1p(everything.R.idxs), xi, 'spline');
	Phi = cos(m * phi);
	Phip = -m * sin(m * phi);
	dv = A * S1 .* R1 .* Phi;
	dgrad = A * [S1p .* R1 .* Phi; S1 .* R1p .* Phi; S1 .* R1 .* Phip];
	max_abs_change = max(abs([dv, dgrad(1, :), dgrad(2, :), dgrad(3, :)]));
end
