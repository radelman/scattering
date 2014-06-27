%
% Calculate a term in the spheroidal wave function expansion for calculating
% the scattered field from an incident field due to a point source striking a
% spheroid.  This code works for both the prolate and oblate cases.
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
%     xi1 - the size of the spheroid
%     alpha - can be 'soft' for sound-soft spheroid, 'hard' for sound-hard
%             spheroid, or a numerical value.  for the latter, on the boundary
%             of the spheroid, v + alpha * dv/dn = 0
% Return Values:
%     dv - the term
%     dgrad - the gradient of the term in spheroidal coordinates
%     max_abs_change - the maximum absolute change
%
function [dv, dgrad, max_abs_change] = calculate_point_source_scat_term(k, a, c, m, n, everything, eta, xi, phi, eta0, xi0, xi1, alpha)
	if (m > 0)
		epsilon = 2.0;
	else
		epsilon = 1.0;
	end
	S1 = interp1(everything.S1.eta(everything.S1.S1_idxs), everything.S1.S1(everything.S1.S1_idxs), eta0, 'spline');
	R3 = interp1(everything.R.xi(everything.R.idxs), everything.R.R3(everything.R.idxs), xi0, 'spline');
	A = ((1i * k) / (2.0 * pi)) * epsilon * S1 * R3;
	R1 = interp1(everything.R.xi(everything.R.idxs), everything.R.R1(everything.R.idxs), xi1, 'spline');
	R1p = interp1(everything.R.xi(everything.R.idxs), everything.R.R1p(everything.R.idxs), xi1, 'spline');
	R3 = interp1(everything.R.xi(everything.R.idxs), everything.R.R3(everything.R.idxs), xi1, 'spline');
	R3p = interp1(everything.R.xi(everything.R.idxs), everything.R.R3p(everything.R.idxs), xi1, 'spline');
	if (ischar(alpha))
		if (strcmp(alpha, 'soft'))
			B = -A * (R1 / R3);
		else
			B = -A * (R1p / R3p);
		end
	else
		B = -A * ((R1 + alpha * R1p) / (R3 + alpha * R3p));
	end
	S1 = interp1(everything.S1.eta(everything.S1.S1_idxs), everything.S1.S1(everything.S1.S1_idxs), eta, 'spline');
	S1p = interp1(everything.S1.eta(everything.S1.S1p_idxs), everything.S1.S1p(everything.S1.S1p_idxs), eta, 'spline');
	R3 = interp1(everything.R.xi(everything.R.idxs), everything.R.R3(everything.R.idxs), xi, 'spline');
	R3p = interp1(everything.R.xi(everything.R.idxs), everything.R.R3p(everything.R.idxs), xi, 'spline');
	Phi = cos(m * phi);
	Phip = -m * sin(m * phi);
	dv = B * S1 .* R3 .* Phi;
	dgrad = B * [S1p .* R3 .* Phi; S1 .* R3p .* Phi; S1 .* R3 .* Phip];
	max_abs_change = max(abs([dv(xi >= xi1), dgrad(1, xi >= xi1), dgrad(2, xi >= xi1), dgrad(3, xi >= xi1)]));
end
