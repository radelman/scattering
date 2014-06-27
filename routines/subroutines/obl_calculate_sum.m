%
% Calculate an oblate spheroidal wave function expansion.  This function
% doesn't contain any code for specific cases; instead, it calls a
% user-specified function to calculate each term in the expansion.  It also
% handles converting the positions of the evaluation points from Cartesian
% coordinates to oblate spheroidal coordinates and loading the precomputed
% oblate spheroidal wave functions.
%
% Arguments:
%     k - the wavenumber
%     a - the interfocal distance divided by two
%     path - the directory in which the oblate spheroidal wave functions have
%            been precomputed and stored
%     x, y, z - the positions of the evaluation points in Cartesian coordinates
%     obl_calculate_term - the function to call to calculate each term, and
%                          should be of the form, [dv, dgrad_obl,
%                          max_abs_change] = obl_calculate_term(k, a, c, m, n,
%                          everything, eta, xi, phi)
% Return Values:
%     v - the expansion
%     grad_cart - the gradient of the expansion in Cartesian coordinates
%     max_abs_change - the maximum absolute change per term
%
function [v, grad_cart, max_abs_change] = obl_calculate_sum(k, a, path, x, y, z, obl_calculate_term)
	c = k * a;
	cart = [x; y; z];
	obl = cart_to_obl(a, cart);
	eta = obl(1, :);
	xi = obl(2, :);
	phi = obl(3, :);
	v = zeros(1, length(x));
	grad_obl = zeros(3, length(x));
	max_abs_change = [];
	for m = 0 : 499
		break_again = 0;
		for n = m : m + 499
			try
				everything = obl_open_everything(path, c, m, n);
				everything.R.R3 = everything.R.R1 + 1i * everything.R.R2;
				everything.R.R3p = everything.R.R1p + 1i * everything.R.R2p;
				everything.S1.S1_idxs = isfinite(everything.S1.S1);
				everything.S1.S1p_idxs = isfinite(everything.S1.S1p);
				everything.R.idxs = isfinite(everything.R.R1) & ...
				                    isfinite(everything.R.R1p) & ...
				                    isfinite(everything.R.R2) & ...
				                    isfinite(everything.R.R2p);
				[ ...
				dv, dgrad_obl, max_abs_change(m + 1, n - m + 1) ...
				] = obl_calculate_term(k, a, c, m, n, everything, eta, xi, phi);
				v = v + dv;
				grad_obl = grad_obl + dgrad_obl;
				fprintf('m = %3d, n = %3d, max_abs_change = %.5e\n', m, n, max_abs_change(m + 1, n - m + 1));
				if (max_abs_change(m + 1, n - m + 1) < eps)
					break_again = break_again + 1;
				else
					break_again = 0;
				end
			catch exception
				break_again = break_again + 1;
			end
			if (break_again == 3)
				break;
			end
		end
	end
	grad_cart = grad_obl_to_cart(a, obl, grad_obl);
end
