%
% Calculate a prolate spheroidal wave function expansion.  This function
% doesn't contain any code for specific cases; instead, it calls a
% user-specified function to calculate each term in the expansion.  It also
% handles converting the positions of the evaluation points from Cartesian
% coordinates to prolate spheroidal coordinates and loading the precomputed
% prolate spheroidal wave functions.
%
% Arguments:
%     k - the wavenumber
%     a - the interfocal distance divided by two
%     path - the directory in which the prolate spheroidal wave function have
%            been precomputed and stored
%     x, y, z - the positions of the evaluation points in Cartesian coordinates
%     pro_calculate_term - the function to call to calculate each term, and
%                          should be of the form, pro_calculate_term(k, a, c,
%                          m, n, everything, eta, xi, phi)
% Return Values:
%     v - the expansion
%     max_abs_change - the maximum absolute change per term over all the
%                      evaluation points
%
function [v, max_abs_change] = pro_calculate_sum(k, a, path, x, y, z, pro_calculate_term)
	c = k * a;
	cart = [x; y; z];
	pro = cart_to_pro(a, cart);
	eta = pro(1, :);
	xi = pro(2, :);
	phi = pro(3, :);
	v = zeros(1, length(x));
	max_abs_change = [];
	for m = 0 : 499
		break_again = 0;
		for n = m : m + 499
			try
				everything = pro_open_everything(path, c, m, n);
				good_idxs = isfinite(everything.R.R1) & isfinite(everything.R.R1p) & isfinite(everything.R.R2) & isfinite(everything.R.R2p);
				everything.R.xi = everything.R.xi(good_idxs);
				everything.R.R1 = everything.R.R1(good_idxs);
				everything.R.R1p = everything.R.R1p(good_idxs);
				everything.R.R2 = everything.R.R2(good_idxs);
				everything.R.R2p = everything.R.R2p(good_idxs);
				[ ...
				change, max_abs_change(m + 1, n - m + 1) ...
				] = pro_calculate_term(k, a, c, m, n, everything, eta, xi, phi);
				v = v + change;
				if (max_abs_change(m + 1, n - m + 1) > 0.0)
					fprintf('m = %3d, n = %3d, n - m = %3d, max_abs_change = %.5e\n', m, n, n - m, max_abs_change(m + 1, n - m + 1));
					if (max_abs_change(m + 1, n - m + 1) < eps)
						break_again = break_again + 1;
					else
						break_again = 0;
					end
				else
					break_again = break_again + 1;
				end
			catch exception
				break_again = break_again + 1;
			end
			if (break_again == 3)
				break;
			end
		end
	end
end
