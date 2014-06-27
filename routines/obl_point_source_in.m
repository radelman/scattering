%
% Calculate the incident field due to a point source using oblate spheroidal
% wave functions.
%
% Arguments:
%     k - the wavenumber
%     a - the interfocal distance divided by two
%     eta0, xi0 - the position of the point source in oblate spheroidal
%                 coordinates
%     path - the directory in which the oblate spheroidal wave functions have
%            been precomputed and stored
%     x, y, z - the positions of the evaluation points in Cartesian coordinates
% Return Values:
%     v_in - the incident field
%     grad_in_cart - the gradient of the incident field in Cartesian
%                    coordinates
%     max_abs_change_in - the maximum absolute change per term
%
function [v_in, grad_in_cart, max_abs_change_in] = obl_point_source_in(k, a, eta0, xi0, path, x, y, z)
	[ ...
	v_in, grad_in_cart, max_abs_change_in ...
	] = obl_calculate_sum(k, a, path, x, y, z, @(k, a, c, m, n, everything, eta, xi, phi)(calculate_point_source_in_term(k, a, c, m, n, everything, eta, xi, phi, eta0, xi0)));
end
