%
% Calculate the incident field due to a plane wave using prolate spheroidal
% wave functions.
%
% Arguments:
%     k - the wavenumber
%     a - the interfocal distance divided by two
%     theta0 - the direction of the plane wave in spherical coordinates
%     path - the directory in which the prolate spheroidal wave functions have
%            been precomputed and stored
%     x, y, z - the positions of the evaluation points in Cartesian coordinates
% Return Values:
%     v_in - the incident field
%     grad_in_cart - the gradient of the incident field in Cartesian
%                    coordinates
%     max_abs_change_in - the maximum absolute change per term
%
function [v_in, grad_in_cart, max_abs_change_in] = pro_plane_wave_in(k, a, theta0, path, x, y, z)
	[ ...
	v_in, grad_in_cart, max_abs_change_in ...
	] = pro_calculate_sum(k, a, path, x, y, z, @(k, a, c, m, n, everything, eta, xi, phi)(calculate_plane_wave_in_term(k, a, c, m, n, everything, eta, xi, phi, theta0)));
end
