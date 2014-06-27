%
% Calculate the scattered field from a plane wave striking an oblate spheroid
% with Robin boundary conditions.
%
% Arguments:
%     k - the wavenumber
%     a - the interfocal distance divided by two
%     theta0 - the direction of the plane wave in spherical coordinates
%     path - the directory in which the oblate spheroidal wave functions have
%            been precomputed and stored
%     xi1 - the size of the oblate spheroid
%     alpha - on the boundary of the spheroid, v + alpha * dv/dn = 0
%     x, y, z - the positions of the evaluation points in Cartesian coordinates
% Return Values:
%     v_scat - the scattered field
%     grad_scat_cart - the gradient of the scattered field in Cartesian
%                      coordinates
%     max_abs_change_scat - the maximum absolute change per term
%
function [v_scat, grad_scat_cart, max_abs_change_scat] = obl_plane_wave_scat_robin(k, a, theta0, path, xi1, alpha, x, y, z)
	[ ...
	v_scat, grad_scat_cart, max_abs_change_scat ...
	] = obl_calculate_sum(k, a, path, x, y, z, @(k, a, c, m, n, everything, eta, xi, phi)(calculate_plane_wave_scat_term(k, a, c, m, n, everything, eta, xi, phi, theta0, xi1, alpha)));
end
