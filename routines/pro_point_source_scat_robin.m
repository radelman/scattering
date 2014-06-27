%
% Calculate the scattered field from an incident field due to a point source
% striking a prolate spheroid with Robin boundary conditions.
%
% Arguments:
%     k - the wavenumber
%     a - the interfocal distance divided by two
%     eta0, xi0 - the position of the point source in prolate spheroidal
%                 coordinates
%     path - the directory in which the prolate spheroidal wave functions have
%            been precomputed and stored
%     xi1 - the size of the prolate spheroid
%     alpha - on the boundary of the spheroid, v + alpha * dv/dn = 0
%     x, y, z - the positions of the evaluation points in Cartesian coordinates
% Return Values:
%     v_scat - the scattered field
%     grad_scat_cart - the gradient of the scattered field in Cartesian
%                      coordinates
%     max_abs_change_scat - the maximum absolute change per term
%
function [v_scat, grad_scat_cart, max_abs_change_scat] = pro_point_source_scat_robin(k, a, eta0, xi0, path, xi1, alpha, x, y, z)
	[ ...
	v_scat, grad_scat_cart, max_abs_change_scat ...
	] = pro_calculate_sum(k, a, path, x, y, z, @(k, a, c, m, n, everything, eta, xi, phi)(calculate_point_source_scat_term(k, a, c, m, n, everything, eta, xi, phi, eta0, xi0, xi1, alpha)));
end
