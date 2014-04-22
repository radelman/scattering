%
% Calculate the scattered field from an incident field due to a point source
% striking a sound-soft prolate spheroid.
%
% Arguments:
%     k - the wavenumber
%     a - the interfocal distance divided by two
%     eta0, xi0 - the position of the point source in prolate spheroidal
%                 coordinates
%     path - the directory in which the prolate spheroidal wave function have
%            been precomputed and stored
%     xi1 - the size of the prolate spheroid
%     x, y, z - the positions of the evaluation points in Cartesian coordinates
% Return Values:
%     v_scat - the scattered field
%     max_abs_change - the maximum absolute change per mode over all the
%                      evaluation points
%
function [v_scat, max_abs_change] = pro_point_source_scat_soft(k, a, eta0, xi0, path, xi1, x, y, z)
	[ ...
	v_scat, max_abs_change ...
	] = pro_calculate_sum(k, a, path, x, y, z, @(k, a, c, m, n, everything, eta, xi, phi)(calculate_point_source_scat_term(k, a, c, m, n, everything, eta, xi, phi, eta0, xi0, xi1, 'soft')));
end
