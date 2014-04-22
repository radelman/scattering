%
% Calculate the incident field due to a point source using oblate spheroidal
% wave functions.
%
% Arguments:
%     k - the wavenumber
%     a - the interfocal distance divided by two
%     eta0, xi0 - the position of the point source in oblate spheroidal
%                 coordinates
%     path - the directory in which the oblate spheroidal wave function have
%            been precomputed and stored
%     x, y, z - the positions of the evaluation points in Cartesian coordinates
% Return Values:
%     v_in - the incident field
%     max_abs_change - the maximum absolute change per mode over all the
%                      evaluation points
%
function [v_in, max_abs_change] = obl_point_source_in(k, a, eta0, xi0, path, x, y, z)
	[ ...
	v_in, max_abs_change ...
	] = obl_calculate_sum(k, a, path, x, y, z, @(k, a, c, m, n, everything, eta, xi, phi)(calculate_point_source_in_term(k, a, c, m, n, everything, eta, xi, phi, eta0, xi0)));
end
