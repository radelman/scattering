%
% Calculate the incident field due to a point source.
%
% Arguments:
%     k - the wavenumber
%     px, py, pz - the position of the point source in Cartesian coordinates
%     x, y, z - the positions of the evaluation points in Cartesian coordinates
% Return Values:
%     v_in - the incident field
%
function v_in = point_source_in(k, px, py, pz, x, y, z)
	r = sqrt((x - px) .^ 2 + (y - py) .^ 2 + (z - pz) .^ 2);
	v_in = exp(1i * k * r) ./ (4.0 * pi * r);
end
