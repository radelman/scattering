function D = obl_calculate_D(a, obl)
	eta = obl(1);
	xi = obl(2);
	phi = obl(3);
	D = zeros(3, 3);
	D(1, 1) = -a * ((1.0 - (eta ^ 2)) ^ (-1.0 / 2.0)) * eta * sqrt((xi ^ 2) + 1.0) * cos(phi);
	D(1, 2) = -a * ((1.0 - (eta ^ 2)) ^ (-1.0 / 2.0)) * eta * sqrt((xi ^ 2) + 1.0) * sin(phi);
	D(1, 3) = a * xi;
	D(2, 1) = a * sqrt(1.0 - (eta ^ 2)) * (((xi ^ 2) + 1.0) ^ (-1.0 / 2.0)) * xi * cos(phi);
	D(2, 2) = a * sqrt(1.0 - (eta ^ 2)) * (((xi ^ 2) + 1.0) ^ (-1.0 / 2.0)) * xi * sin(phi);
	D(2, 3) = a * eta;
	D(3, 1) = -a * sqrt(1.0 - (eta ^ 2)) * sqrt((xi ^ 2) + 1.0) * sin(phi);
	D(3, 2) = a * sqrt(1.0 - (eta ^ 2)) * sqrt((xi ^ 2) + 1.0) * cos(phi);
	D(3, 3) = 0.0;
end
