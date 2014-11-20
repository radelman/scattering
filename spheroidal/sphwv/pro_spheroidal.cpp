//
// Copyright (c) 2014, Ross Adelman, Nail A. Gumerov, and Ramani Duraiswami
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//

#include "adder.hpp"
#include "common_spheroidal.hpp"
#include <iostream>
#include "pro_spheroidal.hpp"
#include "real.hpp"
#include <vector>

static void calculate_Q(std::vector<real> & Q, bool verbose, const real & m0, const real & n1, const real & xi);

real calculate_c_squared(const real & c)
{
	return c * c;
}

void calculate_Rmn1_1(real & R1, real & R1p, bool verbose, const real & c, const real & m, const real & n, const real & n_dr, const std::vector<real> & dr, const real & F, const real & xi)
{
	calculate_Rmn1_1_shared(R1, R1p, verbose, c, m, n, n_dr, dr, xi);
	if (m > real::ZERO)
	{
		R1p = pow(F, -real::ONE) * (m / real::TWO) * pow(real::ONE - real::ONE / (xi * xi), m / real::TWO - real::ONE) * (real::TWO / (xi * xi * xi)) * R1 + pow(F, -real::ONE) * pow(real::ONE - real::ONE / (xi * xi), m / real::TWO) * R1p;
	}
	else
	{
		R1p = pow(F, -real::ONE) * R1p;
	}
	R1 = pow(F, -real::ONE) * pow(real::ONE - real::ONE / (xi * xi), m / real::TWO) * R1;
}

void calculate_Rmn1_2(real & R1, real & R1p, bool verbose, const real & c, const real & m, const real & n, const real & k1, const real & n_c2k, const std::vector<real> & c2k, const real & xi)
{
	adder R1_adder;
	real a;
	real change;
	adder R1p_adder;
	real ap;
	real changep;
	
	R1 = real::ZERO;
	R1_adder.clear();
	for (real k = real::ZERO; k <= n_c2k - real::ONE; k = k + real::ONE)
	{
		if (k > real::ZERO)
		{
			a = a * -real::ONE * (xi * xi - real::ONE);
		}
		else
		{
			a = real::ONE;
		}
		change = a * c2k[gzbi(k)];
		R1 = R1 + change;
		R1_adder.add(change);
		if (k > real::ZERO && abs(change) > real::ZERO && abs(change / R1) < real::SMALL_ENOUGH)
		{
			if (verbose)
			{
				std::cout << "calculate_Rmn1_2: " << (change / R1).get_string(10) << ", " << ((R1 - R1_adder.calculate_sum()) / R1_adder.calculate_sum()).get_string(10) << ", " << (change / R1_adder.calculate_sum()).get_string(10) << std::endl;
			}
			break;
		}
	}
	R1 = R1_adder.calculate_sum();
	R1p = real::ZERO;
	R1p_adder.clear();
	for (real k = real::ONE; k <= n_c2k - real::ONE; k = k + real::ONE)
	{
		if (k > real::ONE)
		{
			ap = ap * -real::ONE * (xi * xi - real::ONE);
		}
		else
		{
			ap = -real::ONE;
		}
		changep = ap * c2k[gzbi(k)] * k * real::TWO * xi;
		R1p = R1p + changep;
		R1p_adder.add(changep);
		if (k > real::ONE && abs(changep) > real::ZERO && abs(changep / R1p) < real::SMALL_ENOUGH)
		{
			if (verbose)
			{
				std::cout << "calculate_Rmn1_2: " << (changep / R1p).get_string(10) << ", " << ((R1p - R1p_adder.calculate_sum()) / R1p_adder.calculate_sum()).get_string(10) << ", " << (changep / R1p_adder.calculate_sum()).get_string(10) << std::endl;
			}
			break;
		}
	}
	R1p = R1p_adder.calculate_sum();
	if (remainder(n - m, real::TWO) == real::ZERO)
	{
		if (m > real::ZERO)
		{
			R1p = pow(k1, -real::ONE) * (m / real::TWO) * pow(xi * xi - real::ONE, m / real::TWO - real::ONE) * real::TWO * xi * R1 + pow(k1, -real::ONE) * pow(xi * xi - real::ONE, m / real::TWO) * R1p;
		}
		else
		{
			R1p = pow(k1, -real::ONE) * R1p;
		}
		R1 = pow(k1, -real::ONE) * pow(xi * xi - real::ONE, m / real::TWO) * R1;
	}
	else
	{
		if (m > real::ZERO)
		{
			R1p = pow(k1, -real::ONE) * pow(xi * xi - real::ONE, m / real::TWO) * R1 + pow(k1, -real::ONE) * xi * (m / real::TWO) * pow(xi * xi - real::ONE, m / real::TWO - real::ONE) * real::TWO * xi * R1 + pow(k1, -real::ONE) * xi * pow(xi * xi - real::ONE, m / real::TWO) * R1p;
		}
		else
		{
			R1p = pow(k1, -real::ONE) * R1 + pow(k1, -real::ONE) * xi * R1p;
		}
		R1 = pow(k1, -real::ONE) * xi * pow(xi * xi - real::ONE, m / real::TWO) * R1;
	}
}

void calculate_Rmn2_1(real & R2, real & R2p, bool verbose, const real & c, const real & m, const real & n, const real & n_dr, const std::vector<real> & dr, const real & F, const real & xi)
{
	calculate_Rmn2_1_shared(R2, R2p, verbose, c, m, n, n_dr, dr, xi);
	if (m > real::ZERO)
	{
		R2p = pow(F, -real::ONE) * (m / real::TWO) * pow(real::ONE - real::ONE / (xi * xi), m / real::TWO - real::ONE) * (real::TWO / (xi * xi * xi)) * R2 + pow(F, -real::ONE) * pow(real::ONE - real::ONE / (xi * xi), m / real::TWO) * R2p;
	}
	else
	{
		R2p = pow(F, -real::ONE) * R2p;
	}
	R2 = pow(F, -real::ONE) * pow(real::ONE - real::ONE / (xi * xi), m / real::TWO) * R2;
}

static void calculate_Q(std::vector<real> & Q, bool verbose, const real & m0, const real & n1, const real & xi)
{
	real Q0;
	real Q1;
	real Q2;
	real b0;
	std::vector<real> a;
	std::vector<real> b;
	real prev_N;
	real N;
	
	Q.clear();
	for (real n = real::ZERO; n <= n1; n = n + real::ONE)
	{
		Q.push_back(real::ZERO);
	}
	Q0 = (real::ONE / real::TWO) * log((xi + real::ONE) / (xi - real::ONE));
	Q1 = -pow(xi * xi - real::ONE, -real::ONE / real::TWO);
	if (m0 > real::ONE)
	{
		for (real m = real::TWO; m <= m0; m = m + real::ONE)
		{
			Q2 = (m - real::ONE) * (-m + real::TWO) * Q0 - ((real::TWO * (m - real::ONE)) / pow(xi * xi - real::ONE, real::ONE / real::TWO)) * xi * Q1;
			Q0 = Q1;
			Q1 = Q2;
		}
		Q[gzbi(real::ZERO)] = Q2;
	}
	else if (m0 > real::ZERO)
	{
		Q[gzbi(real::ZERO)] = Q1;
	}
	else
	{
		Q[gzbi(real::ZERO)] = Q0;
	}
	for (real n = n1; n >= real::ONE; n = n - real::ONE)
	{
		if (n < n1)
		{
			N = (n + m0) / ((real::TWO * n + real::ONE) * xi - (n - m0 + real::ONE) * N);
		}
		else
		{
			b0 = real::ZERO;
			a.clear();
			a.push_back(real::ONE);
			b.clear();
			prev_N = real::NAN;
			for (real i = n; i <= n + real("8388608.0"); i = i + real::ONE)
			{
				if (remainder(i - n, real("10000.0")) == real::ZERO)
				{
					if (verbose)
					{
						std::cout << "calculate_Q: " << (i - n).get_int() << std::endl;
					}
				}
				a.back() = a.back() * (i + m0);
				b.push_back((real::TWO * i + real::ONE) * xi);
				a.push_back(i - m0 + real::ONE);
				if (i > n && pow(real::TWO, round(log(i - n) / log(real::TWO))) == i - n)
				{
					N = -calculate_continued_fraction(b0, a, b);
					if (prev_N == prev_N)
					{
						if (verbose)
						{
							std::cout << "calculate_Q: " << (i - n).get_int() << ": " << ((N - prev_N) / prev_N).get_string(10) << std::endl;
						}
						if (abs((N - prev_N) / prev_N) < real::SMALL_ENOUGH)
						{
							break;
						}
					}
					prev_N = N;
				}
			}
		}
		Q[gzbi(n)] = N;
	}
	for (real n = real::ONE; n <= n1; n = n + real::ONE)
	{
		Q[gzbi(n)] = Q[gzbi(n)] * Q[gzbi(n - real::ONE)];
	}
}

void calculate_Rmn2_2(real & R2, real & R2p, bool verbose, const real & c, const real & m, const real & n, const real & n_dr, const std::vector<real> & dr, const real & n_dr_neg, const std::vector<real> & dr_neg, const real & k2, const real & xi)
{
	real v_max;
	std::vector<real> P;
	std::vector<real> Pp;
	std::vector<real> Q;
	std::vector<real> Qp;
	real Q1;
	real Q2;
	adder R2_adder;
	adder R2p_adder;
	real change;
	real changep;
	
	if (xi == real::ONE)
	{
		R2 = real::NAN;
		R2p = real::NAN;
		return;
	}
	// The following code assumes there are at least m + 2 entries in P, so
	// make sure there are.
	v_max = max(m + real::ONE, n_dr_neg - m);
	P.clear();
	for (real v = real::ZERO; v <= v_max; v = v + real::ONE)
	{
		P.push_back(real::ZERO);
	}
	P[gzbi(real::ZERO)] = real::ONE;
	for (real v = real::ONE; v <= m; v = v + real::ONE)
	{
		P[gzbi(v)] = (real::TWO * v - real::ONE) * pow(xi * xi - real::ONE, real::ONE / real::TWO) * P[gzbi(v - real::ONE)];
		P[gzbi(v - real::ONE)] = real::ZERO;
	}
	P[gzbi(m + real::ONE)] = (real::TWO * m + real::ONE) * xi * P[gzbi(m)];
	for (real v = m + real::TWO; v <= v_max; v = v + real::ONE)
	{
		P[gzbi(v)] = (real::ONE / (v - m)) * (-(v + m - real::ONE) * P[gzbi(v - real::TWO)] + (real::TWO * v - real::ONE) * xi * P[gzbi(v - real::ONE)]);
	}
	Pp.clear();
	for (real v = real::ZERO; v <= v_max - real::ONE; v = v + real::ONE)
	{
		Pp.push_back(real::ZERO);
	}
	for (real v = m; v <= v_max - real::ONE; v = v + real::ONE)
	{
		Pp[gzbi(v)] = (real::ONE / (xi * xi - real::ONE)) * (-(v + real::ONE) * xi * P[gzbi(v)] + (v - m + real::ONE) * P[gzbi(v + real::ONE)]);
	}
	calculate_Q(Q, verbose, m, m + n_dr, xi);
	Qp.clear();
	for (real v = real::ZERO; v <= m + n_dr - real::ONE; v = v + real::ONE)
	{
		Qp.push_back(real::ZERO);
	}
	for (real v = real::ZERO; v <= m + n_dr - real::ONE; v = v + real::ONE)
	{
		Qp[gzbi(v)] = (real::ONE / (xi * xi - real::ONE)) * (-(v + 1) * xi * Q[gzbi(v)] + (v - m + real::ONE) * Q[gzbi(v + real::ONE)]);
	}
	Q1 = Q[gzbi(real::ZERO)];
	Q2 = Q[gzbi(real::ONE)];
	for (real v = -real::ONE; v >= -m; v = v - real::ONE)
	{
		P[gnobi(v)] = (real::ONE / (v + m + real::ONE)) * ((real::TWO * v + real::THREE) * xi * Q1 - (v - m + real::TWO) * Q2);
		Pp[gnobi(v)] = (real::ONE / (xi * xi - real::ONE)) * (-(v + real::ONE) * xi * P[gnobi(v)] + (v - m + real::ONE) * Q1);
		Q2 = Q1;
		Q1 = P[gnobi(v)];
	}
	R2 = real::ZERO;
	R2_adder.clear();
	R2p = real::ZERO;
	R2p_adder.clear();
	if (remainder(n - m, real::TWO) == real::ZERO)
	{
		for (real r = -real::TWO * m - real::TWO; r >= -n_dr_neg; r = r - real::TWO)
		{
			change = dr_neg[gnobi(r)] * P[gzbi(-r - m - real::ONE)];
			R2 = R2 + change;
			R2_adder.add(change);
			changep = dr_neg[gnobi(r)] * Pp[gzbi(-r - m - real::ONE)];
			R2p = R2p + changep;
			R2p_adder.add(changep);
			if (r < -real::TWO * m - real::TWO && abs(change) > real::ZERO && abs(change / R2) < real::SMALL_ENOUGH && abs(changep) > real::ZERO && abs(changep / R2p) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Rmn2_2: " << (change / R2).get_string(10) << ", " << (changep / R2p).get_string(10) << ", " << ((R2 - R2_adder.calculate_sum()) / R2_adder.calculate_sum()).get_string(10) << ", " << ((R2p - R2p_adder.calculate_sum()) / R2p_adder.calculate_sum()).get_string(10) << ", " << (change / R2_adder.calculate_sum()).get_string(10) << ", " << (changep / R2p_adder.calculate_sum()).get_string(10) << std::endl;
				}
				break;
			}
		}
		// There is a max here because n_dr_neg may be smaller than m.
		for (real r = -real::TWO; r >= max(-n_dr_neg, -real::TWO * m); r = r - real::TWO)
		{
			if (m + r >= real::ZERO)
			{
				change = dr_neg[gnobi(r)] * Q[gzbi(m + r)];
				changep = dr_neg[gnobi(r)] * Qp[gzbi(m + r)];
			}
			else
			{
				change = dr_neg[gnobi(r)] * P[gnobi(m + r)];
				changep = dr_neg[gnobi(r)] * Pp[gnobi(m + r)];
			}
			R2 = R2 + change;
			R2_adder.add(change);
			R2p = R2p + changep;
			R2p_adder.add(changep);
			if (r < -real::TWO && abs(change) > real::ZERO && abs(change / R2) < real::SMALL_ENOUGH && abs(changep) > real::ZERO && abs(changep / R2p) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Rmn2_2: " << (change / R2).get_string(10) << ", " << (changep / R2p).get_string(10) << ", " << ((R2 - R2_adder.calculate_sum()) / R2_adder.calculate_sum()).get_string(10) << ", " << ((R2p - R2p_adder.calculate_sum()) / R2p_adder.calculate_sum()).get_string(10) << ", " << (change / R2_adder.calculate_sum()).get_string(10) << ", " << (changep / R2p_adder.calculate_sum()).get_string(10) << std::endl;
				}
				break;
			}
		}
		for (real r = real::ZERO; r <= n_dr - real::TWO; r = r + real::TWO)
		{
			change = dr[gzbi(r)] * Q[gzbi(m + r)];
			R2 = R2 + change;
			R2_adder.add(change);
			changep = dr[gzbi(r)] * Qp[gzbi(m + r)];
			R2p = R2p + changep;
			R2p_adder.add(changep);
			if (r > real::ZERO && abs(change) > real::ZERO && abs(change / R2) < real::SMALL_ENOUGH && abs(changep) > real::ZERO && abs(changep / R2p) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Rmn2_2: " << (change / R2).get_string(10) << ", " << (changep / R2p).get_string(10) << ", " << ((R2 - R2_adder.calculate_sum()) / R2_adder.calculate_sum()).get_string(10) << ", " << ((R2p - R2p_adder.calculate_sum()) / R2p_adder.calculate_sum()).get_string(10) << ", " << (change / R2_adder.calculate_sum()).get_string(10) << ", " << (changep / R2p_adder.calculate_sum()).get_string(10) << std::endl;
				}
				break;
			}
		}
	}
	else
	{
		for (real r = -real::TWO * m - real::ONE; r >= -n_dr_neg + real::ONE; r = r - real::TWO)
		{
			change = dr_neg[gnobi(r)] * P[gzbi(-r - m - real::ONE)];
			R2 = R2 + change;
			R2_adder.add(change);
			changep = dr_neg[gnobi(r)] * Pp[gzbi(-r - m - real::ONE)];
			R2p = R2p + changep;
			R2p_adder.add(changep);
			if (r < -real::TWO * m - real::ONE && abs(change) > real::ZERO && abs(change / R2) < real::SMALL_ENOUGH && abs(changep) > real::ZERO && abs(changep / R2p) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Rmn2_2: " << (change / R2).get_string(10) << ", " << (changep / R2p).get_string(10) << ", " << ((R2 - R2_adder.calculate_sum()) / R2_adder.calculate_sum()).get_string(10) << ", " << ((R2p - R2p_adder.calculate_sum()) / R2p_adder.calculate_sum()).get_string(10) << ", " << (change / R2_adder.calculate_sum()).get_string(10) << ", " << (changep / R2p_adder.calculate_sum()).get_string(10) << std::endl;
				}
				break;
			}
		}
		// See the corresponding comment when n - m = even.
		for (real r = -real::ONE; r >= max(-n_dr_neg + real::ONE, -real::TWO * m + real::ONE); r = r - real::TWO)
		{
			if (m + r >= real::ZERO)
			{
				change = dr_neg[gnobi(r)] * Q[gzbi(m + r)];
				changep = dr_neg[gnobi(r)] * Qp[gzbi(m + r)];
			}
			else
			{
				change = dr_neg[gnobi(r)] * P[gnobi(m + r)];
				changep = dr_neg[gnobi(r)] * Pp[gnobi(m + r)];
			}
			R2 = R2 + change;
			R2_adder.add(change);
			R2p = R2p + changep;
			R2p_adder.add(changep);
			if (r < -real::ONE && abs(change) > real::ZERO && abs(change / R2) < real::SMALL_ENOUGH && abs(changep) > real::ZERO && abs(changep / R2p) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Rmn2_2: " << (change / R2).get_string(10) << ", " << (changep / R2p).get_string(10) << ", " << ((R2 - R2_adder.calculate_sum()) / R2_adder.calculate_sum()).get_string(10) << ", " << ((R2p - R2p_adder.calculate_sum()) / R2p_adder.calculate_sum()).get_string(10) << ", " << (change / R2_adder.calculate_sum()).get_string(10) << ", " << (changep / R2p_adder.calculate_sum()).get_string(10) << std::endl;
				}
				break;
			}
		}
		for (real r = real::ONE; r <= n_dr - real::ONE; r = r + real::TWO)
		{
			change = dr[gzbi(r)] * Q[gzbi(m + r)];
			R2 = R2 + change;
			R2_adder.add(change);
			changep = dr[gzbi(r)] * Qp[gzbi(m + r)];
			R2p = R2p + changep;
			R2p_adder.add(changep);
			if (r > real::ONE && abs(change) > real::ZERO && abs(change / R2) < real::SMALL_ENOUGH && abs(changep) > real::ZERO && abs(changep / R2p) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Rmn2_2: " << (change / R2).get_string(10) << ", " << (changep / R2p).get_string(10) << ", " << ((R2 - R2_adder.calculate_sum()) / R2_adder.calculate_sum()).get_string(10) << ", " << ((R2p - R2p_adder.calculate_sum()) / R2p_adder.calculate_sum()).get_string(10) << ", " << (change / R2_adder.calculate_sum()).get_string(10) << ", " << (changep / R2p_adder.calculate_sum()).get_string(10) << std::endl;
				}
				break;
			}
		}
	}
	R2 = R2_adder.calculate_sum();
	R2p = R2p_adder.calculate_sum();
	R2 = pow(k2, -real::ONE) * R2;
	R2p = pow(k2, -real::ONE) * R2p;
}
