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
#include <algorithm>
#include "common_spheroidal.hpp"
#include <iostream>
#include "obl_spheroidal.hpp"
#include "real.hpp"
#include <vector>

static real get_c2k(bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr, real & n_c2k, std::vector<real> & c2k, const real & k);
static void calculate_B2rmn_coefficients(std::vector<real> & alpha, std::vector<real> & beta, std::vector<real> & gamma, std::vector<real> & h, bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr, const real & k1, real & n_c2k, std::vector<real> & c2k, const real & Q, std::vector<real> & h_saved, const real & r0, const real & r1);
static void calculate_B2rmn_forward(std::vector<real> & alpha, std::vector<real> & beta, std::vector<real> & gamma, std::vector<real> & h, std::vector<real> & B2r, bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr, const real & k1, real & n_c2k, std::vector<real> & c2k, const real & Q, std::vector<real> & h_saved, const real & n_B2r);
static void calculate_B2rmn_backward(std::vector<real> & B2r, bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr, const real & k1, real & n_c2k, std::vector<real> & c2k, const real & Q, std::vector<real> & h_saved, const real & n_B2r, std::vector<real> & alpha, std::vector<real> & beta, std::vector<real> & gamma, std::vector<real> & h, const real & r0, const real & B0);
static void calculate_B2rmn_once(std::vector<real> & B2r, bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr, const real & k1, real & n_c2k, std::vector<real> & c2k, const real & Q, std::vector<real> & h_saved, const real & n_B2r);
static complex calculate_continued_fraction(const complex & b0, const std::vector<complex> & a, const std::vector<complex> & b);
static void calculate_Q(std::vector<complex> & Q, bool verbose, const real & m0, const real & n1, const real & xi);

real calculate_c_squared(const real & c)
{
	return -c * c;
}

static real get_c2k(bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr, real & n_c2k, std::vector<real> & c2k, const real & k)
{
	if (k >= n_c2k)
	{
		while (k >= n_c2k)
		{
			n_c2k = real::TWO * n_c2k;
		}
		calculate_c2kmn(c2k, verbose, c, m, n, lambda, n_dr, dr, n_c2k, real::ZERO);
	}
	return c2k[gzbi(k)];
}

real calculate_Qmn(bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr, const real & k1, real & n_c2k, std::vector<real> & c2k)
{
	std::vector<real> Ck;
	std::vector<real> Bk;
	std::vector<real> Ak;
	real Q;
	
	Ck.clear();
	for (real k = real::ZERO; k <= m; k = k + real::ONE)
	{
		Ck.push_back(real::ZERO);
	}
	for (real k = real::ZERO; k <= m; k = k + real::ONE)
	{
		Ck[gzbi(k)] = get_c2k(verbose, c, m, n, lambda, n_dr, dr, n_c2k, c2k, k);
	}
	Bk.clear();
	for (real k = real::ZERO; k <= m; k = k + real::ONE)
	{
		Bk.push_back(real::ZERO);
	}
	for (real i = real::ZERO; i <= m; i = i + real::ONE)
	{
		Bk[gzbi(i)] = real::ZERO;
		for (real k = real::ZERO; k <= i; k = k + real::ONE)
		{
			Bk[gzbi(i)] = Bk[gzbi(i)] + Ck[gzbi(k)] * Ck[gzbi(i - k)];
		}
	}
	Ak.clear();
	for (real k = real::ZERO; k <= m; k = k + real::ONE)
	{
		Ak.push_back(real::ZERO);
	}
	Ak[gzbi(real::ZERO)] = real::ONE / Bk[gzbi(real::ZERO)];
	for (real i = real::ONE; i <= m; i = i + real::ONE)
	{
		Ak[gzbi(i)] = real::ZERO;
		for (real k = real::ZERO; k <= i - real::ONE; k = k + real::ONE)
		{
			Ak[gzbi(i)] = Ak[gzbi(i)] + Ak[gzbi(k)] * Bk[gzbi(i - k)];
		}
		Ak[gzbi(i)] = -(real::ONE / Bk[0]) * Ak[gzbi(i)];
	}
	if (remainder(n - m, real::TWO) == real::ZERO)
	{
		Q = real::ZERO;
		for (real k = real::ZERO; k <= m; k = k + real::ONE)
		{
			Q = Q + Ak[gzbi(k)] * (factorial(real::TWO * m - real::TWO * k) / (pow(real::TWO, m - k) * factorial(m - k) * pow(real::TWO, m - k) * factorial(m - k)));
		}
		Q = k1 * k1 * (real::ONE / c) * Q;
	}
	else
	{
		Q = real::ZERO;
		for (real k = real::ZERO; k <= m; k = k + real::ONE)
		{
			Q = Q + Ak[gzbi(k)] * (factorial(real::TWO * m - real::TWO * k + real::ONE) / (pow(real::TWO, m - k) * factorial(m - k) * pow(real::TWO, m - k) * factorial(m - k)));
		}
		Q = -k1 * k1 * (real::ONE / c) * Q;
	}
	return Q;
}

static void calculate_B2rmn_coefficients(std::vector<real> & alpha, std::vector<real> & beta, std::vector<real> & gamma, std::vector<real> & h, bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr, const real & k1, real & n_c2k, std::vector<real> & c2k, const real & Q, std::vector<real> & h_saved, const real & r0, const real & r1)
{
	real h0;
	adder h0_adder;
	real k0;
	real a;
	real change;
	real h1;
	adder h1_adder;
	
	while (r1 >= real((int)h_saved.size()))
	{
		h_saved.push_back(real::NAN);
	}
	alpha.clear();
	beta.clear();
	gamma.clear();
	h.clear();
	for (real r = r0; r <= r1; r = r + real::ONE)
	{
		alpha.push_back(real::ZERO);
		beta.push_back(real::ZERO);
		gamma.push_back(real::ZERO);
		h.push_back(real::ZERO);
	}
	if (remainder(n - m, real::TWO) == real::ZERO)
	{
		for (real r = r0; r <= r1; r = r + real::ONE)
		{
			alpha[gzbi(r - r0)] = (real::TWO * r + real::TWO) * (real::TWO * r + real::THREE);
			beta[gzbi(r - r0)] = (real::TWO * r + real::ONE) * (real::TWO * r - real::TWO * m + real::TWO) + m * (m - real::ONE) - lambda;
			gamma[gzbi(r - r0)] = c * c;
			if (h_saved[gzbi(r)] == h_saved[gzbi(r)])
			{
				h[gzbi(r - r0)] = h_saved[gzbi(r)];
			}
			else
			{
				h0 = real::ZERO;
				h0_adder.clear();
				k0 = max(real::ZERO, r - m + real::ONE);
				for (real k = k0; ; k = k + real::ONE)
				{
					if (k > k0)
					{
						a = a * ((m + k - real::ONE) / (m + k - r - real::ONE));
					}
					else
					{
						if (k == r - m + real::ONE)
						{
							a = real::ONE;
						}
						else
						{
							a = factorial(m - real::ONE) / (factorial(m - real::ONE - r) * factorial(r));
						}
					}
					change = get_c2k(verbose, c, m, n, lambda, n_dr, dr, n_c2k, c2k, k) * (m + real::TWO * k) * a;
					h0 = h0 + change;
					h0_adder.add(change);
					if (k > k0 && abs(change) > real::ZERO && abs(change / h0) < real::SMALL_ENOUGH)
					{
						if (verbose)
						{
							std::cout << "calculate_B2rmn_coefficients: r = " << r.get_int() << ": " << (change / h0).get_string(10) << ", " << ((h0 - h0_adder.calculate_sum()) / h0_adder.calculate_sum()).get_string(10) << ", " << (change / h0_adder.calculate_sum()).get_string(10) << std::endl;
						}
						break;
					}
				}
				h0 = h0_adder.calculate_sum();
				h[gzbi(r - r0)] = -((real::TWO * Q) / k1) * h0;
				h_saved[gzbi(r)] = h[gzbi(r - r0)];
			}
		}
	}
	else
	{
		for (real r = r0; r <= r1; r = r + real::ONE)
		{
			alpha[gzbi(r - r0)] = (real::TWO * r + real::ONE) * (real::TWO * r + real::TWO);
			beta[gzbi(r - r0)] = real::TWO * r * (real::TWO * r - real::TWO * m + real::ONE) + m * (m - real::ONE) - lambda;
			gamma[gzbi(r - r0)] = c * c;
			if (h_saved[gzbi(r)] == h_saved[gzbi(r)])
			{
				h[gzbi(r - r0)] = h_saved[gzbi(r)];
			}
			else
			{
				h0 = real::ZERO;
				h0_adder.clear();
				k0 = max(real::ZERO, r - m);
				for (real k = k0; ; k = k + real::ONE)
				{
					if (k > k0)
					{
						a = a * ((m + k) / (m + k - r));
					}
					else
					{
						if (k == r - m)
						{
							a = real::ONE;
						}
						else
						{
							a = factorial(m) / (factorial(m - r) * factorial(r));
						}
					}
					change = get_c2k(verbose, c, m, n, lambda, n_dr, dr, n_c2k, c2k, k) * (m + real::TWO * k + real::ONE) * a;
					h0 = h0 + change;
					h0_adder.add(change);
					if (k > k0 && abs(change) > real::ZERO && abs(change / h0) < real::SMALL_ENOUGH)
					{
						if (verbose)
						{
							std::cout << "calculate_B2rmn_coefficients: r = " << r.get_int() << ": " << (change / h0).get_string(10) << ", " << ((h0 - h0_adder.calculate_sum()) / h0_adder.calculate_sum()).get_string(10) << ", " << (change / h0_adder.calculate_sum()).get_string(10) << std::endl;
						}
						break;
					}
				}
				h0 = h0_adder.calculate_sum();
				h1 = real::ZERO;
				h1_adder.clear();
				k0 = max(real::ZERO, r - m + real::ONE);
				for (real k = k0; ; k = k + real::ONE)
				{
					if (k > k0)
					{
						a = a * ((m + k - real::ONE) / (m + k - r - real::ONE));
					}
					else
					{
						if (k == r - m + real::ONE)
						{
							a = real::ONE;
						}
						else
						{
							a = factorial(m - real::ONE) / (factorial(m - real::ONE - r) * factorial(r));
						}
					}
					change = get_c2k(verbose, c, m, n, lambda, n_dr, dr, n_c2k, c2k, k) * (m + real::TWO * k) * a;
					h1 = h1 + change;
					h1_adder.add(change);
					if (k > k0 && abs(change) > real::ZERO && abs(change / h1) < real::SMALL_ENOUGH)
					{
						if (verbose)
						{
							std::cout << "calculate_B2rmn_coefficients: r = " << r.get_int() << ": " << (change / h1).get_string(10) << ", " << ((h1 - h1_adder.calculate_sum()) / h1_adder.calculate_sum()).get_string(10) << ", " << (change / h1_adder.calculate_sum()).get_string(10) << std::endl;
						}
						break;
					}
				}
				h1 = h1_adder.calculate_sum();
				h[gzbi(r - r0)] = -((real::TWO * Q) / k1) * (h0 - h1);
				h_saved[gzbi(r)] = h[gzbi(r - r0)];
			}
		}
	}
}

static void calculate_B2rmn_forward(std::vector<real> & alpha, std::vector<real> & beta, std::vector<real> & gamma, std::vector<real> & h, std::vector<real> & B2r, bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr, const real & k1, real & n_c2k, std::vector<real> & c2k, const real & Q, std::vector<real> & h_saved, const real & n_B2r)
{
	real R1;
	real R1p;
	
	calculate_B2rmn_coefficients(alpha, beta, gamma, h, verbose, c, m, n, lambda, n_dr, dr, k1, n_c2k, c2k, Q, h_saved, real::ZERO, n_B2r - real::ONE);
	B2r.clear();
	for (real r = real::ZERO; r <= n_B2r - real::ONE; r = r + real::ONE)
	{
		B2r.push_back(real::ZERO);
	}
	calculate_Rmn1_2(R1, R1p, verbose, c, m, n, k1, n_c2k, c2k, real::ZERO);
	if (remainder(n - m, real::TWO) == real::ZERO)
	{
		B2r[gzbi(real::ZERO)] = pow(c * R1, -real::ONE) - Q * R1;
	}
	else
	{
		B2r[gzbi(real::ZERO)] = -pow(c * R1p, -real::ONE);
	}
	for (real r = real::ZERO; r <= n_B2r - real::TWO; r = r + real::ONE)
	{
		if (r > real::ZERO)
		{
			B2r[gzbi(r + real::ONE)] = (h[gzbi(r)] - beta[gzbi(r)] * B2r[gzbi(r)] - gamma[gzbi(r)] * B2r[gzbi(r - real::ONE)]) / alpha[gzbi(r)];
		}
		else
		{
			B2r[gzbi(r + real::ONE)] = (h[gzbi(r)] - beta[gzbi(r)] * B2r[gzbi(r)]) / alpha[gzbi(r)];
		}
	}
}

static void calculate_B2rmn_backward(std::vector<real> & B2r, bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr, const real & k1, real & n_c2k, std::vector<real> & c2k, const real & Q, std::vector<real> & h_saved, const real & n_B2r, std::vector<real> & alpha, std::vector<real> & beta, std::vector<real> & gamma, std::vector<real> & h, const real & r0, const real & B0)
{
	std::vector<real> p;
	std::vector<real> e;
	std::vector<real> prev_B2r;
	real prev_n_B2r_more;
	std::vector<real> alpha_more;
	std::vector<real> beta_more;
	std::vector<real> gamma_more;
	std::vector<real> h_more;
	real max_abs_change;
	
	for (real n_B2r_more = n_B2r; ; n_B2r_more = n_B2r_more + real("100"))
	{
		if (n_B2r_more > n_B2r)
		{
			calculate_B2rmn_coefficients(alpha_more, beta_more, gamma_more, h_more, verbose, c, m, n, lambda, n_dr, dr, k1, n_c2k, c2k, Q, h_saved, prev_n_B2r_more, n_B2r_more - real::ONE);
			alpha.insert(alpha.end(), alpha_more.begin(), alpha_more.end());
			beta.insert(beta.end(), beta_more.begin(), beta_more.end());
			gamma.insert(gamma.end(), gamma_more.begin(), gamma_more.end());
			h.insert(h.end(), h_more.begin(), h_more.end());
			for (real r = prev_n_B2r_more; r <= n_B2r_more - real::ONE; r = r + real::ONE)
			{
				p.push_back(real::ZERO);
				e.push_back(real::ZERO);
			}
			for (real r = prev_n_B2r_more; r <= n_B2r_more - real::ONE; r = r + real::ONE)
			{
				p[gzbi(r)] = alpha[gzbi(r)] / (beta[gzbi(r)] - gamma[gzbi(r)] * p[gzbi(r - real::ONE)]);
				e[gzbi(r)] = (h[gzbi(r)] - gamma[gzbi(r)] * e[gzbi(r - real::ONE)]) / (beta[gzbi(r)] - gamma[gzbi(r)] * p[gzbi(r - real::ONE)]);
			}
			for (real r = prev_n_B2r_more; r <= n_B2r_more - real::ONE; r = r + real::ONE)
			{
				B2r.push_back(real::ZERO);
			}
		}
		else
		{
			for (real r = real::ZERO; r <= n_B2r - real::ONE; r = r + real::ONE)
			{
				p.push_back(real::ZERO);
				e.push_back(real::ZERO);
			}
			p[gzbi(r0 + real::ONE)] = alpha[gzbi(r0 + real::ONE)] / beta[gzbi(r0 + real::ONE)];
			e[gzbi(r0 + real::ONE)] = (h[gzbi(r0 + real::ONE)] - gamma[gzbi(r0 + real::ONE)] * B0) / beta[gzbi(r0 + real::ONE)];
			for (real r = r0 + real::TWO; r <= n_B2r - real::ONE; r = r + real::ONE)
			{
				p[gzbi(r)] = alpha[gzbi(r)] / (beta[gzbi(r)] - gamma[gzbi(r)] * p[gzbi(r - real::ONE)]);
				e[gzbi(r)] = (h[gzbi(r)] - gamma[gzbi(r)] * e[gzbi(r - real::ONE)]) / (beta[gzbi(r)] - gamma[gzbi(r)] * p[gzbi(r - real::ONE)]);
			}
			B2r.clear();
			prev_B2r.clear();
			for (real r = real::ZERO; r <= n_B2r - real::ONE; r = r + real::ONE)
			{
				B2r.push_back(real::ZERO);
				prev_B2r.push_back(real::ZERO);
			}
		}
		B2r[gzbi(n_B2r_more - real::ONE)] = e[gzbi(n_B2r_more - real::ONE)];
		for (real r = n_B2r_more - real::TWO; r >= r0 + real::ONE; r = r - real::ONE)
		{
			B2r[gzbi(r)] = e[gzbi(r)] - p[gzbi(r)] * B2r[gzbi(r + real::ONE)];
		}
		B2r[gzbi(r0)] = B0;
		if (n_B2r_more > n_B2r)
		{
			max_abs_change = real::ZERO;
			for (real r = r0; r <= n_B2r - real::ONE; r = r + real::ONE)
			{
				max_abs_change = max(max_abs_change, abs((B2r[gzbi(r)] - prev_B2r[gzbi(r)]) / prev_B2r[gzbi(r)]));
			}
			if (verbose)
			{
				std:: cout << "calculate_B2rmn_backward: " << max_abs_change.get_string(10) << std::endl;
			}
			if (max_abs_change < real::SMALL_ENOUGH)
			{
				break;
			}
		}
		std::copy(B2r.begin() + gzbi(r0), B2r.begin() + gzbi(n_B2r), prev_B2r.begin() + gzbi(r0));
		prev_n_B2r_more = n_B2r_more;
	}
	B2r.erase(B2r.begin() + gzbi(n_B2r), B2r.end());
}

static void calculate_B2rmn_once(std::vector<real> & B2r, bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr, const real & k1, real & n_c2k, std::vector<real> & c2k, const real & Q, std::vector<real> & h_saved, const real & n_B2r)
{
	std::vector<real> alpha;
	std::vector<real> beta;
	std::vector<real> gamma;
	std::vector<real> h;
	std::vector<real> B2r_forward;
	real max_r;
	std::vector<real> B2r_backward;
	
	B2r.clear();
	for (real r = real::ZERO; r <= n_B2r - real::ONE; r = r + real::ONE)
	{
		B2r.push_back(real::ZERO);
	}
	calculate_B2rmn_forward(alpha, beta, gamma, h, B2r_forward, verbose, c, m, n, lambda, n_dr, dr, k1, n_c2k, c2k, Q, h_saved, n_B2r);
	max_r = -real::ONE;
	for (real r = real::ZERO; r <= n_B2r - real::ONE; r = r + real::ONE)
	{
		if (max_r == -real::ONE || abs(B2r_forward[gzbi(r)]) > abs(B2r_forward[gzbi(max_r)]))
		{
			max_r = r;
		}
	}
	std::copy(B2r_forward.begin(), B2r_forward.begin() + gzbi(max_r + real::ONE), B2r.begin());
	if (max_r < n_B2r - 1)
	{
		calculate_B2rmn_backward(B2r_backward, verbose, c, m, n, lambda, n_dr, dr, k1, n_c2k, c2k, Q, h_saved, n_B2r, alpha, beta, gamma, h, max_r, B2r[gzbi(max_r)]);
		std::copy(B2r_backward.begin() + gzbi(max_r + real::ONE), B2r_backward.end(), B2r.begin() + gzbi(max_r + real::ONE));
	}
}

void calculate_B2rmn(std::vector<real> & B2r, bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr, const real & k1, real & n_c2k, std::vector<real> & c2k, const real & Q, real & n_B2r, const real & B2r_min)
{
	real n_B2r_orig;
	std::vector<real> h_saved;
	real remove_where;
	
	n_B2r_orig = n_B2r;
	h_saved.clear();
	for ( ; ; n_B2r = real::TWO * n_B2r)
	{
		calculate_B2rmn_once(B2r, verbose, c, m, n, lambda, n_dr, dr, k1, n_c2k, c2k, Q, h_saved, n_B2r);
		if (B2r_min == real::ZERO || abs(B2r[gzbi(n_B2r - real::ONE)]) < B2r_min)
		{
			break;
		}
	}
	if (B2r_min > real::ZERO)
	{
		for (real r = n_B2r - real::ONE; r >= real::ZERO; r = r - real::ONE)
		{
			if (abs(B2r[gzbi(r)]) >= B2r_min)
			{
				remove_where = r + real::TWO;
				if (remove_where >= n_B2r_orig && remove_where <= n_B2r - real::ONE)
				{
					B2r.erase(B2r.begin() + gzbi(remove_where), B2r.end());
					n_B2r = real((int)B2r.size());
				}
				break;
			}
		}
	}
}

void calculate_Rmn1_1(real & R1, real & R1p, bool verbose, const real & c, const real & m, const real & n, const real & n_dr, const std::vector<real> & dr, const real & F, const real & xi)
{
	calculate_Rmn1_1_shared(R1, R1p, verbose, c, m, n, n_dr, dr, xi);
	if (m > real::ZERO)
	{
		R1p = pow(F, -real::ONE) * (m / real::TWO) * pow(real::ONE + real::ONE / (xi * xi), m / real::TWO - real::ONE) * (-real::TWO / (xi * xi * xi)) * R1 + pow(F, -real::ONE) * pow(real::ONE + real::ONE / (xi * xi), m / real::TWO) * R1p;
		R1 = pow(F, -real::ONE) * pow(real::ONE + real::ONE / (xi * xi), m / real::TWO) * R1;
	}
	else
	{
		R1p = pow(F, -real::ONE) * R1p;
		R1 = pow(F, -real::ONE) * R1;
	}
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
			a = a * (xi * xi + real::ONE);
		}
		else
		{
			a = real::ONE;
		}
		change = c2k[gzbi(k)] * a;
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
			ap = ap * (xi * xi + real::ONE);
		}
		else
		{
			ap = real::ONE;
		}
		changep = c2k[gzbi(k)] * k * ap * real::TWO * xi;
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
		R1p = pow(k1, -real::ONE) * (m / real::TWO) * pow(xi * xi + real::ONE, m / real::TWO - real::ONE) * real::TWO * xi * R1 + pow(k1, -real::ONE) * pow(xi * xi + real::ONE, m / real::TWO) * R1p;
		R1 = pow(k1, -real::ONE) * pow(xi * xi + real::ONE, m / real::TWO) * R1;
	}
	else
	{
		R1p = pow(k1, -real::ONE) * pow(xi * xi + real::ONE, m / real::TWO) * R1 + pow(k1, -real::ONE) * xi * (m / real::TWO) * pow(xi * xi + real::ONE, m / real::TWO - real::ONE) * real::TWO * xi * R1 + pow(k1, -real::ONE) * xi * pow(xi * xi + real::ONE, m / real::TWO) * R1p;
		R1 = pow(k1, -real::ONE) * xi * pow(xi * xi + real::ONE, m / real::TWO) * R1;
	}
}

void calculate_Rmn2_1(real & R2, real & R2p, bool verbose, const real & c, const real & m, const real & n, const real & n_dr, const std::vector<real> & dr, const real & F, const real & xi)
{
	calculate_Rmn2_1_shared(R2, R2p, verbose, c, m, n, n_dr, dr, xi);
	if (m > real::ZERO)
	{
		R2p = pow(F, -real::ONE) * (m / real::TWO) * pow(real::ONE + real::ONE / (xi * xi), m / real::TWO - real::ONE) * (-real::TWO / (xi * xi * xi)) * R2 + pow(F, -real::ONE) * pow(real::ONE + real::ONE / (xi * xi), m / real::TWO) * R2p;
		R2 = pow(F, -real::ONE) * pow(real::ONE + real::ONE / (xi * xi), m / real::TWO) * R2;
	}
	else
	{
		R2p = pow(F, -real::ONE) * R2p;
		R2 = pow(F, -real::ONE) * R2;
	}
}

static complex calculate_continued_fraction(const complex & b0, const std::vector<complex> & a, const std::vector<complex> & b)
{
	complex x;
	
	x = real::ZERO;
	if ((int)a.size() < (int)b.size())
	{
		for (int i = (int)a.size() - 1; i >= 0; --i)
		{
			x = a[i] / (b[i] - x);
		}
	}
	else
	{
		for (int i = (int)b.size() - 1; i >= 0; --i)
		{
			x = a[i] / (b[i] - x);
		}
	}
	x = b0 - x;
	return x;
}

static void calculate_Q(std::vector<complex> & Q, bool verbose, const real & m0, const real & n1, const real & xi)
{
	complex x;
	complex Q0;
	complex Q1;
	complex Q2;
	complex b0;
	std::vector<complex> a;
	std::vector<complex> b;
	complex prev_N;
	complex N;
	
	x = complex::I * xi;
	Q.clear();
	for (real n = real::ZERO; n <= n1; n = n + real::ONE)
	{
		Q.push_back(real::ZERO);
	}
	Q0 = (real::ONE / real::TWO) * log((x + real::ONE) / (x - real::ONE));
	Q1 = -pow(x * x - real::ONE, -real::ONE / real::TWO);
	if (m0 > real::ONE)
	{
		for (real m = real::TWO; m <= m0; m = m + real::ONE)
		{
			Q2 = (m - real::ONE) * (-m + real::TWO) * Q0 - ((real::TWO * (m - real::ONE)) / pow(x * x - real::ONE, real::ONE / real::TWO)) * x * Q1;
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
			N = (n + m0) / ((real::TWO * n + real::ONE) * x - (n - m0 + real::ONE) * N);
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
				b.push_back((real::TWO * i + real::ONE) * x);
				a.push_back(i - m0 + real::ONE);
				if (i > n && pow(real::TWO, round(log(i - n) / log(real::TWO))) == i - n)
				{
					N = -calculate_continued_fraction(b0, a, b);
					if (prev_N.a == prev_N.a)
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
	complex x;
	real v_max;
	std::vector<complex> P;
	std::vector<complex> Pp;
	std::vector<complex> Q;
	std::vector<complex> Qp;
	complex Q1;
	complex Q2;
	complex R2_complex;
	complex R2p_complex;
	complex_adder R2_adder;
	complex_adder R2p_adder;
	complex change;
	complex changep;
	
	if (xi > real::ZERO)
	{
		x = complex::I * xi;
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
			P[gzbi(v)] = (real::TWO * v - real::ONE) * pow(x * x - real::ONE, real::ONE / real::TWO) * P[gzbi(v - real::ONE)];
			P[gzbi(v - real::ONE)] = real::ZERO;
		}
		P[gzbi(m + real::ONE)] = (real::TWO * m + real::ONE) * x * P[gzbi(m)];
		for (real v = m + real::TWO; v <= v_max; v = v + real::ONE)
		{
			P[gzbi(v)] = (real::ONE / (v - m)) * (-(v + m - real::ONE) * P[gzbi(v - real::TWO)] + (real::TWO * v - real::ONE) * x * P[gzbi(v - real::ONE)]);
		}
		Pp.clear();
		for (real v = real::ZERO; v <= v_max - real::ONE; v = v + real::ONE)
		{
			Pp.push_back(real::ZERO);
		}
		for (real v = m; v <= v_max - real::ONE; v = v + real::ONE)
		{
			Pp[gzbi(v)] = (real::ONE / (x * x - real::ONE)) * (-(v + real::ONE) * x * P[gzbi(v)] + (v - m + real::ONE) * P[gzbi(v + real::ONE)]);
		}
		calculate_Q(Q, verbose, m, m + n_dr, xi);
		Qp.clear();
		for (real v = real::ZERO; v <= m + n_dr - real::ONE; v = v + real::ONE)
		{
			Qp.push_back(real::ZERO);
		}
		for (real v = real::ZERO; v <= m + n_dr - real::ONE; v = v + real::ONE)
		{
			Qp[gzbi(v)] = (real::ONE / (x * x - real::ONE)) * (-(v + 1) * x * Q[gzbi(v)] + (v - m + real::ONE) * Q[gzbi(v + real::ONE)]);
		}
		Q1 = Q[gzbi(real::ZERO)];
		Q2 = Q[gzbi(real::ONE)];
		for (real v = -real::ONE; v >= -m; v = v - real::ONE)
		{
			P[gnobi(v)] = (real::ONE / (v + m + real::ONE)) * ((real::TWO * v + real::THREE) * x * Q1 - (v - m + real::TWO) * Q2);
			Pp[gnobi(v)] = (real::ONE / (x * x - real::ONE)) * (-(v + real::ONE) * x * P[gnobi(v)] + (v - m + real::ONE) * Q1);
			Q2 = Q1;
			Q1 = P[gnobi(v)];
		}
	}
	else
	{
		// For the code at the end of this else section to work properly, there
		// need to be at least as many entries in P as there are in Q.
		v_max = max(m + n_dr - real::ONE, n_dr_neg - m - real::ONE);
		P.clear();
		Pp.clear();
		for (real v = real::ZERO; v <= v_max; v = v + real::ONE)
		{
			P.push_back(real::ZERO);
			Pp.push_back(real::ZERO);
		}
		for (real v = m; v <= v_max; v = v + real::ONE)
		{
			if (remainder(v - m, real::TWO) == real::ZERO)
			{
				P[gzbi(v)] = (pow(-real::ONE, (v - m) / real::TWO) * factorial(v + m)) / (pow(real::TWO, v) * factorial((v + m) / real::TWO) * factorial((v - m) / real::TWO));
			}
			else
			{
				Pp[gzbi(v)] = (pow(-real::ONE, (v - m - real::ONE) / real::TWO) * factorial(v + m + real::ONE)) / (pow(real::TWO, v) * factorial((v + m + real::ONE) / real::TWO) * factorial((v - m - real::ONE) / real::TWO));
			}
		}
		Q.clear();
		Qp.clear();
		for (real v = real::ZERO; v <= m + n_dr - real::ONE; v = v + real::ONE)
		{
			Q.push_back(real::ZERO);
			Qp.push_back(real::ZERO);
		}
		for (real v = m; v <= m + n_dr - real::ONE; v = v + real::ONE)
		{
			if (remainder(v - m, real::TWO) == real::ZERO)
			{
				Qp[gzbi(v)] = (pow(-real::ONE, (v - m) / real::TWO) * pow(real::TWO, v - real::ONE) * factorial((v + m) / real::TWO) * factorial((v - m - real::TWO) / real::TWO)) / factorial(v - m - real::ONE);
			}
			else
			{
				Q[gzbi(v)] = (pow(-real::ONE, (v - m + real::ONE) / real::TWO) * pow(real::TWO, v - real::ONE) * factorial((v + m - real::ONE) / real::TWO) * factorial((v - m - real::ONE) / real::TWO)) / factorial(v - m);
			}
		}
		for (real v = m - real::ONE; v >= real::ZERO; v = v - real::ONE)
		{
			Q[gzbi(v)] = -((v - m + real::TWO) / (v + m + real::ONE)) * Q[gzbi(v + real::TWO)];
			Qp[gzbi(v)] = -(v - m + real::ONE) * Q[gzbi(v + real::ONE)];
		}
		Q1 = Q[gzbi(real::ZERO)];
		Q2 = Q[gzbi(real::ONE)];
		for (real v = -real::ONE; v >= -m; v = v - real::ONE)
		{
			P[gnobi(v)] = -((v - m + real::TWO) / (v + m + real::ONE)) * Q2;
			Pp[gnobi(v)] = -(v - m + real::ONE) * Q1;
			Q2 = Q1;
			Q1 = P[gnobi(v)];
		}
		for (real v = real::ZERO; v <= m - real::ONE; v = v + real::ONE)
		{
			Q[gzbi(v)] = pow(complex::I, m) * Q[gzbi(v)];
			Qp[gzbi(v)] = pow(complex::I, m) * Qp[gzbi(v)];
		}
		for (real v = m; v <= m + n_dr - real::ONE; v = v + real::ONE)
		{
			Q[gzbi(v)] = pow(complex::I, m) * (Q[gzbi(v)] - (real::ONE / real::TWO) * real::PI * complex::I * P[gzbi(v)]);
			Qp[gzbi(v)] = pow(complex::I, m) * (Qp[gzbi(v)] - (real::ONE / real::TWO) * real::PI * complex::I * Pp[gzbi(v)]);
		}
		for (real v = -real::ONE; v >= -m; v = v - real::ONE)
		{
			P[gnobi(v)] = pow(complex::I, m) * P[gnobi(v)];
			Pp[gnobi(v)] = pow(complex::I, m) * Pp[gnobi(v)];
		}
		for (real v = m; v <= v_max; v = v + real::ONE)
		{
			P[gzbi(v)] = pow(complex::I, m) * P[gzbi(v)];
			Pp[gzbi(v)] = pow(complex::I, m) * Pp[gzbi(v)];
		}
	}
	R2_complex = real::ZERO;
	R2_adder.clear();
	R2p_complex = real::ZERO;
	R2p_adder.clear();
	if (remainder(n - m, real::TWO) == real::ZERO)
	{
		for (real r = -real::TWO * m - real::TWO; r >= -n_dr_neg; r = r - real::TWO)
		{
			change = dr_neg[gnobi(r)] * P[gzbi(-r - m - real::ONE)];
			R2_complex = R2_complex + change;
			R2_adder.add(change);
			changep = dr_neg[gnobi(r)] * Pp[gzbi(-r - m - real::ONE)] * complex::I;
			R2p_complex = R2p_complex + changep;
			R2p_adder.add(changep);
			if (r < -real::TWO * m - real::TWO && abs(change) > real::ZERO && abs(change / R2_complex) < real::SMALL_ENOUGH && abs(changep) > real::ZERO && abs(changep / R2p_complex) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Rmn2_2: " << (change / R2_complex).get_string(10) << ", " << (changep / R2p_complex).get_string(10) << ", " << ((R2_complex - R2_adder.calculate_sum()) / R2_adder.calculate_sum()).get_string(10) << ", " << ((R2p_complex - R2p_adder.calculate_sum()) / R2p_adder.calculate_sum()).get_string(10) << ", " << (change / R2_adder.calculate_sum()).get_string(10) << ", " << (changep / R2p_adder.calculate_sum()).get_string(10) << std::endl;
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
				changep = dr_neg[gnobi(r)] * Qp[gzbi(m + r)] * complex::I;
			}
			else
			{
				change = dr_neg[gnobi(r)] * P[gnobi(m + r)];
				changep = dr_neg[gnobi(r)] * Pp[gnobi(m + r)] * complex::I;
			}
			R2_complex = R2_complex + change;
			R2_adder.add(change);
			R2p_complex = R2p_complex + changep;
			R2p_adder.add(changep);
			if (r < -real::TWO && abs(change) > real::ZERO && abs(change / R2_complex) < real::SMALL_ENOUGH && abs(changep) > real::ZERO && abs(changep / R2p_complex) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Rmn2_2: " << (change / R2_complex).get_string(10) << ", " << (changep / R2p_complex).get_string(10) << ", " << ((R2_complex - R2_adder.calculate_sum()) / R2_adder.calculate_sum()).get_string(10) << ", " << ((R2p_complex - R2p_adder.calculate_sum()) / R2p_adder.calculate_sum()).get_string(10) << ", " << (change / R2_adder.calculate_sum()).get_string(10) << ", " << (changep / R2p_adder.calculate_sum()).get_string(10) << std::endl;
				}
				break;
			}
		}
		for (real r = real::ZERO; r <= n_dr - real::TWO; r = r + real::TWO)
		{
			change = dr[gzbi(r)] * Q[gzbi(m + r)];
			R2_complex = R2_complex + change;
			R2_adder.add(change);
			changep = dr[gzbi(r)] * Qp[gzbi(m + r)] * complex::I;
			R2p_complex = R2p_complex + changep;
			R2p_adder.add(changep);
			if (r > real::ZERO && abs(change) > real::ZERO && abs(change / R2_complex) < real::SMALL_ENOUGH && abs(changep) > real::ZERO && abs(changep / R2p_complex) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Rmn2_2: " << (change / R2_complex).get_string(10) << ", " << (changep / R2p_complex).get_string(10) << ", " << ((R2_complex - R2_adder.calculate_sum()) / R2_adder.calculate_sum()).get_string(10) << ", " << ((R2p_complex - R2p_adder.calculate_sum()) / R2p_adder.calculate_sum()).get_string(10) << ", " << (change / R2_adder.calculate_sum()).get_string(10) << ", " << (changep / R2p_adder.calculate_sum()).get_string(10) << std::endl;
				}
				break;
			}
		}
		R2_complex = R2_adder.calculate_sum();
		R2p_complex = R2p_adder.calculate_sum();
		R2_complex = pow(k2 / pow(-complex::I, m - real::ONE), -real::ONE) * R2_complex;
		R2p_complex = pow(k2 / pow(-complex::I, m - real::ONE), -real::ONE) * R2p_complex;
	}
	else
	{
		for (real r = -real::TWO * m - real::ONE; r >= -n_dr_neg + real::ONE; r = r - real::TWO)
		{
			change = dr_neg[gnobi(r)] * P[gzbi(-r - m - real::ONE)];
			R2_complex = R2_complex + change;
			R2_adder.add(change);
			changep = dr_neg[gnobi(r)] * Pp[gzbi(-r - m - real::ONE)] * complex::I;
			R2p_complex = R2p_complex + changep;
			R2p_adder.add(changep);
			if (r < -real::TWO * m - real::ONE && abs(change) > real::ZERO && abs(change / R2_complex) < real::SMALL_ENOUGH && abs(changep) > real::ZERO && abs(changep / R2p_complex) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Rmn2_2: " << (change / R2_complex).get_string(10) << ", " << (changep / R2p_complex).get_string(10) << ", " << ((R2_complex - R2_adder.calculate_sum()) / R2_adder.calculate_sum()).get_string(10) << ", " << ((R2p_complex - R2p_adder.calculate_sum()) / R2p_adder.calculate_sum()).get_string(10) << ", " << (change / R2_adder.calculate_sum()).get_string(10) << ", " << (changep / R2p_adder.calculate_sum()).get_string(10) << std::endl;
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
				changep = dr_neg[gnobi(r)] * Qp[gzbi(m + r)] * complex::I;
			}
			else
			{
				change = dr_neg[gnobi(r)] * P[gnobi(m + r)];
				changep = dr_neg[gnobi(r)] * Pp[gnobi(m + r)] * complex::I;
			}
			R2_complex = R2_complex + change;
			R2_adder.add(change);
			R2p_complex = R2p_complex + changep;
			R2p_adder.add(changep);
			if (r < -real::ONE && abs(change) > real::ZERO && abs(change / R2_complex) < real::SMALL_ENOUGH && abs(changep) > real::ZERO && abs(changep / R2p_complex) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Rmn2_2: " << (change / R2_complex).get_string(10) << ", " << (changep / R2p_complex).get_string(10) << ", " << ((R2_complex - R2_adder.calculate_sum()) / R2_adder.calculate_sum()).get_string(10) << ", " << ((R2p_complex - R2p_adder.calculate_sum()) / R2p_adder.calculate_sum()).get_string(10) << ", " << (change / R2_adder.calculate_sum()).get_string(10) << ", " << (changep / R2p_adder.calculate_sum()).get_string(10) << std::endl;
				}
				break;
			}
		}
		for (real r = real::ONE; r <= n_dr - real::ONE; r = r + real::TWO)
		{
			change = dr[gzbi(r)] * Q[gzbi(m + r)];
			R2_complex = R2_complex + change;
			R2_adder.add(change);
			changep = dr[gzbi(r)] * Qp[gzbi(m + r)] * complex::I;
			R2p_complex = R2p_complex + changep;
			R2p_adder.add(changep);
			if (r > real::ONE && abs(change) > real::ZERO && abs(change / R2_complex) < real::SMALL_ENOUGH && abs(changep) > real::ZERO && abs(changep / R2p_complex) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Rmn2_2: " << (change / R2_complex).get_string(10) << ", " << (changep / R2p_complex).get_string(10) << ", " << ((R2_complex - R2_adder.calculate_sum()) / R2_adder.calculate_sum()).get_string(10) << ", " << ((R2p_complex - R2p_adder.calculate_sum()) / R2p_adder.calculate_sum()).get_string(10) << ", " << (change / R2_adder.calculate_sum()).get_string(10) << ", " << (changep / R2p_adder.calculate_sum()).get_string(10) << std::endl;
				}
				break;
			}
		}
		R2_complex = R2_adder.calculate_sum();
		R2p_complex = R2p_adder.calculate_sum();
		R2_complex = pow(k2 / pow(-complex::I, m - real::TWO), -real::ONE) * R2_complex;
		R2p_complex = pow(k2 / pow(-complex::I, m - real::TWO), -real::ONE) * R2p_complex;
	}
	R2 = R2_complex.a;
	R2p = R2p_complex.a;
}

void calculate_Rmn2_3(real & R2, real & R2p, bool verbose, const real & c, const real & m, const real & n, const real & Q, const real & n_B2r, const std::vector<real> & B2r, const real & xi, const real & R1, const real & R1p)
{
	adder R2_adder;
	real a;
	real change;
	adder R2p_adder;
	real ap;
	real changep;
	
	R2 = real::ZERO;
	R2_adder.clear();
	for (real r = real::ZERO; r <= n_B2r - real::ONE; r = r + real::ONE)
	{
		if (r > real::ZERO)
		{
			a = a * xi * xi;
		}
		else
		{
			a = real::ONE;
		}
		change = B2r[gzbi(r)] * a;
		R2 = R2 + change;
		R2_adder.add(change);
		if (r > real::ZERO && abs(change) > real::ZERO && abs(change / R2) < real::SMALL_ENOUGH)
		{
			if (verbose)
			{
				std::cout << "calculate_Rmn2_3: " << (change / R2).get_string(10) << ", " << ((R2 - R2_adder.calculate_sum()) / R2_adder.calculate_sum()).get_string(10) << ", " << (change / R2_adder.calculate_sum()).get_string(10) << std::endl;
			}
			break;
		}
	}
	R2 = R2_adder.calculate_sum();
	R2p = real::ZERO;
	R2p_adder.clear();
	for (real r = real::ONE; r <= n_B2r - real::ONE; r = r + real::ONE)
	{
		if (r > real::ONE)
		{
			ap = ap * xi * xi;
		}
		else
		{
			ap = xi;
		}
		changep = B2r[gzbi(r)] * real::TWO * r * ap;
		R2p = R2p + changep;
		R2p_adder.add(changep);
		if (r > real::ONE && abs(changep) > real::ZERO && abs(changep / R2p) < real::SMALL_ENOUGH)
		{
			if (verbose)
			{
				std::cout << "calculate_Rmn2_3: " << (changep / R2p).get_string(10) << ", " << ((R2p - R2p_adder.calculate_sum()) / R2p_adder.calculate_sum()).get_string(10) << ", " << (changep / R2p_adder.calculate_sum()).get_string(10) << std::endl;
			}
			break;
		}
	}
	R2p = R2p_adder.calculate_sum();
	if (remainder(n - m, real::TWO) == real::ZERO)
	{
		R2p = pow(xi * xi + real::ONE, -m / real::TWO) * R2 + xi * (-m / real::TWO) * pow(xi * xi + real::ONE, -m / real::TWO - real::ONE) * real::TWO * xi * R2 + xi * pow(xi * xi + real::ONE, -m / real::TWO) * R2p;
		R2 = xi * pow(xi * xi + real::ONE, -m / real::TWO) * R2;
	}
	else
	{
		R2p = (-m / real::TWO) * pow(xi * xi + real::ONE, -m / real::TWO - real::ONE) * real::TWO * xi * R2 + pow(xi * xi + real::ONE, -m / real::TWO) * R2p;
		R2 = pow(xi * xi + real::ONE, -m / real::TWO) * R2;
	}
	R2 = Q * R1 * (atan(xi) - real::PI / real::TWO) + R2;
	R2p = Q * R1p * (atan(xi) - real::PI / real::TWO) + Q * R1 * (real::ONE / (xi * xi + real::ONE)) + R2p;
}
