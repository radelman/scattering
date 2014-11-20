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
#include "real.hpp"
#include <vector>

static real calculate_alphar(const real & c, const real & m, const real & r);
static real calculate_betar(const real & c, const real & m, const real & r);
static real calculate_gammar(const real & c, const real & m, const real & r);
static real calculate_betarm(const real & c, const real & m, const real & r);
static real calculate_gammarm(const real & c, const real & m, const real & r);
static real calculate_U(bool verbose, const real & c, const real & m, const real & n, const real & lambda);
static void calculate_zero(real & x, real & Ux, bool verbose, const real & c, const real & m, const real & n, real a, real Ua, real b, real Ub);
static real calculate_Nrm(bool verbose, const real & c, const real & m, const real & r, const real & lambda);
static real calculate_Arm(const real & c, const real & m, const real & r);
static real calculate_Brm(const real & c, const real & m, const real & lambda, const real & r);
static real calculate_Crm(const real & c, const real & m, const real & r);
static real get_dr(bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr, const real & r);
static real get_dr_neg(bool verbose, const real & c, const real & m, const real & n, const real & lambda, const real & n_dr, const std::vector<real> & dr, real & n_dr_neg, std::vector<real> & dr_neg, const real & r);

static real calculate_alphar(const real & c, const real & m, const real & r)
{
	return (((real::TWO * m + r + real::TWO) * (real::TWO * m + r + real::ONE)) / ((real::TWO * m + real::TWO * r + real::FIVE) * (real::TWO * m + real::TWO * r + real::THREE))) * calculate_c_squared(c);
}

static real calculate_betar(const real & c, const real & m, const real & r)
{
	return (m + r) * (m + r + real::ONE) + ((real::TWO * (m + r) * (m + r + real::ONE) - real::TWO * m * m - real::ONE) / ((real::TWO * m + real::TWO * r - real::ONE) * (real::TWO * m + real::TWO * r + real::THREE))) * calculate_c_squared(c);
}

static real calculate_gammar(const real & c, const real & m, const real & r)
{
	return ((r * (r - real::ONE)) / ((real::TWO * m + real::TWO * r - real::THREE) * (real::TWO * m + real::TWO * r - real::ONE))) * calculate_c_squared(c);
}

static real calculate_betarm(const real & c, const real & m, const real & r)
{
	return calculate_gammar(c, m, r) * calculate_alphar(c, m, r - real::TWO);
}

static real calculate_gammarm(const real & c, const real & m, const real & r)
{
	return calculate_betar(c, m, r);
}

real calculate_continued_fraction(const real & b0, const std::vector<real> & a, const std::vector<real> & b)
{
	real x;
	
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

static real calculate_U(bool verbose, const real & c, const real & m, const real & n, const real & lambda)
{
	real r;
	real b0;
	std::vector<real> a;
	std::vector<real> b;
	real U1;
	real prev_U2;
	real U2;
	real U;
	
	r = n - m;
	b0 = calculate_gammarm(c, m, r) - lambda;
	a.clear();
	b.clear();
	if (remainder(n - m, real::TWO) == real::ZERO)
	{
		for (real i = r; i >= real::TWO; i = i - real::TWO)
		{
			a.push_back(calculate_betarm(c, m, i));
			b.push_back(calculate_gammarm(c, m, i - real::TWO) - lambda);
		}
	}
	else
	{
		for (real i = r; i >= real::THREE; i = i - real::TWO)
		{
			a.push_back(calculate_betarm(c, m, i));
			b.push_back(calculate_gammarm(c, m, i - real::TWO) - lambda);
		}
	}
	U1 = calculate_continued_fraction(b0, a, b);
	b0 = real::ZERO;
	a.clear();
	b.clear();
	prev_U2 = real::NAN;
	for (real i = r + real::TWO; ; i = i + real::TWO)
	{
		a.push_back(calculate_betarm(c, m, i));
		b.push_back(calculate_gammarm(c, m, i) - lambda);
		if (remainder(i - (r + real::TWO), real("100.0")) == real::ZERO)
		{
			U2 = calculate_continued_fraction(b0, a, b);
			if (prev_U2 == prev_U2)
			{
				if (abs((U2 - prev_U2) / prev_U2) < real::SMALL_ENOUGH)
				{
					if (verbose)
					{
						std::cout << "calculate_U: " << ((U2 - prev_U2) / prev_U2).get_string(10) << std::endl;
					}
					break;
				}
			}
			prev_U2 = U2;
		}
	}
	U = U1 + U2;
	return U;
}

static void calculate_zero(real & x, real & Ux, bool verbose, const real & c, const real & m, const real & n, real a, real Ua, real b, real Ub)
{
	real prev_x;
	
	for (int i = 0; ; ++i)
	{
		x = a - (Ua / (Ub - Ua)) * (b - a);
		Ux = calculate_U(verbose, c, m, n, x);
		if (Ux != real::ZERO)
		{
			if (Ua < real::ZERO && Ub > real::ZERO)
			{
				if (Ux < real::ZERO)
				{
					a = x;
					Ua = Ux;
				}
				else
				{
					b = x;
					Ub = Ux;
				}
			}
			else
			{
				if (Ux > real::ZERO)
				{
					a = x;
					Ua = Ux;
				}
				else
				{
					b = x;
					Ub = Ux;
				}
			}
		}
		else
		{
			std::cout << "calculate_zero: U = 0.0" << std::endl;
			return;
		}
		if (i > 0)
		{
			if (verbose)
			{
				std::cout << "calculate_zero: " << ((x - prev_x) / prev_x).get_string(10) << std::endl;
			}
			if (abs((x - prev_x) / prev_x) < real::SMALL_ENOUGH)
			{
				break;
			}
		}
		prev_x = x;
	}
}

void calculate_lambdamn(real & lambda, bool verbose, const real & c, const real & m, const real & n, const real & lambda_approx)
{
	real x;
	real Ux;
	real d;
	real a;
	real Ua;
	real b;
	real Ub;
	
	x = lambda_approx;
	Ux = calculate_U(verbose, c, m, n, x);
	d = pow(real::TWO, -real("100.0")) * x;
	while (true)
	{
		a = x - d;
		Ua = calculate_U(verbose, c, m, n, a);
		b = x + d;
		Ub = calculate_U(verbose, c, m, n, b);
		if (verbose)
		{
			std::cout << "calculate_lambdamn: " << Ua.get_string(10) << ", " << Ub.get_string(10) << std::endl;
		}
		if ((Ua < real::ZERO && Ub > real::ZERO) || (Ua > real::ZERO && Ub < real::ZERO))
		{
			break;
		}
		d = real::TWO * d;
	}
	calculate_zero(x, Ux, verbose, c, m, n, a, Ua, b, Ub);
	lambda = x;
}

static real calculate_Nrm(bool verbose, const real & c, const real & m, const real & r, const real & lambda)
{
	real b0;
	std::vector<real> a;
	std::vector<real> b;
	real N;
	real prev_N;
	
	b0 = real::ZERO;
	a.clear();
	b.clear();
	for (real i = r; ; i = i + real::TWO)
	{
		a.push_back(calculate_betarm(c, m, i));
		b.push_back(calculate_gammarm(c, m, i) - lambda);
		N = -calculate_continued_fraction(b0, a, b);
		if (i > r)
		{
			if (abs((N - prev_N) / prev_N) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Nrm: " << ((N - prev_N) / prev_N).get_string(10) << std::endl;
				}
				break;
			}
		}
		prev_N = N;
	}
	return N;
}

void calculate_drmn(std::vector<real> & dr, bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, const real & dr_min)
{
	real n_dr_orig;
	real N;
	bool converged;
	real x;
	adder x_adder;
	real a;
	real change;
	real s;
	real remove_where;
	
	n_dr_orig = n_dr;
	for ( ; ; n_dr = real::TWO * n_dr)
	{
		if (n_dr > n_dr_orig)
		{
			for (real r = n_dr / real::TWO; r <= n_dr - real::ONE; r = r + real::ONE)
			{
				dr.push_back(real::ZERO);
			}
		}
		else
		{
			dr.clear();
			for (real r = real::ZERO; r <= n_dr_orig - real::ONE; r = r + real::ONE)
			{
				dr.push_back(real::ZERO);
			}
		}
		if (remainder(n - m, real::TWO) == real::ZERO)
		{
			dr[gzbi(n_dr - real::TWO)] = real::ONE;
			for (real r = n_dr - real::TWO; r >= real::TWO; r = r - real::TWO)
			{
				if (r < n_dr - real::TWO)
				{
					N = calculate_betarm(c, m, r) / (calculate_gammarm(c, m, r) - lambda - N);
				}
				else
				{
					N = calculate_Nrm(verbose, c, m, r, lambda);
				}
				dr[gzbi(r - real::TWO)] = -(calculate_alphar(c, m, r - real::TWO) / N) * dr[gzbi(r)];
			}
			converged = false;
			x = real::ZERO;
			x_adder.clear();
			for (real r = real::ZERO; r <= n_dr - real::TWO; r = r + real::TWO)
			{
				if (r > real::ZERO)
				{
					a = a * ((-real::ONE * (real::TWO * m + r - real::ONE) * (real::TWO * m + r)) / (real::FOUR * (m + r / real::TWO) * (r / real::TWO)));
				}
				else
				{
					a = factorial(real::TWO * m) / factorial(m);
				}
				change = dr[gzbi(r)] * a;
				x = x + change;
				x_adder.add(change);
				if (r > real::ZERO && abs(change) > real::ZERO && abs(change / x) < real::SMALL_ENOUGH)
				{
					converged = true;
					if (verbose)
					{
						std::cout << "calculate_drmn: " << (change / x).get_string(10) << ", " << ((x - x_adder.calculate_sum()) / x_adder.calculate_sum()).get_string(10) << ", " << (change / x_adder.calculate_sum()).get_string(10) << std::endl;
					}
					break;
				}
			}
			x = x_adder.calculate_sum();
			if (!converged)
			{
				if (verbose)
				{
					std::cout << "calculate_drmn: warning: x did not converge" << std::endl;
				}
			}
			s = ((pow(-real::ONE, (n - m) / real::TWO) * factorial(n + m)) / (pow(real::TWO, n - m) * factorial((n + m) / real::TWO) * factorial((n - m) / real::TWO))) / x;
			for (real r = real::ZERO; r <= n_dr - real::TWO; r = r + real::TWO)
			{
				dr[gzbi(r)] = s * dr[gzbi(r)];
			}
			if (converged && (dr_min == real::ZERO || abs(dr[gzbi(n_dr - real::TWO)]) < dr_min))
			{
				break;
			}
		}
		else
		{
			dr[gzbi(n_dr - real::ONE)] = real::ONE;
			for (real r = n_dr - real::ONE; r >= real::THREE; r = r - real::TWO)
			{
				if (r < n_dr - real::ONE)
				{
					N = calculate_betarm(c, m, r) / (calculate_gammarm(c, m, r) - lambda - N);
				}
				else
				{
					N = calculate_Nrm(verbose, c, m, r, lambda);
				}
				dr[gzbi(r - real::TWO)] = -(calculate_alphar(c, m, r - real::TWO) / N) * dr[gzbi(r)];
			}
			converged = false;
			x = real::ZERO;
			x_adder.clear();
			for (real r = real::ONE; r <= n_dr - real::ONE; r = r + real::TWO)
			{
				if (r > real::ONE)
				{
					a = a * ((-real::ONE * (real::TWO * m + r) * (real::TWO * m + r + real::ONE)) / (real::FOUR * (m + r / real::TWO + real::ONE / real::TWO) * (r / real::TWO - real::ONE / real::TWO)));
				}
				else
				{
					a = factorial(real::TWO * m + real::TWO) / (real::TWO * factorial(m + real::ONE));
				}
				change = dr[gzbi(r)] * a;
				x = x + change;
				x_adder.add(change);
				if (r > real::ZERO && abs(change) > real::ZERO && abs(change / x) < real::SMALL_ENOUGH)
				{
					converged = true;
					if (verbose)
					{
						std::cout << "calculate_drmn: " << (change / x).get_string(10) << ", " << ((x - x_adder.calculate_sum()) / x_adder.calculate_sum()).get_string(10) << ", " << (change / x_adder.calculate_sum()).get_string(10) << std::endl;
					}
					break;
				}
			}
			x = x_adder.calculate_sum();
			if (!converged)
			{
				if (verbose)
				{
					std::cout << "calculate_drmn: warning: x did not converge" << std::endl;
				}
			}
			s = ((pow(-real::ONE, (n - m - real::ONE) / real::TWO) * factorial(n + m + real::ONE)) / (pow(real::TWO, n - m) * factorial((n + m + real::ONE) / real::TWO) * factorial((n - m - real::ONE) / real::TWO))) / x;
			for (real r = real::ONE; r <= n_dr - real::ONE; r = r + real::TWO)
			{
				dr[gzbi(r)] = s * dr[gzbi(r)];
			}
			if (converged && (dr_min == real::ZERO || abs(dr[gzbi(n_dr - real::ONE)]) < dr_min))
			{
				break;
			}
		}
	}
	if (dr_min > real::ZERO)
	{
		if (remainder(n - m, real::TWO) == real::ZERO)
		{
			for (real r = n_dr - real::TWO; r >= real::ZERO; r = r - real::TWO)
			{
				if (abs(dr[gzbi(r)]) >= dr_min)
				{
					remove_where = r + real::FOUR;
					if (remove_where >= n_dr_orig && remove_where <= n_dr - real::TWO)
					{
						dr.erase(dr.begin() + gzbi(remove_where), dr.end());
						n_dr = real((int)dr.size());
					}
					break;
				}
			}
		}
		else
		{
			for (real r = n_dr - real::ONE; r >= real::ONE; r = r - real::TWO)
			{
				if (abs(dr[gzbi(r)]) >= dr_min)
				{
					remove_where = r + real::THREE;
					if (remove_where >= n_dr_orig && remove_where <= n_dr - real::TWO)
					{
						dr.erase(dr.begin() + gzbi(remove_where), dr.end());
						n_dr = real((int)dr.size());
					}
					break;
				}
			}
		}
	}
}

static real calculate_Arm(const real & c, const real & m, const real & r)
{
	return calculate_alphar(c, m, r - real::TWO);
}

static real calculate_Brm(const real & c, const real & m, const real & lambda, const real & r)
{
	return calculate_betar(c, m, r) - lambda;
}

static real calculate_Crm(const real & c, const real & m, const real & r)
{
	return calculate_gammar(c, m, r + real::TWO);
}

void calculate_drmn_neg(std::vector<real> & dr_neg, bool verbose, const real & c, const real & m, const real & n, const real & lambda, const real & n_dr, const std::vector<real> & dr, real & n_dr_neg, const real & dr_neg_min)
{
	real n_dr_neg_orig;
	real b0;
	std::vector<real> a;
	std::vector<real> b;
	real N;
	real prev_N;
	real s;
	real remove_where;
	
	n_dr_neg_orig = n_dr_neg;
	for ( ; ; n_dr_neg = real::TWO * n_dr_neg)
	{
		if (n_dr_neg > n_dr_neg_orig)
		{
			for (real r = -n_dr_neg / real::TWO - real::ONE; r >= -n_dr_neg; r = r - real::ONE)
			{
				dr_neg.push_back(real::ZERO);
			}
		}
		else
		{
			dr_neg.clear();
			for (real r = -real::ONE; r >= -n_dr_neg_orig; r = r - real::ONE)
			{
				dr_neg.push_back(real::ZERO);
			}
		}
		if (remainder(n - m, real::TWO) == real::ZERO)
		{
			dr_neg[gnobi(-n_dr_neg)] = real::ONE;
			for (real r = -n_dr_neg; r <= -real::TWO; r = r + real::TWO)
			{
				if (r > -n_dr_neg)
				{
					if (r != -real::TWO * m - real::TWO)
					{
						N = -calculate_Arm(c, m, r + real::TWO) / (calculate_Brm(c, m, lambda, r) + calculate_Crm(c, m, r - real::TWO) * N);
					}
					else
					{
						N = (calculate_c_squared(c) / ((real::TWO * m - real::ONE) * (real::TWO * m + real::ONE))) / (calculate_Brm(c, m, lambda, r) + calculate_Crm(c, m, r - real::TWO) * N);
					}
				}
				else
				{
					b0 = real::ZERO;
					a.clear();
					b.clear();
					if (r != -real::TWO * m - real::TWO)
					{
						a.push_back(calculate_Arm(c, m, r + real::TWO));
					}
					else
					{
						a.push_back(-calculate_c_squared(c) / ((real::TWO * m - real::ONE) * (real::TWO * m + real::ONE)));
					}
					b.push_back(calculate_Brm(c, m, lambda, r));
					for (real i = r - real::TWO; ; i = i - real::TWO)
					{
						a.push_back(calculate_Crm(c, m, i) * calculate_Arm(c, m, i + real::TWO));
						b.push_back(calculate_Brm(c, m, lambda, i));
						N = calculate_continued_fraction(b0, a, b);
						if (i < r - real::TWO)
						{
							if (abs((N - prev_N) / prev_N) < real::SMALL_ENOUGH)
							{
								if (verbose)
								{
									std::cout << "calculate_drmn_neg: " << ((N - prev_N) / prev_N).get_string(10) << std::endl;
								}
								break;
							}
						}
						prev_N = N;
					}
				}
				if (r < -real::TWO)
				{
					dr_neg[gnobi(r + real::TWO)] = dr_neg[gnobi(r)] / N;
					if (r == -real::TWO * m - real::TWO)
					{
						N = real::ZERO;
					}
				}
			}
			s = dr[gzbi(real::ZERO)] / (dr_neg[gnobi(-real::TWO)] / N);
			for (real r = -n_dr_neg; r <= -real::TWO; r = r + real::TWO)
			{
				dr_neg[gnobi(r)] = s * dr_neg[gnobi(r)];
			}
			if (dr_neg_min == real::ZERO || abs(dr_neg[gnobi(-n_dr_neg)]) < dr_neg_min)
			{
				break;
			}
		}
		else
		{
			dr_neg[gnobi(-n_dr_neg + real::ONE)] = real::ONE;
			for (real r = -n_dr_neg + real::ONE; r <= -real::ONE; r = r + real::TWO)
			{
				if (r > -n_dr_neg + real::ONE)
				{
					if (r != -real::TWO * m - real::ONE)
					{
						N = -calculate_Arm(c, m, r + real::TWO) / (calculate_Brm(c, m, lambda, r) + calculate_Crm(c, m, r - real::TWO) * N);
					}
					else
					{
						N = -(calculate_c_squared(c) / ((real::TWO * m - real::ONE) * (real::TWO * m - real::THREE))) / (calculate_Brm(c, m, lambda, r) + calculate_Crm(c, m, r - real::TWO) * N);
					}
				}
				else
				{
					b0 = real::ZERO;
					a.clear();
					b.clear();
					if (r != -real::TWO * m - real::ONE)
					{
						a.push_back(calculate_Arm(c, m, r + real::TWO));
					}
					else
					{
						a.push_back(calculate_c_squared(c) / ((real::TWO * m - real::ONE) * (real::TWO * m - real::THREE)));
					}
					b.push_back(calculate_Brm(c, m, lambda, r));
					for (real i = r - real::TWO; ; i = i - real::TWO)
					{
						a.push_back(calculate_Crm(c, m, i) * calculate_Arm(c, m, i + real::TWO));
						b.push_back(calculate_Brm(c, m, lambda, i));
						N = calculate_continued_fraction(b0, a, b);
						if (i < r - real::TWO)
						{
							if (abs((N - prev_N) / prev_N) < real::SMALL_ENOUGH)
							{
								if (verbose)
								{
									std::cout << "calculate_drmn_neg: " << ((N - prev_N) / prev_N).get_string(10) << std::endl;
								}
								break;
							}
						}
						prev_N = N;
					}
				}
				if (r < -real::ONE)
				{
					dr_neg[gnobi(r + real::TWO)] = dr_neg[gnobi(r)] / N;
					if (r == -real::TWO * m - real::ONE)
					{
						N = real::ZERO;
					}
				}
			}
			s = dr[gzbi(real::ONE)] / (dr_neg[gnobi(-real::ONE)] / N);
			for (real r = -n_dr_neg + real::ONE; r <= -real::ONE; r = r + real::TWO)
			{
				dr_neg[gnobi(r)] = s * dr_neg[gnobi(r)];
			}
			if (dr_neg_min == real::ZERO || abs(dr_neg[gnobi(-n_dr_neg + real::ONE)]) < dr_neg_min)
			{
				break;
			}
		}
	}
	if (dr_neg_min > real::ZERO)
	{
		if (remainder(n - m, real::TWO) == real::ZERO)
		{
			for (real r = -n_dr_neg; r <= -real::TWO; r = r + real::TWO)
			{
				if (abs(dr_neg[gnobi(r)]) >= dr_neg_min)
				{
					remove_where = r - real::THREE;
					if (remove_where <= -n_dr_neg_orig - real::ONE && remove_where >= -n_dr_neg + real::ONE)
					{
						dr_neg.erase(dr_neg.begin() + gnobi(remove_where), dr_neg.end());
						n_dr_neg = real((int)dr_neg.size());
					}
					break;
				}
			}
		}
		else
		{
			for (real r = -n_dr_neg + real::ONE; r <= -real::ONE; r = r + real::TWO)
			{
				if (abs(dr_neg[gnobi(r)]) >= dr_neg_min)
				{
					remove_where = r - real::FOUR;
					if (remove_where <= -n_dr_neg_orig - real::ONE && remove_where >= -n_dr_neg + real::ONE)
					{
						dr_neg.erase(dr_neg.begin() + gnobi(remove_where), dr_neg.end());
						n_dr_neg = real((int)dr_neg.size());
					}
					break;
				}
			}
		}
	}
}

static real get_dr(bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr, const real & r)
{
	if (r >= n_dr)
	{
		while (r >= n_dr)
		{
			n_dr = real::TWO * n_dr;
		}
		calculate_drmn(dr, verbose, c, m, n, lambda, n_dr, real::ZERO);
	}
	return dr[gzbi(r)];
}

real calculate_Nmn(bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr)
{
	real N;
	adder N_adder;
	real a;
	real change;
	
	N = real::ZERO;
	N_adder.clear();
	if (remainder(n - m, real::TWO) == real::ZERO)
	{
		for (real r = real::ZERO; ; r = r + real::TWO)
		{
			if (r > real::ZERO)
			{
				a = a * (((real::TWO * m + r - real::ONE) * (real::TWO * m + r)) / ((r - real::ONE) * r));
			}
			else
			{
				a = factorial(real::TWO * m);
			}
			change = get_dr(verbose, c, m, n, lambda, n_dr, dr, r) * get_dr(verbose, c, m, n, lambda, n_dr, dr, r) * (a / (real::TWO * m + real::TWO * r + real::ONE));
			N = N + change;
			N_adder.add(change);
			if (r > real::ZERO && abs(change) > real::ZERO && abs(change / N) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Nmn: " << (change / N).get_string(10) << ", " << ((N - N_adder.calculate_sum()) / N_adder.calculate_sum()).get_string(10) << ", " << (change / N_adder.calculate_sum()).get_string(10) << std::endl;
				}
				break;
			}
		}
	}
	else
	{
		for (real r = real::ONE; ; r = r + real::TWO)
		{
			if (r > real::ONE)
			{
				a = a * (((real::TWO * m + r - real::ONE) * (real::TWO * m + r)) / ((r - real::ONE) * r));
			}
			else
			{
				a = factorial(real::TWO * m + real::ONE);
			}
			change = get_dr(verbose, c, m, n, lambda, n_dr, dr, r) * get_dr(verbose, c, m, n, lambda, n_dr, dr, r) * (a / (real::TWO * m + real::TWO * r + real::ONE));
			N = N + change;
			N_adder.add(change);
			if (r > real::ONE && abs(change) > real::ZERO && abs(change / N) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Nmn: " << (change / N).get_string(10) << ", " << ((N - N_adder.calculate_sum()) / N_adder.calculate_sum()).get_string(10) << ", " << (change / N_adder.calculate_sum()).get_string(10) << std::endl;
				}
				break;
			}
		}
	}
	N = N_adder.calculate_sum();
	N = real::TWO * N;
	return N;
}

real calculate_Fmn(bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr)
{
	real F;
	adder F_adder;
	real a;
	real change;
	
	F = real::ZERO;
	F_adder.clear();
	if (remainder(n - m, real::TWO) == real::ZERO)
	{
		for (real r = real::ZERO; ; r = r + real::TWO)
		{
			if (r > real::ZERO)
			{
				a = a * (((real::TWO * m + r - real::ONE) * (real::TWO * m + r)) / ((r - real::ONE) * r));
			}
			else
			{
				a = factorial(real::TWO * m);
			}
			change = get_dr(verbose, c, m, n, lambda, n_dr, dr, r) * a;
			F = F + change;
			F_adder.add(change);
			if (r > real::ZERO && abs(change) > real::ZERO && abs(change / F) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Fmn: " << (change / F).get_string(10) << ", " << ((F - F_adder.calculate_sum()) / F_adder.calculate_sum()).get_string(10) << ", " << (change / F_adder.calculate_sum()).get_string(10) << std::endl;
				}
				break;
			}
		}
	}
	else
	{
		for (real r = real::ONE; ; r = r + real::TWO)
		{
			if (r > real::ONE)
			{
				a = a * (((real::TWO * m + r - real::ONE) * (real::TWO * m + r)) / ((r - real::ONE) * r));
			}
			else
			{
				a = factorial(real::TWO * m + real::ONE);
			}
			change = get_dr(verbose, c, m, n, lambda, n_dr, dr, r) * a;
			F = F + change;
			F_adder.add(change);
			if (r > real::ONE && abs(change) > real::ZERO && abs(change / F) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Fmn: " << (change / F).get_string(10) << ", " << ((F - F_adder.calculate_sum()) / F_adder.calculate_sum()).get_string(10) << ", " << (change / F_adder.calculate_sum()).get_string(10) << std::endl;
				}
				break;
			}
		}
	}
	F = F_adder.calculate_sum();
	return F;
}

//
// In the oblate case, there's an implicit factor of 1i ^ m when n - m is even
// and 1i ^ (m + 1) when n - m is odd.
//
real calculate_kmn1(bool verbose, const real & c, const real & m, const real & n, const real & n_dr, const std::vector<real> & dr, const real & F)
{
	real k1;
	
	if (remainder(n - m, real::TWO) == real::ZERO)
	{
		k1 = ((real::TWO * m + real::ONE) * factorial(m + n) * F) / (pow(real::TWO, m + n) * dr[gzbi(real::ZERO)] * pow(c, m) * factorial(m) * factorial((n - m) / real::TWO) * factorial((m + n) / real::TWO));
	}
	else
	{
		k1 = ((real::TWO * m + real::THREE) * factorial(m + n + real::ONE) * F) / (pow(real::TWO, m + n) * dr[gzbi(real::ONE)] * pow(c, m + real::ONE) * factorial(m) * factorial((n - m - real::ONE) / real::TWO) * factorial((m + n + real::ONE) / real::TWO));
	}
	return k1;
}

static real get_dr_neg(bool verbose, const real & c, const real & m, const real & n, const real & lambda, const real & n_dr, const std::vector<real> & dr, real & n_dr_neg, std::vector<real> & dr_neg, const real & r)
{
	if (r <= -n_dr_neg - real::ONE)
	{
		while (r <= -n_dr_neg - real::ONE)
		{
			n_dr_neg = real::TWO * n_dr_neg;
		}
		calculate_drmn_neg(dr_neg, verbose, c, m, n, lambda, n_dr, dr, n_dr_neg, real::ZERO);
	}
	return dr_neg[gnobi(r)];
}

real calculate_kmn2(bool verbose, const real & c, const real & m, const real & n, const real & lambda, const real & n_dr, const std::vector<real> & dr, real & n_dr_neg, std::vector<real> & dr_neg, const real & F)
{
	real dr1;
	real k2;
	
	if (remainder(n - m, real::TWO) == real::ZERO)
	{
		if (-real::TWO * m < real::ZERO)
		{
			dr1 = get_dr_neg(verbose, c, m, n, lambda, n_dr, dr, n_dr_neg, dr_neg, -real::TWO * m);
		}
		else
		{
			dr1 = dr[gzbi(real::ZERO)];
		}
		k2 = (pow(real::TWO, n - m) * factorial(real::TWO * m) * factorial((n - m) / real::TWO) * factorial((m + n) / real::TWO) * dr1 * F) / ((real::TWO * m - real::ONE) * factorial(m) * factorial(m + n) * pow(c, m - real::ONE));
	}
	else
	{
		if (-real::TWO * m + real::ONE < real::ONE)
		{
			dr1 = get_dr_neg(verbose, c, m, n, lambda, n_dr, dr, n_dr_neg, dr_neg, -real::TWO * m + real::ONE);
		}
		else
		{
			dr1 = dr[gzbi(real::ONE)];
		}
		k2 = -((pow(real::TWO, n - m) * factorial(real::TWO * m) * factorial((n - m - real::ONE) / real::TWO) * factorial((m + n + real::ONE) / real::TWO) * dr1 * F) / ((real::TWO * m - real::THREE) * (real::TWO * m - real::ONE) * factorial(m) * factorial(m + n + real::ONE) * pow(c, m - real::TWO)));
	}
	return k2;
}

void calculate_c2kmn(std::vector<real> & c2k, bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr, real & n_c2k, const real & c2k_min)
{
	real n_c2k_orig;
	real prev_n_c2k;
	adder c2k_adder;
	real a;
	real a0;
	real change;
	real remove_where;
	
	n_c2k_orig = n_c2k;
	for (prev_n_c2k = real((int)c2k.size()); ; prev_n_c2k = n_c2k, n_c2k = real::TWO * n_c2k)
	{
		for (real k = prev_n_c2k; k <= n_c2k - real::ONE; k = k + real::ONE)
		{
			c2k.push_back(real::ZERO);
		}
		for (real k = prev_n_c2k; k <= n_c2k - real::ONE; k = k + real::ONE)
		{
			c2k[gzbi(k)] = real::ZERO;
			c2k_adder.clear();
			if (remainder(n - m, real::TWO) == real::ZERO)
			{
				for (real r = real::TWO * k; ; r = r + real::TWO)
				{
					if (r > real::TWO * k)
					{
						a = a * (((real::TWO * m + r - real::ONE) * (real::TWO * m + r) * (-r / real::TWO) * (m + r / real::TWO + k - real::ONE / real::TWO)) / ((r - real::ONE) * r * (-r / real::TWO + k) * (m + r / real::TWO - real::ONE / real::TWO)));
						if (r == real::TWO * k + real::TWO)
						{
							a0 = a;
						}
					}
					else
					{
						if (k > prev_n_c2k)
						{
							a = a0 * (-real::ONE) * (m + real::TWO * k - real::ONE / real::TWO);
						}
						else
						{
							a = (factorial(real::TWO * m + r) / factorial(r)) * pochhammer(-r / real::TWO, prev_n_c2k) * pochhammer(m + r / real::TWO + real::ONE / real::TWO, prev_n_c2k);
						}
					}
					change = get_dr(verbose, c, m, n, lambda, n_dr, dr, r) * a;
					c2k[gzbi(k)] = c2k[gzbi(k)] + change;
					c2k_adder.add(change);
					if (r > real::TWO * k && abs(change) > real::ZERO && abs(change / c2k[gzbi(k)]) < real::SMALL_ENOUGH)
					{
						if (verbose)
						{
							std::cout << "calculate_c2kmn: k = " << k.get_int() << ": " << (change / c2k[gzbi(k)]).get_string(10) << ", " << ((c2k[gzbi(k)] - c2k_adder.calculate_sum()) / c2k_adder.calculate_sum()).get_string(10) << ", " << (change / c2k_adder.calculate_sum()).get_string(10) << std::endl;
						}
						break;
					}
				}
			}
			else
			{
				for (real r = real::TWO * k + real::ONE; ; r = r + real::TWO)
				{
					if (r > real::TWO * k + real::ONE)
					{
						a = a * (((real::TWO * m + r - real::ONE) * (real::TWO * m + r) * (-r / real::TWO + real::ONE / real::TWO) * (m + r / real::TWO + k)) / ((r - real::ONE) * r * (-r / real::TWO + k + real::ONE / real::TWO) * (m + r / real::TWO)));
						if (r == real::TWO * k + real::THREE)
						{
							a0 = a;
						}
					}
					else
					{
						if (k > prev_n_c2k)
						{
							a = a0 * (-real::ONE) * (m + real::TWO * k + real::ONE / real::TWO);
						}
						else
						{
							a = (factorial(real::TWO * m + r) / factorial(r)) * pochhammer(-(r - real::ONE) / real::TWO, prev_n_c2k) * pochhammer(m + r / real::TWO + real::ONE, prev_n_c2k);
						}
					}
					change = get_dr(verbose, c, m, n, lambda, n_dr, dr, r) * a;
					c2k[gzbi(k)] = c2k[gzbi(k)] + change;
					c2k_adder.add(change);
					if (r > real::TWO * k + real::ONE && abs(change) > real::ZERO && abs(change / c2k[gzbi(k)]) < real::SMALL_ENOUGH)
					{
						if (verbose)
						{
							std::cout << "calculate_c2kmn: k = " << k.get_int() << ": " << (change / c2k[gzbi(k)]).get_string(10) << ", " << ((c2k[gzbi(k)] - c2k_adder.calculate_sum()) / c2k_adder.calculate_sum()).get_string(10) << ", " << (change / c2k_adder.calculate_sum()).get_string(10) << std::endl;
						}
						break;
					}
				}
			}
			c2k[gzbi(k)] = c2k_adder.calculate_sum();
			c2k[gzbi(k)] = (real::ONE / (pow(real::TWO, m) * factorial(m + k) * factorial(k))) * c2k[gzbi(k)];
		}
		if (c2k_min == real::ZERO || abs(c2k[gzbi(n_c2k - real::ONE)]) < c2k_min)
		{
			break;
		}
	}
	if (c2k_min > real::ZERO)
	{
		for (real k = n_c2k - real::ONE; k >= real::ZERO; k = k - real::ONE)
		{
			if (abs(c2k[gzbi(k)]) >= c2k_min)
			{
				remove_where = k + real::TWO;
				if (remove_where >= n_c2k_orig && remove_where <= n_c2k - real::ONE)
				{
					c2k.erase(c2k.begin() + gzbi(remove_where), c2k.end());
					n_c2k = real((int)c2k.size());
				}
				break;
			}
		}
	}
}

void calculate_Smn1_1(real & S1, real & S1p, bool verbose, const real & c, const real & m, const real & n, const real & n_dr, const std::vector<real> & dr, const real & eta)
{
	adder S1_adder;
	adder S1p_adder;
	real P0;
	real P1;
	real P0p;
	real change;
	real changep;
	real P1p;
	
	S1 = real::ZERO;
	S1_adder.clear();
	S1p = real::ZERO;
	S1p_adder.clear();
	P0 = real::ONE;
	for (real v = real::ONE; v <= m; v = v + real::ONE)
	{
		P0 = -(real::TWO * v - real::ONE) * pow(real::ONE - eta * eta, real::ONE / real::TWO) * P0;
	}
	P1 = (real::TWO * m + real::ONE) * eta * P0;
	if (remainder(n - m, real::TWO) == real::ZERO)
	{
		for (real r = real::ZERO; r <= n_dr - real::TWO; r = r + real::TWO)
		{
			if (r > real::ZERO)
			{
				P0 = (real::ONE / r) * (-(real::TWO * m + r - real::ONE) * P0 + (real::TWO * m + real::TWO * r - real::ONE) * eta * P1);
				P1 = (real::ONE / (r + real::ONE)) * (-(real::TWO * m + r) * P1 + (real::TWO * m + real::TWO * r + real::ONE) * eta * P0);
			}
			if (abs(eta) < real::ONE)
			{
				P0p = (real::ONE / (real::ONE - eta * eta)) * ((m + r + real::ONE) * eta * P0 - (r + real::ONE) * P1);
			}
			else
			{
				if (m > real::TWO)
				{
					P0p = real::ZERO;
				}
				else if (m > real::ONE)
				{
					P0p = -((m + r - real::ONE) * (m + r) * (m + r + real::ONE) * (m + r + real::TWO)) / real::FOUR;
				}
				else if (m > real::ZERO)
				{
					P0p = real::INF;
				}
				else
				{
					P0p = ((m + r) * (m + r + real::ONE)) / real::TWO;
				}
				if (eta == -real::ONE)
				{
					P0p = -P0p;
				}
			}
			change = dr[gzbi(r)] * P0;
			S1 = S1 + change;
			S1_adder.add(change);
			changep = dr[gzbi(r)] * P0p;
			S1p = S1p + changep;
			S1p_adder.add(changep);
			if (r > real::ZERO && abs(change) > real::ZERO && abs(change / S1) < real::SMALL_ENOUGH && abs(changep) > real::ZERO && abs(changep / S1p) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Smn1_1: " << (change / S1).get_string(10) << ", " << (changep / S1p).get_string(10) << ", " << ((S1 - S1_adder.calculate_sum()) / S1_adder.calculate_sum()).get_string(10) << ", " << ((S1p - S1p_adder.calculate_sum()) / S1p_adder.calculate_sum()).get_string(10) << ", " << (change / S1_adder.calculate_sum()).get_string(10) << ", " << (changep / S1p_adder.calculate_sum()).get_string(10) << std::endl;
				}
				break;
			}
		}
	}
	else
	{
		for (real r = real::ONE; r <= n_dr - real::ONE; r = r + real::TWO)
		{
			if (r > real::ONE)
			{
				P0 = (real::ONE / (r - real::ONE)) * (-(real::TWO * m + r - real::TWO) * P0 + (real::TWO * m + real::TWO * r - real::THREE) * eta * P1);
				P1 = (real::ONE / r) * (-(real::TWO * m + r - real::ONE) * P1 + (real::TWO * m + real::TWO * r - real::ONE) * eta * P0);
			}
			if (abs(eta) < real::ONE)
			{
				P1p = (real::ONE / (real::ONE - eta * eta)) * ((real::TWO * m + r) * P0 - (m + r) * eta * P1);
			}
			else
			{
				if (m > real::TWO)
				{
					P1p = real::ZERO;
				}
				else if (m > real::ONE)
				{
					P1p = -((m + r - real::ONE) * (m + r) * (m + r + real::ONE) * (m + r + real::TWO)) / real::FOUR;
				}
				else if (m > real::ZERO)
				{
					P1p = real::INF;
				}
				else
				{
					P1p = ((m + r) * (m + r + real::ONE)) / real::TWO;
				}
			}
			change = dr[gzbi(r)] * P1;
			S1 = S1 + change;
			S1_adder.add(change);
			changep = dr[gzbi(r)] * P1p;
			S1p = S1p + changep;
			S1p_adder.add(changep);
			if (r > real::ONE && abs(change) > real::ZERO && abs(change / S1) < real::SMALL_ENOUGH && abs(changep) > real::ZERO && abs(changep / S1p) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Smn1_1: " << (change / S1).get_string(10) << ", " << (changep / S1p).get_string(10) << ", " << ((S1 - S1_adder.calculate_sum()) / S1_adder.calculate_sum()).get_string(10) << ", " << ((S1p - S1p_adder.calculate_sum()) / S1p_adder.calculate_sum()).get_string(10) << ", " << (change / S1_adder.calculate_sum()).get_string(10) << ", " << (changep / S1p_adder.calculate_sum()).get_string(10) << std::endl;
				}
				break;
			}
		}
	}
	S1 = S1_adder.calculate_sum();
	S1p = S1p_adder.calculate_sum();
}

void calculate_Smn1_2(real & S1, real & S1p, bool verbose, const real & c, const real & m, const real & n, const real & n_c2k, const std::vector<real> & c2k, const real & eta)
{
	adder S1_adder;
	real a;
	real change;
	adder S1p_adder;
	real ap;
	real changep;
	
	S1 = real::ZERO;
	S1_adder.clear();
	for (real k = real::ZERO; k <= n_c2k - real::ONE; k = k + real::ONE)
	{
		if (k > real::ZERO)
		{
			a = a * (real::ONE - eta * eta);
		}
		else
		{
			a = real::ONE;
		}
		change = c2k[gzbi(k)] * a;
		S1 = S1 + change;
		S1_adder.add(change);
		if (k > real::ZERO && abs(change) > real::ZERO && abs(change / S1) < real::SMALL_ENOUGH)
		{
			if (verbose)
			{
				std::cout << "calculate_Smn1_2: " << (change / S1).get_string(10) << ", " << ((S1 - S1_adder.calculate_sum()) / S1_adder.calculate_sum()).get_string(10) << ", " << (change / S1_adder.calculate_sum()).get_string(10) << std::endl;
			}
			break;
		}
	}
	S1 = S1_adder.calculate_sum();
	S1p = real::ZERO;
	S1p_adder.clear();
	for (real k = real::ONE; k <= n_c2k - real::ONE; k = k + real::ONE)
	{
		if (k > real::ONE)
		{
			ap = ap * (real::ONE - eta * eta);
		}
		else
		{
			ap = real::ONE;
		}
		changep = c2k[gzbi(k)] * k * ap * (-real::TWO * eta);
		S1p = S1p + changep;
		S1p_adder.add(changep);
		if (k > real::ONE && abs(changep) > real::ZERO && abs(changep / S1p) < real::SMALL_ENOUGH)
		{
			if (verbose)
			{
				std::cout << "calculate_Smn1_2: " << (changep / S1p).get_string(10) << ", " << ((S1p - S1p_adder.calculate_sum()) / S1p_adder.calculate_sum()).get_string(10) << ", " << (changep / S1p_adder.calculate_sum()).get_string(10) << std::endl;
			}
			break;
		}
	}
	S1p = S1p_adder.calculate_sum();
	if (remainder(n - m, real::TWO) == real::ZERO)
	{
		if (m > real::ZERO)
		{
			S1p = pow(-real::ONE, m) * (m / real::TWO) * pow(real::ONE - eta * eta, m / real::TWO - real::ONE) * (-real::TWO * eta) * S1 + pow(-real::ONE, m) * pow(real::ONE - eta * eta, m / real::TWO) * S1p;
		}
		else
		{
			S1p = pow(-real::ONE, m) * S1p;
		}
		S1 = pow(-real::ONE, m) * pow(real::ONE - eta * eta, m / real::TWO) * S1;
	}
	else
	{
		if (m > real::ZERO)
		{
			S1p = pow(-real::ONE, m) * pow(real::ONE - eta * eta, m / real::TWO) * S1 + pow(-real::ONE, m) * eta * (m / real::TWO) * pow(real::ONE - eta * eta, m / real::TWO - real::ONE) * (-real::TWO * eta) * S1 + pow(-real::ONE, m) * eta * pow(real::ONE - eta * eta, m / real::TWO) * S1p;
		}
		else
		{
			S1p = pow(-real::ONE, m) * S1 + pow(-real::ONE, m) * eta * S1p;
		}
		S1 = pow(-real::ONE, m) * eta * pow(real::ONE - eta * eta, m / real::TWO) * S1;
	}
}

void calculate_Rmn1_1_shared(real & R1, real & R1p, bool verbose, const real & c, const real & m, const real & n, const real & n_dr, const std::vector<real> & dr, const real & xi)
{
	std::vector<real> jn;
	real b0;
	std::vector<real> a;
	std::vector<real> b;
	real N;
	real prev_N;
	real s;
	std::vector<real> jnp;
	adder R1_adder;
	adder R1p_adder;
	real d;
	real change;
	real changep;
	
	if (xi > real::ZERO)
	{
		jn.clear();
		for (real v = real::ZERO; v <= m + n_dr; v = v + real::ONE)
		{
			jn.push_back(real::ZERO);
		}
		jn[gzbi(m + n_dr)] = real::ONE;
		for (real v = m + n_dr; v >= real::ONE; v = v - real::ONE)
		{
			if (v < m + n_dr)
			{
				N = real::ONE / ((real::TWO * v + real::ONE) / (c * xi) - N);
			}
			else
			{
				b0 = real::ZERO;
				a.clear();
				b.clear();
				for (real i = v; i < v + real("10000.0"); i = i + real::ONE)
				{
					a.push_back(real::ONE);
					b.push_back((real::TWO * i + real::ONE) / (c * xi));
					N = -calculate_continued_fraction(b0, a, b);
					if (i > v)
					{
						if (abs((N - prev_N) / prev_N) < real::SMALL_ENOUGH)
						{
							if (verbose)
							{
								std::cout << "calculate_Rmn1_1: " << ((N - prev_N) / prev_N).get_string(10) << std::endl;
							}
							break;
						}
					}
					prev_N = N;
				}
			}
			jn[gzbi(v - real::ONE)] = jn[gzbi(v)] / N;
		}
		s = (sin(c * xi) / (c * xi)) / jn[gzbi(real::ZERO)];
		for (real v = real::ZERO; v <= m + n_dr; v = v + real::ONE)
		{
			jn[gzbi(v)] = s * jn[gzbi(v)];
		}
		jnp.clear();
		for (real v = real::ZERO; v <= m + n_dr - real::ONE; v = v + real::ONE)
		{
			jnp.push_back(real::ZERO);
		}
		for (real v = real::ZERO; v <= m + n_dr - real::ONE; v = v + real::ONE)
		{
			jnp[gzbi(v)] = (v / (c * xi)) * jn[gzbi(v)] - jn[gzbi(v + real::ONE)];
		}
	}
	else
	{
		jn.clear();
		jnp.clear();
		for (real v = real::ZERO; v <= m + n_dr - real::ONE; v = v + real::ONE)
		{
			jn.push_back(real::ZERO);
			jnp.push_back(real::ZERO);
		}
		jn[gzbi(real::ZERO)] = real::ONE;
		jnp[gzbi(real::ONE)] = real::ONE / real::THREE;
	}
	R1 = real::ZERO;
	R1_adder.clear();
	R1p = real::ZERO;
	R1p_adder.clear();
	if (remainder(n - m, real::TWO) == real::ZERO)
	{
		for (real r = real::ZERO; r <= n_dr - real::TWO; r = r + real::TWO)
		{
			if (r > real::ZERO)
			{
				d = d * -real::ONE * (((real::TWO * m + r - real::ONE) * (real::TWO * m + r)) / ((r - real::ONE) * r));
			}
			else
			{
				d = pow(-real::ONE, -(n - m) / real::TWO) * factorial(real::TWO * m);
			}
			change = d * dr[gzbi(r)] * jn[gzbi(m + r)];
			R1 = R1 + change;
			R1_adder.add(change);
			changep = d * dr[gzbi(r)] * jnp[gzbi(m + r)] * c;
			R1p = R1p + changep;
			R1p_adder.add(changep);
			if (r > real::ZERO && abs(change) > real::ZERO && abs(change / R1) < real::SMALL_ENOUGH && abs(changep) > real::ZERO && abs(changep / R1p) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Rmn1_1: " << (change / R1).get_string(10) << ", " << (changep / R1p).get_string(10) << ", " << ((R1 - R1_adder.calculate_sum()) / R1_adder.calculate_sum()).get_string(10) << ", " << ((R1p - R1p_adder.calculate_sum()) / R1p_adder.calculate_sum()).get_string(10) << ", " << (change / R1_adder.calculate_sum()).get_string(10) << ", " << (changep / R1p_adder.calculate_sum()).get_string(10) << std::endl;
				}
				break;
			}
		}
	}
	else
	{
		for (real r = real::ONE; r <= n_dr - real::ONE; r = r + real::TWO)
		{
			if (r > real::ONE)
			{
				d = d * -real::ONE * (((real::TWO * m + r - real::ONE) * (real::TWO * m + r)) / ((r - real::ONE) * r));
			}
			else
			{
				d = pow(-real::ONE, (real::ONE - (n - m)) / real::TWO) * factorial(real::TWO * m + real::ONE);
			}
			change = d * dr[gzbi(r)] * jn[gzbi(m + r)];
			R1 = R1 + change;
			R1_adder.add(change);
			changep = d * dr[gzbi(r)] * jnp[gzbi(m + r)] * c;
			R1p = R1p + changep;
			R1p_adder.add(changep);
			if (r > real::ONE && abs(change) > real::ZERO && abs(change / R1) < real::SMALL_ENOUGH && abs(changep) > real::ZERO && abs(changep / R1p) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Rmn1_1: " << (change / R1).get_string(10) << ", " << (changep / R1p).get_string(10) << ", " << ((R1 - R1_adder.calculate_sum()) / R1_adder.calculate_sum()).get_string(10) << ", " << ((R1p - R1p_adder.calculate_sum()) / R1p_adder.calculate_sum()).get_string(10) << ", " << (change / R1_adder.calculate_sum()).get_string(10) << ", " << (changep / R1p_adder.calculate_sum()).get_string(10) << std::endl;
				}
				break;
			}
		}
	}
	R1 = R1_adder.calculate_sum();
	R1p = R1p_adder.calculate_sum();
}

void calculate_Rmn2_1_shared(real & R2, real & R2p, bool verbose, const real & c, const real & m, const real & n, const real & n_dr, const std::vector<real> & dr, const real & xi)
{
	real y0;
	real y1;
	real y2;
	adder R2_adder;
	adder R2p_adder;
	real a;
	real y0p;
	real change;
	real changep;
	real y1p;
	
	y0 = -cos(c * xi) / (c * xi);
	y1 = -cos(c * xi) / ((c * xi) * (c * xi)) - sin(c * xi) / (c * xi);
	for (real v = real::ZERO; v <= m - real::ONE; v = v + real::ONE)
	{
		y2 = -y0 + ((real::TWO * v + real::THREE) / (c * xi)) * y1;
		y0 = y1;
		y1 = y2;
	}
	R2 = real::ZERO;
	R2_adder.clear();
	R2p = real::ZERO;
	R2p_adder.clear();
	if (remainder(n - m, real::TWO) == real::ZERO)
	{
		for (real r = real::ZERO; r <= n_dr - real::TWO; r = r + real::TWO)
		{
			if (r > real::ZERO)
			{
				a = a * -real::ONE * (((real::TWO * m + r - real::ONE) * (real::TWO * m + r)) / ((r - real::ONE) * r));
			}
			else
			{
				a = pow(-real::ONE, -(n - m) / real::TWO) * factorial(real::TWO * m);
			}
			y0p = ((m + r) / (c * xi)) * y0 - y1;
			change = a * dr[gzbi(r)] * y0;
			R2 = R2 + change;
			R2_adder.add(change);
			changep = a * dr[gzbi(r)] * y0p * c;
			R2p = R2p + changep;
			R2p_adder.add(changep);
			if (r > real::ZERO && abs(change) > real::ZERO && abs(change / R2) < real::SMALL_ENOUGH && abs(changep) > real::ZERO && abs(changep / R2p) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Rmn2_1: " << (change / R2).get_string(10) << ", " << (changep / R2p).get_string(10) << ", " << ((R2 - R2_adder.calculate_sum()) / R2_adder.calculate_sum()).get_string(10) << ", " << ((R2p - R2p_adder.calculate_sum()) / R2p_adder.calculate_sum()).get_string(10) << ", " << (change / R2_adder.calculate_sum()).get_string(10) << ", " << (changep / R2p_adder.calculate_sum()).get_string(10) << std::endl;
				}
				break;
			}
			y0 = -y0 + ((real::TWO * (m + r) + real::THREE) / (c * xi)) * y1;
			y1 = -y1 + ((real::TWO * (m + r + real::ONE) + real::THREE) / (c * xi)) * y0;
		}
	}
	else
	{
		for (real r = real::ONE; r <= n_dr - real::ONE; r = r + real::TWO)
		{
			if (r > real::ONE)
			{
				a = a * -real::ONE * (((real::TWO * m + r - real::ONE) * (real::TWO * m + r)) / ((r - real::ONE) * r));
			}
			else
			{
				a = pow(-real::ONE, (real::ONE - (n - m)) / real::TWO) * factorial(real::TWO * m + real::ONE);
			}
			y1p = y0 - ((m + r + real::ONE) / (c * xi)) * y1;
			change = a * dr[gzbi(r)] * y1;
			R2 = R2 + change;
			R2_adder.add(change);
			changep = a * dr[gzbi(r)] * y1p * c;
			R2p = R2p + changep;
			R2p_adder.add(changep);
			if (r > real::ONE && abs(change) > real::ZERO && abs(change / R2) < real::SMALL_ENOUGH && abs(changep) > real::ZERO && abs(changep / R2p) < real::SMALL_ENOUGH)
			{
				if (verbose)
				{
					std::cout << "calculate_Rmn2_1: " << (change / R2).get_string(10) << ", " << (changep / R2p).get_string(10) << ", " << ((R2 - R2_adder.calculate_sum()) / R2_adder.calculate_sum()).get_string(10) << ", " << ((R2p - R2p_adder.calculate_sum()) / R2p_adder.calculate_sum()).get_string(10) << ", " << (change / R2_adder.calculate_sum()).get_string(10) << ", " << (changep / R2p_adder.calculate_sum()).get_string(10) << std::endl;
				}
				break;
			}
			y0 = -y0 + ((real::TWO * (m + r - real::ONE) + real::THREE) / (c * xi)) * y1;
			y1 = -y1 + ((real::TWO * (m + r) + real::THREE) / (c * xi)) * y0;
		}
	}
	R2 = R2_adder.calculate_sum();
	R2p = R2p_adder.calculate_sum();
}
