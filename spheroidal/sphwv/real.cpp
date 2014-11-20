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

#include <cstdio>
#include <iostream>
#include <mpfr.h>
#include "real.hpp"
#include <string>

//
// max_reals is the maximum number of reals that can be active at the same
// time.  When there are too many reals in use, the program stops.  This
// serves as a way to prevent the program from using too much memory, and
// causing the machine to thrash or whatever.  n_reals is the number of reals
// being used, and is incremented or decremented in the real's constructor or
// destructor, respectively.
//
int real::max_reals = -1;
int real::n_reals = 0;

real real::NAN;
real real::INF;
real real::ZERO;
real real::ONE;
real real::TWO;
real real::THREE;
real real::FOUR;
real real::FIVE;
real real::PI;
real real::EPS;
real real::SMALL_ENOUGH;
complex complex::I;

void real::begin(int precision, int mmax_reals)
{
	max_reals = mmax_reals;
	n_reals = 0;
	mpfr_set_default_prec(precision);
	mpfr_clear(NAN.r);
	mpfr_init(NAN.r);
	mpfr_set_nan(NAN.r);
	mpfr_clear(INF.r);
	mpfr_init(INF.r);
	mpfr_set_inf(INF.r, 0);
	mpfr_clear(ZERO.r);
	mpfr_init(ZERO.r);
	ZERO = real("0.0");
	mpfr_clear(ONE.r);
	mpfr_init(ONE.r);
	ONE = real("1.0");
	mpfr_clear(TWO.r);
	mpfr_init(TWO.r);
	TWO = real("2.0");
	mpfr_clear(THREE.r);
	mpfr_init(THREE.r);
	THREE = real("3.0");
	mpfr_clear(FOUR.r);
	mpfr_init(FOUR.r);
	FOUR = real("4.0");
	mpfr_clear(FIVE.r);
	mpfr_init(FIVE.r);
	FIVE = real("5.0");
	mpfr_clear(PI.r);
	mpfr_init(PI.r);
	mpfr_const_pi(PI.r, MPFR_RNDN);
	mpfr_clear(EPS.r);
	mpfr_init(EPS.r);
	EPS = ONE;
	while (ONE + EPS != ONE)
	{
		EPS = EPS / TWO;
	}
	mpfr_clear(SMALL_ENOUGH.r);
	mpfr_init(SMALL_ENOUGH.r);
	SMALL_ENOUGH = EPS;
}

real::real()
{
	mpfr_init(r);
	mpfr_set_str(r, "0.0", 10, MPFR_RNDN);
	++n_reals;
	if (max_reals != -1 && n_reals > max_reals)
	{
		std::cout << "error: program has reached memory capacity" << std::endl;
		exit(1);
	}
}

real::real(const real & a)
{
	mpfr_init(r);
	mpfr_set(r, a.r, MPFR_RNDN);
	++n_reals;
	if (max_reals != -1 && n_reals > max_reals)
	{
		std::cout << "error: program has reached memory capacity" << std::endl;
		exit(1);
	}
}

real::real(int a)
{
	mpfr_init(r);
	mpfr_set_si(r, a, MPFR_RNDZ);
	++n_reals;
	if (max_reals != -1 && n_reals > max_reals)
	{
		std::cout << "error: program has reached memory capacity" << std::endl;
		exit(1);
	}
}

real::real(double a)
{
	mpfr_init(r);
	mpfr_set_d(r, a, MPFR_RNDN);
	++n_reals;
	if (max_reals != -1 && n_reals > max_reals)
	{
		std::cout << "error: program has reached memory capacity" << std::endl;
		exit(1);
	}
}

real::real(const std::string & a)
{
	mpfr_init(r);
	mpfr_set_str(r, a.c_str(), 10, MPFR_RNDN);
	++n_reals;
	if (max_reals != -1 && n_reals > max_reals)
	{
		std::cout << "error: program has reached memory capacity" << std::endl;
		exit(1);
	}
}

real::~real()
{
	mpfr_clear(r);
	--n_reals;
}

real & real::operator =(const real & a)
{
	mpfr_set(r, a.r, MPFR_RNDN);
	return *this;
}

real & real::operator =(int a)
{
	mpfr_set_si(r, a, MPFR_RNDN);
	return *this;
}

real & real::operator =(double a)
{
	mpfr_set_d(r, a, MPFR_RNDN);
	return *this;
}

real & real::operator =(const std::string & a)
{
	mpfr_set_str(r, a.c_str(), 10, MPFR_RNDN);
	return *this;
}

int real::get_int() const
{
	return mpfr_get_si(r, MPFR_RNDZ);
}

double real::get_double() const
{
	return mpfr_get_d(r, MPFR_RNDN);
}

std::string real::get_string(int p) const
{
	std::string string;
	char *raw_string;
	mpfr_exp_t exp;
	
	if (mpfr_nan_p(r) != 0)
	{
		string = "nan";
	}
	else if (mpfr_inf_p(r) != 0)
	{
		if (mpfr_sgn(r) < 0)
		{
			string = "-inf";
		}
		else
		{
			string = "inf";
		}
	}
	else
	{
		// Dynamically allocate a character array to hold the base-10
		// representation of the real.
		raw_string = new char[(int)mpfr_get_default_prec()];
		mpfr_get_str(raw_string, &exp, 10, p, r, MPFR_RNDN);
		string = std::string(raw_string);
		if (mpfr_zero_p(r) == 0)
		{
			// The cast is to prevent a warning in case mpfr_exp_t is not
			// defined as an int, but as a short or a long.
			std::sprintf(raw_string, "%d", (int)(exp - 1));
		}
		else
		{
			std::sprintf(raw_string, "0");
		}
		if (string[0] == '-')
		{
			string = string.substr(0, 2) + std::string(".") + string.substr(2) + std::string("e") + std::string(raw_string);
		}
		else
		{
			string = string.substr(0, 1) + std::string(".") + string.substr(1) + std::string("e") + std::string(raw_string);
		}
		// Free the array.
		delete[] raw_string;
	}
	return string;
}

real operator +(const real & a, const real & b)
{
	real x;
	
	mpfr_add(x.r, a.r, b.r, MPFR_RNDN);
	return x;
}

real operator -(const real & a, const real & b)
{
	real x;
	
	mpfr_sub(x.r, a.r, b.r, MPFR_RNDN);
	return x;
}

real operator *(const real & a, const real & b)
{
	real x;
	
	mpfr_mul(x.r, a.r, b.r, MPFR_RNDN);
	return x;
}

real operator /(const real & a, const real & b)
{
	real x;
	
	mpfr_div(x.r, a.r, b.r, MPFR_RNDN);
	return x;
}

real operator -(const real & a)
{
	real x;
	
	mpfr_neg(x.r, a.r, MPFR_RNDN);
	return x;
}

bool operator >(const real & a, const real & b)
{
	return mpfr_greater_p(a.r, b.r) != 0;
}

bool operator >=(const real & a, const real & b)
{
	return mpfr_greaterequal_p(a.r, b.r) != 0;
}

bool operator <(const real & a, const real & b)
{
	return mpfr_less_p(a.r, b.r) != 0;
}

bool operator <=(const real & a, const real & b)
{
	return mpfr_lessequal_p(a.r, b.r) != 0;
}

bool operator ==(const real & a, const real & b)
{
	return mpfr_equal_p(a.r, b.r) != 0;
}

bool operator !=(const real & a, const real & b)
{
	return mpfr_equal_p(a.r, b.r) == 0;
}

real abs(const real & a)
{
	real x;
	
	mpfr_abs(x.r, a.r, MPFR_RNDN);
	return x;
}

real atan(const real & a)
{
	real x;
	
	mpfr_atan(x.r, a.r, MPFR_RNDN);
	return x;
}

real atan2(const real & b, const real & a)
{
	real x;
	
	mpfr_atan2(x.r, b.r, a.r, MPFR_RNDN);
	return x;
}

real cos(const real & a)
{
	real x;
	
	mpfr_cos(x.r, a.r, MPFR_RNDN);
	return x;
}

real factorial(const real & a)
{
	real x;
	
	mpfr_fac_ui(x.r, mpfr_get_ui(a.r, MPFR_RNDZ), MPFR_RNDN);
	return x;
}

real log(const real & a)
{
	real x;
	
	mpfr_log(x.r, a.r, MPFR_RNDN);
	return x;
}

real max(const real & a, const real & b)
{
	real x;
	
	mpfr_max(x.r, a.r, b.r, MPFR_RNDN);
	return x;
}

real pochhammer(const real & a, const real & k)
{
	real x;
	
	x = real::ONE;
	for (real i = real::ZERO; i <= k - real::ONE; i = i + real::ONE)
	{
		x = x * (a + i);
	}
	return x;
}

real pow(const real & a, const real & b)
{
	real x;
	
	mpfr_pow(x.r, a.r, b.r, MPFR_RNDN);
	return x;
}

real remainder(const real & a, const real & b)
{
	real x;
	
	mpfr_remainder(x.r, a.r, b.r, MPFR_RNDN);
	return x;
}

real round(const real & a)
{
	real x;
	
	mpfr_rint(x.r, a.r, MPFR_RNDN);
	return x;
}

real sin(const real & a)
{
	real x;
	
	mpfr_sin(x.r, a.r, MPFR_RNDN);
	return x;
}

int gzbi(const real & i)
{
	return i.get_int();
}

int gnobi(const real & i)
{
	return -i.get_int() - 1;
}

void complex::begin()
{
	mpfr_clear(I.a.r);
	mpfr_init(I.a.r);
	mpfr_clear(I.b.r);
	mpfr_init(I.b.r);
	I = complex(real::ZERO, real::ONE);
}

complex::complex()
{
	a = real::ZERO;
	b = real::ZERO;
}

complex::complex(const real & aa)
{
	a = aa;
	b = real::ZERO;
}

complex::complex(const real & aa, const real & bb)
{
	a = aa;
	b = bb;
}

complex & complex::operator =(const real & aa)
{
	a = aa;
	b = real::ZERO;
	return *this;
}

complex & complex::operator =(const complex & aa)
{
	a = aa.a;
	b = aa.b;
	return *this;
}

std::string complex::get_string(int p) const
{
	return a.get_string(p) + " + i * " + b.get_string(p);
}

complex operator +(const complex & a, const complex & b)
{
	complex x;
	
	x.a = a.a + b.a;
	x.b = a.b + b.b;
	return x;
}

complex operator -(const complex & a, const complex & b)
{
	complex x;
	
	x.a = a.a - b.a;
	x.b = a.b - b.b;
	return x;
}

complex operator *(const complex & a, const complex & b)
{
	complex x;
	
	x.a = a.a * b.a - a.b * b.b;
	x.b = a.a * b.b + a.b * b.a;
	return x;
}

complex operator /(const complex & a, const complex & b)
{
	real r;
	complex x;
	
	r = b.a * b.a + b.b * b.b;
	x.a = (a.a * b.a + a.b * b.b) / r;
	x.b = (a.b * b.a - a.a * b.b) / r;
	return x;
}

complex operator -(const complex & a)
{
	complex x;
	
	x.a = -a.a;
	x.b = -a.b;
	return x;
}

real abs(const complex & a)
{
	return pow(a.a * a.a + a.b * a.b, real::ONE / real::TWO);
}

//
// Returns the value on the principal branch.
//
complex log(const complex & a)
{
	real r;
	real angle;
	complex x;
	
	r = abs(a);
	angle = atan2(a.b, a.a);
	x = log(r) + complex::I * angle;
	return x;
}

//
// Returns the value on the principal branch.
//
complex pow(const complex & a, const real & b)
{
	real r;
	real angle;
	complex x;
	
	r = abs(a);
	angle = atan2(a.b, a.a);
	x = pow(r, b) * (cos(b * angle) + complex::I * sin(b * angle));
	return x;
}
