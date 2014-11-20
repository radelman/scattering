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

#ifndef REAL_HPP
#define REAL_HPP

#include <mpfr.h>
#include <string>

class real
{
public:
	static int max_reals;
	static int n_reals;
	static real NAN;
	static real INF;
	static real ZERO;
	static real ONE;
	static real TWO;
	static real THREE;
	static real FOUR;
	static real FIVE;
	static real PI;
	static real EPS;
	static real SMALL_ENOUGH;
	
	static void begin(int precision, int mmax_reals);
	
	mpfr_t r;
	
	real();
	real(const real & a);
	real(int a);
	real(double a);
	real(const std::string & a);
	~real();
	real & operator =(const real & a);
	real & operator =(int a);
	real & operator =(double a);
	real & operator =(const std::string & a);
	int get_int() const;
	double get_double() const;
	std::string get_string(int p) const;
};

class complex
{
public:
	static complex I;
	
	static void begin();
	
	real a;
	real b;
	
	complex();
	complex(const real & aa);
	complex(const real & aa, const real & bb);
	complex & operator =(const real & aa);
	complex & operator =(const complex & aa);
	std::string get_string(int p) const;
};

real operator +(const real & a, const real & b);
real operator -(const real & a, const real & b);
real operator *(const real & a, const real & b);
real operator /(const real & a, const real & b);
real operator -(const real & a);
bool operator >(const real & a, const real & b);
bool operator >=(const real & a, const real & b);
bool operator <(const real & a, const real & b);
bool operator <=(const real & a, const real & b);
bool operator ==(const real & a, const real & b);
bool operator !=(const real & a, const real & b);
real abs(const real & a);
real atan(const real & a);
real atan2(const real & b, const real & a);
real cos(const real & a);
real factorial(const real & a);
real log(const real & a);
real max(const real & a, const real & b);
real pochhammer(const real & a, const real & k);
real pow(const real & a, const real & b);
real remainder(const real & a, const real & b);
real round(const real & a);
real sin(const real & a);
int gzbi(const real & i);
int gnobi(const real & i);
complex operator +(const complex & a, const complex & b);
complex operator -(const complex & a, const complex & b);
complex operator *(const complex & a, const complex & b);
complex operator /(const complex & a, const complex & b);
complex operator -(const complex & a);
real abs(const complex & a);
complex log(const complex & a);
complex pow(const complex & a, const real & b);

#endif
