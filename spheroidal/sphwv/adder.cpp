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
#include "real.hpp"
#include <vector>

static bool calculate_sum_comparator(const real & a, const real & b);

adder::adder()
{
	clear();
}

void adder::clear()
{
	addends.clear();
	prev_sum = real::NAN;
}

void adder::add(const real & a)
{
	addends.push_back(a);
	prev_sum = real::NAN;
}

static bool calculate_sum_comparator(const real & a, const real & b)
{
	return abs(a) < abs(b);
}

//
// Calculate the sum by sorting the vector of addends by magnitude in
// ascending order, and then using pairwise addition.
//
real adder::calculate_sum()
{
	std::vector<int> nums;
	std::vector<real> sums;
	real sum;
	
	if (prev_sum == prev_sum)
	{
		sum = prev_sum;
	}
	else
	{
		std::sort(addends.begin(), addends.end(), calculate_sum_comparator);
		nums.clear();
		sums.clear();
		for (int i = 0; i < (int)addends.size(); ++i)
		{
			nums.push_back(1);
			sums.push_back(addends[i]);
			while ((int)nums.size() > 1 && nums[(int)nums.size() - 2] == nums[(int)nums.size() - 1])
			{
				nums[(int)nums.size() - 2] = nums[(int)nums.size() - 2] + nums[(int)nums.size() - 1];
				sums[(int)nums.size() - 2] = sums[(int)nums.size() - 2] + sums[(int)nums.size() - 1];
				nums.pop_back();
				sums.pop_back();
			}
		}
		std::sort(sums.begin(), sums.end(), calculate_sum_comparator);
		sum = real::ZERO;
		for (int i = 0; i < (int)sums.size(); ++i)
		{
			sum = sum + sums[i];
		}
	}
	prev_sum = sum;
	return sum;
}

complex_adder::complex_adder()
{
	real_adder.clear();
	imag_adder.clear();
}

void complex_adder::clear()
{
	real_adder.clear();
	imag_adder.clear();
}

void complex_adder::add(const complex & a)
{
	real_adder.add(a.a);
	imag_adder.add(a.b);
}

complex complex_adder::calculate_sum()
{
	return real_adder.calculate_sum() + complex::I * imag_adder.calculate_sum();
}
