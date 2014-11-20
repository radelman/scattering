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

#ifndef COMMON_SPHEROIDAL_HPP
#define COMMON_SPHEROIDAL_HPP

#include "real.hpp"
#include <vector>

//
// A lot of the code in pro_sphwv and obl_sphwv that calculates the
// coefficients are the exact same.  The only difference is whether c ^ 2 or
// -c ^ 2 is used.  For that, there's a function called calculate_c_squared,
// which is defined differently in pro_ and obl_spheroidal.cpp.
//
real calculate_c_squared(const real & c);

real calculate_continued_fraction(const real & b0, const std::vector<real> & a, const std::vector<real> & b);
void calculate_lambdamn(real & lambda, bool verbose, const real & c, const real & m, const real & n, const real & lambda_approx);
void calculate_drmn(std::vector<real> & dr, bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, const real & dr_min);
void calculate_drmn_neg(std::vector<real> & dr_neg, bool verbose, const real & c, const real & m, const real & n, const real & lambda, const real & n_dr, const std::vector<real> & dr, real & n_dr_neg, const real & dr_neg_min);
real calculate_Nmn(bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr);
real calculate_Fmn(bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr);
real calculate_kmn1(bool verbose, const real & c, const real & m, const real & n, const real & n_dr, const std::vector<real> & dr, const real & F);
real calculate_kmn2(bool verbose, const real & c, const real & m, const real & n, const real & lambda, const real & n_dr, const std::vector<real> & dr, real & n_dr_neg, std::vector<real> & dr_neg, const real & F);
void calculate_c2kmn(std::vector<real> & c2k, bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr, real & n_c2k, const real & c2k_min);
void calculate_Smn1_1(real & S1, real & S1p, bool verbose, const real & c, const real & m, const real & n, const real & n_dr, const std::vector<real> & dr, const real & eta);
void calculate_Smn1_2(real & S1, real & S1p, bool verbose, const real & c, const real & m, const real & n, const real & n_c2k, const std::vector<real> & c2k, const real & eta);
void calculate_Rmn1_1_shared(real & R1, real & R1p, bool verbose, const real & c, const real & m, const real & n, const real & n_dr, const std::vector<real> & dr, const real & xi);
void calculate_Rmn2_1_shared(real & R2, real & R2p, bool verbose, const real & c, const real & m, const real & n, const real & n_dr, const std::vector<real> & dr, const real & xi);

#endif
