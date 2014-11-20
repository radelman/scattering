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

#ifndef COMMON_MAIN_HPP
#define COMMON_MAIN_HPP

#include "real.hpp"
#include <string>
#include <vector>

std::string generate_name(const real & c, const real & m, const real & n, const std::string & name);
bool save_lambdamn(bool verbose, const real & c, const real & m, const real & n);
bool open_lambdamn(real & lambda, const real & c, const real & m, const real & n);
bool save_drmn(bool verbose, const real & c, const real & m, const real & n, real & n_dr, const real & dr_min);
bool open_drmn(real & n_dr, std::vector<real> & dr, const real & c, const real & m, const real & n);
bool save_drmn_neg(bool verbose, const real & c, const real & m, const real & n, real & n_dr_neg, const real & dr_neg_min);
bool open_drmn_neg(real & n_dr_neg, std::vector<real> & dr_neg, const real & c, const real & m, const real & n);
bool save_Nmn(bool verbose, const real & c, const real & m, const real & n);
bool open_Nmn(real & N, const real & c, const real & m, const real & n);
bool save_Fmn(bool verbose, const real & c, const real & m, const real & n);
bool open_Fmn(real & N, const real & c, const real & m, const real & n);
bool save_kmn1(bool verbose, const real & c, const real & m, const real & n);
bool open_kmn1(real & k1, const real & c, const real & m, const real & n);
bool save_kmn2(bool verbose, const real & c, const real & m, const real & n);
bool open_kmn2(real & k2, const real & c, const real & m, const real & n);
bool save_c2kmn(bool verbose, const real & c, const real & m, const real & n, real & n_c2k, const real & c2k_min);
bool open_c2kmn(real & n_c2k, std::vector<real> & c2k, const real & c, const real & m, const real & n);
bool save_Smn1(bool verbose, const real & c, const real & m, const real & n, const real & a, const real & b, const real & d, const std::string & arg_type, int p);
int parse_args(int argc, char **argv);

#endif
