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

#ifndef OBL_SPHEROIDAL_HPP
#define OBL_SPHEROIDAL_HPP

#include "real.hpp"
#include <vector>

real calculate_Qmn(bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr, const real & k1, real & n_c2k, std::vector<real> & c2k);
void calculate_B2rmn(std::vector<real> & B2r, bool verbose, const real & c, const real & m, const real & n, const real & lambda, real & n_dr, std::vector<real> & dr, const real & k1, real & n_c2k, std::vector<real> & c2k, const real & Q, real & n_B2r, const real & B2r_min);
void calculate_Rmn1_1(real & R1, real & R1p, bool verbose, const real & c, const real & m, const real & n, const real & n_dr, const std::vector<real> & dr, const real & F, const real & xi);
void calculate_Rmn1_2(real & R1, real & R1p, bool verbose, const real & c, const real & m, const real & n, const real & k1, const real & n_c2k, const std::vector<real> & c2k, const real & xi);
void calculate_Rmn2_1(real & R2, real & R2p, bool verbose, const real & c, const real & m, const real & n, const real & n_dr, const std::vector<real> & dr, const real & F, const real & xi);
void calculate_Rmn2_2(real & R2, real & R2p, bool verbose, const real & c, const real & m, const real & n, const real & n_dr, const std::vector<real> & dr, const real & n_dr_neg, const std::vector<real> & dr_neg, const real & k2, const real & xi);
void calculate_Rmn2_3(real & R2, real & R2p, bool verbose, const real & c, const real & m, const real & n, const real & Q, const real & n_B2r, const std::vector<real> & B2r, const real & xi, const real & R1, const real & R1p);

#endif
