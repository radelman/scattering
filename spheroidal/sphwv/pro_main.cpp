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

#include "common_main.hpp"
#include "common_spheroidal.hpp"
#include <cstdio>
#include <iostream>
#include "pro_spheroidal.hpp"
#include "real.hpp"
#include <string>
#include <vector>

static bool save_Rmn(bool verbose, const real & c, const real & m, const real & n, const real & a, const real & b, const real & d, const std::string & arg_type, const std::string & which, int p);

std::string generate_name(const real & c, const real & m, const real & n, const std::string & name)
{
	char raw_string[4096];
	
	std::sprintf(raw_string, "data/pro_%08d_%03d_%03d_%s.txt", (real("1000.0") * c).get_int(), m.get_int(), n.get_int(), name.c_str());
	return std::string(raw_string);
}

static bool save_Rmn(bool verbose, const real & c, const real & m, const real & n, const real & a, const real & b, const real & d, const std::string & arg_type, const std::string & which, int p)
{
	real n_dr;
	std::vector<real> dr;
	real n_dr_neg;
	std::vector<real> dr_neg;
	real F;
	real k1;
	real k2;
	real n_c2k;
	std::vector<real> c2k;
	real xi;
	real log_xi;
	real R1_1;
	real R1p_1;
	real R1_2;
	real R1p_2;
	real R1_log_abs_difference;
	real R1p_log_abs_difference;
	real R2_1;
	real R2p_1;
	real R2_2;
	real R2p_2;
	real R2_log_abs_difference;
	real R2p_log_abs_difference;
	real W;
	real log_W;
	real W_1_1_log_abs_error;
	real W_1_2_log_abs_error;
	real W_2_1_log_abs_error;
	real W_2_2_log_abs_error;
	
	if (!open_drmn(n_dr, dr, c, m, n) ||
	    !open_drmn_neg(n_dr_neg, dr_neg, c, m, n) ||
	    !open_Fmn(F, c, m, n) ||
	    !open_kmn1(k1, c, m, n) ||
	    !open_kmn2(k2, c, m, n) ||
	    !open_c2kmn(n_c2k, c2k, c, m, n))
	{
		return false;
	}
	for (real i = a; i <= b; i = i + d)
	{
		if (arg_type == "xi")
		{
			xi = i;
		}
		else
		{
			xi = pow(i * i + real::ONE, real::ONE / real::TWO);
		}
		log_xi = log(xi - real::ONE);
		if (which.find("R1_1") != std::string::npos)
		{
			calculate_Rmn1_1(R1_1, R1p_1, verbose, c, m, n, n_dr, dr, F, xi);
		}
		else
		{
			R1_1 = real::NAN;
			R1p_1 = real::NAN;
		}
		if (which.find("R1_2") != std::string::npos)
		{
			calculate_Rmn1_2(R1_2, R1p_2, verbose, c, m, n, k1, n_c2k, c2k, xi);
		}
		else
		{
			R1_2 = real::NAN;
			R1p_2 = real::NAN;
		}
		R1_log_abs_difference = log(abs(R1_1 - R1_2));
		R1p_log_abs_difference = log(abs(R1p_1 - R1p_2));
		if (which.find("R2_1") != std::string::npos)
		{
			calculate_Rmn2_1(R2_1, R2p_1, verbose, c, m, n, n_dr, dr, F, xi);
		}
		else
		{
			R2_1 = real::NAN;
			R2p_1 = real::NAN;
		}
		if (which.find("R2_2") != std::string::npos)
		{
			calculate_Rmn2_2(R2_2, R2p_2, verbose, c, m, n, n_dr, dr, n_dr_neg, dr_neg, k2, xi);
		}
		else
		{
			R2_2 = real::NAN;
			R2p_2 = real::NAN;
		}
		R2_log_abs_difference = log(abs(R2_1 - R2_2));
		R2p_log_abs_difference = log(abs(R2p_1 - R2p_2));
		W = real::ONE / (c * (xi * xi - real::ONE));
		log_W = log(W);
		W_1_1_log_abs_error = log(abs(R1_1 * R2p_1 - R1p_1 * R2_1 - W));
		W_1_2_log_abs_error = log(abs(R1_1 * R2p_2 - R1p_1 * R2_2 - W));
		W_2_1_log_abs_error = log(abs(R1_2 * R2p_1 - R1p_2 * R2_1 - W));
		W_2_2_log_abs_error = log(abs(R1_2 * R2p_2 - R1p_2 * R2_2 - W));
		std::cout << i.get_string(p) << ","
		          << xi.get_string(p) << ","
		          << log_xi.get_string(p) << ","
		          << R1_1.get_string(p) << "," << R1p_1.get_string(p) << ","
		          << R1_2.get_string(p) << "," << R1p_2.get_string(p) << ","
		          << R1_log_abs_difference.get_string(p) << "," << R1p_log_abs_difference.get_string(p) << ","
		          << R2_1.get_string(p) << "," << R2p_1.get_string(p) << ","
		          << R2_2.get_string(p) << "," << R2p_2.get_string(p) << ","
		          << R2_log_abs_difference.get_string(p) << "," << R2p_log_abs_difference.get_string(p) << ","
		          << W.get_string(p) << ","
		          << log_W.get_string(p) << ","
		          << W_1_1_log_abs_error.get_string(p) << ","
		          << W_1_2_log_abs_error.get_string(p) << ","
		          << W_2_1_log_abs_error.get_string(p) << ","
		          << W_2_2_log_abs_error.get_string(p) << std::endl;
	}
	return true;
}

int parse_args(int argc, char **argv)
{
	std::string argument;
	std::string value;
	bool verbose;
	bool verbose_entered;
	real c;
	bool c_entered;
	real m;
	bool m_entered;
	real n;
	bool n_entered;
	std::string w;
	bool w_entered;
	real n_dr;
	bool n_dr_entered;
	real dr_min;
	bool dr_min_entered;
	real n_dr_neg;
	bool n_dr_neg_entered;
	real dr_neg_min;
	bool dr_neg_min_entered;
	real n_c2k;
	bool n_c2k_entered;
	real c2k_min;
	bool c2k_min_entered;
	real a;
	bool a_entered;
	real b;
	bool b_entered;
	real d;
	bool d_entered;
	std::string arg_type;
	bool arg_type_entered;
	std::string which;
	bool which_entered;
	int p;
	bool p_entered;

	verbose_entered = false;
	c_entered = false;
	m_entered = false;
	n_entered = false;
	w_entered = false;
	n_dr_entered = false;
	dr_min_entered = false;
	n_dr_neg_entered = false;
	dr_neg_min_entered = false;
	n_c2k_entered = false;
	c2k_min_entered = false;
	a_entered = false;
	b_entered = false;
	d_entered = false;
	arg_type_entered = false;
	which_entered = false;
	p_entered = false;
	for (int i = 1; i < argc; i = i + 2)
	{
		argument = std::string(argv[i]);
		value = std::string(argv[i + 1]);
		if (argument == "-verbose")
		{
			verbose = value == "y";
			verbose_entered = true;
		}
		else if (argument == "-c")
		{
			c = real(value);
			c_entered = true;
		}
		else if (argument == "-m")
		{
			m = real(value);
			m_entered = true;
		}
		else if (argument == "-n")
		{
			n = real(value);
			n_entered = true;
		}
		else if (argument == "-w")
		{
			w = value;
			w_entered = true;
		}
		else if (argument == "-n_dr")
		{
			n_dr = real(value);
			n_dr_entered = true;
		}
		else if (argument == "-dr_min")
		{
			dr_min = real(value);
			dr_min_entered = true;
		}
		else if (argument == "-n_dr_neg")
		{
			n_dr_neg = real(value);
			n_dr_neg_entered = true;
		}
		else if (argument == "-dr_neg_min")
		{
			dr_neg_min = real(value);
			dr_neg_min_entered = true;
		}
		else if (argument == "-n_c2k")
		{
			n_c2k = real(value);
			n_c2k_entered = true;
		}
		else if (argument == "-c2k_min")
		{
			c2k_min = real(value);
			c2k_min_entered = true;
		}
		else if (argument == "-a")
		{
			if (value.substr(0, 4) == "1+2^")
			{
				a = real::ONE + pow(real::TWO, real(value.substr(4, value.length() - 4)));
			}
			else if (value.substr(0, 2) == "2^")
			{
				a = pow(real::TWO, real(value.substr(2, value.length() - 2)));
			}
			else
			{
				a = real(value);
			}
			a_entered = true;
		}
		else if (argument == "-b")
		{
			if (value.substr(0, 4) == "1+2^")
			{
				b = real::ONE + pow(real::TWO, real(value.substr(4, value.length() - 4)));
			}
			else if (value.substr(0, 2) == "2^")
			{
				b = pow(real::TWO, real(value.substr(2, value.length() - 2)));
			}
			else
			{
				b = real(value);
			}
			b_entered = true;
		}
		else if (argument == "-d")
		{
			d = real(value);
			d_entered = true;
		}
		else if (argument == "-arg_type")
		{
			arg_type = value;
			arg_type_entered = true;
		}
		else if (argument == "-which")
		{
			which = value;
			which_entered = true;
		}
		else if (argument == "-p")
		{
			p = std::atoi(value.c_str());
			p_entered = true;
		}
	}
	if (!verbose_entered)
	{
		std::cout << "no value of verbose was entered..." << std::endl;
		return 1;
	}
	else if (!c_entered || !m_entered || !n_entered || !w_entered)
	{
		std::cout << "no value of c, m, n, and/or w was entered..." << std::endl;
		return 1;
	}
	else if (w == "lambda")
	{
		save_lambdamn(verbose, c, m, n);
	}
	else if (w == "dr")
	{
		if (!n_dr_entered || !dr_min_entered)
		{
			std::cout << "no value of n_dr and/or dr_min was entered..." << std::endl;
			return 1;
		}
		save_drmn(verbose, c, m, n, n_dr, dr_min);
	}
	else if (w == "dr_neg")
	{
		if (!n_dr_neg_entered || !dr_neg_min_entered)
		{
			std::cout << "no value of n_dr_neg and/or dr_neg_min was entered..." << std::endl;
			return 1;
		}
		save_drmn_neg(verbose, c, m, n, n_dr_neg, dr_neg_min);
	}
	else if (w == "N")
	{
		save_Nmn(verbose, c, m, n);
	}
	else if (w == "F")
	{
		save_Fmn(verbose, c, m, n);
	}
	else if (w == "k1")
	{
		save_kmn1(verbose, c, m, n);
	}
	else if (w == "k2")
	{
		save_kmn2(verbose, c, m, n);
	}
	else if (w == "c2k")
	{
		if (!n_c2k_entered || !c2k_min_entered)
		{
			std::cout << "no value of n_c2k and/or c2k_min was entered..." << std::endl;
			return 1;
		}
		save_c2kmn(verbose, c, m, n, n_c2k, c2k_min);
	}
	else if (w == "everything")
	{
		if (!n_dr_entered || !dr_min_entered || !n_dr_neg_entered || !dr_neg_min_entered || !n_c2k_entered || !c2k_min_entered)
		{
			std::cout << "no value of n_dr, dr_min, n_dr_neg, dr_neg_min, n_c2k and/or c2k_min was entered..." << std::endl;
			return 1;
		}
		save_lambdamn(verbose, c, m, n);
		save_drmn(verbose, c, m, n, n_dr, dr_min);
		save_drmn_neg(verbose, c, m, n, n_dr_neg, dr_neg_min);
		save_Nmn(verbose, c, m, n);
		save_Fmn(verbose, c, m, n);
		save_kmn1(verbose, c, m, n);
		save_kmn2(verbose, c, m, n);
		save_c2kmn(verbose, c, m, n, n_c2k, c2k_min);
	}
	else if (!a_entered || !b_entered || !d_entered || !arg_type_entered || !p_entered)
	{
		std::cout << "no value of a, b, d, arg_type, and/or p was entered..." << std::endl;
		return 1;
	}
	else if (w == "S1")
	{
		save_Smn1(verbose, c, m, n, a, b, d, arg_type, p);
	}
	else if (!which_entered)
	{
		std::cout << "no value of which was entered..." << std::endl;
		return 1;
	}
	else if (w == "R")
	{
		save_Rmn(verbose, c, m, n, a, b, d, arg_type, which, p);
	}
	return 0;
}
