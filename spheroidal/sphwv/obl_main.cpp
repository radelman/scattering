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
#include "io.hpp"
#include <iostream>
#include "obl_spheroidal.hpp"
#include "real.hpp"
#include <string>
#include <vector>

static bool save_Qmn(bool verbose, const real & c, const real & m, const real & n);
static bool open_Qmn(real & Q, const real & c, const real & m, const real & n);
static bool save_B2rmn(bool verbose, const real & c, const real & m, const real & n, real & n_B2r, const real & B2r_min);
static bool open_B2rmn(real & n_B2r, std::vector<real> & B2r, const real & c, const real & m, const real & n);
static bool save_Rmn(bool verbose, const real & c, const real & m, const real & n, const real & a, const real & b, const real & d, const std::string & arg_type, const std::string & which, int p);

std::string generate_name(const real & c, const real & m, const real & n, const std::string & name)
{
	char raw_string[4096];
	
	std::sprintf(raw_string, "data/obl_%08d_%03d_%03d_%s.txt", (real("1000.0") * c).get_int(), m.get_int(), n.get_int(), name.c_str());
	return std::string(raw_string);
}

static bool save_Qmn(bool verbose, const real & c, const real & m, const real & n)
{
	real lambda;
	real n_dr;
	std::vector<real> dr;
	real k1;
	real n_c2k;
	std::vector<real> c2k;
	real Q;
	
	if (!open_lambdamn(lambda, c, m, n) ||
	    !open_drmn(n_dr, dr, c, m, n) ||
	    !open_kmn1(k1, c, m, n) ||
	    !open_c2kmn(n_c2k, c2k, c, m, n))
	{
		return false;
	}
	Q = calculate_Qmn(verbose, c, m, n, lambda, n_dr, dr, k1, n_c2k, c2k);
	if (!save_data(generate_name(c, m, n, "Q"), Q))
	{
		std::cout << "can't save Q..." << std::endl;
		return false;
	}
	if (!save_log_abs_data(generate_name(c, m, n, "log_abs_Q"), Q))
	{
		std::cout << "can't save log_abs_Q..." << std::endl;
		return false;
	}
	return true;
}

static bool open_Qmn(real & Q, const real & c, const real & m, const real & n)
{
	if (!open_data(Q, generate_name(c, m, n, "Q")))
	{
		std::cout << "can't open Q..." << std::endl;
		return false;
	}
	return true;
}

static bool save_B2rmn(bool verbose, const real & c, const real & m, const real & n, real & n_B2r, const real & B2r_min)
{
	real lambda;
	real n_dr;
	std::vector<real> dr;
	real k1;
	real n_c2k;
	std::vector<real> c2k;
	real Q;
	std::vector<real> B2r;
	
	if (!open_lambdamn(lambda, c, m, n) ||
	    !open_drmn(n_dr, dr, c, m, n) ||
	    !open_kmn1(k1, c, m, n) ||
	    !open_c2kmn(n_c2k, c2k, c, m, n) ||
	    !open_Qmn(Q, c, m, n))
	{
		return false;
	}
	calculate_B2rmn(B2r, verbose, c, m, n, lambda, n_dr, dr, k1, n_c2k, c2k, Q, n_B2r, B2r_min);
	if (!save_data(generate_name(c, m, n, "B2r"), B2r))
	{
		std::cout << "can't save B2r..." << std::endl;
		return false;
	}
	if (!save_log_abs_data(generate_name(c, m, n, "log_abs_B2r"), B2r))
	{
		std::cout << "can't save log_abs_B2r..." << std::endl;
		return false;
	}
	return true;
}

static bool open_B2rmn(real & n_B2r, std::vector<real> & B2r, const real & c, const real & m, const real & n)
{
	if (!open_data(B2r, generate_name(c, m, n, "B2r")))
	{
		std::cout << "can't open B2r..." << std::endl;
		return false;
	}
	n_B2r = real((int)B2r.size());
	return true;
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
	real Q;
	real n_B2r;
	std::vector<real> B2r;
	real xi;
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
	real R2_31;
	real R2p_31;
	real R2_32;
	real R2p_32;
	real R2_log_abs_difference_1_2;
	real R2p_log_abs_difference_1_2;
	real R2_log_abs_difference_1_31;
	real R2p_log_abs_difference_1_31;
	real R2_log_abs_difference_1_32;
	real R2p_log_abs_difference_1_32;
	real R2_log_abs_difference_2_31;
	real R2p_log_abs_difference_2_31;
	real R2_log_abs_difference_2_32;
	real R2p_log_abs_difference_2_32;
	real R2_log_abs_difference_31_32;
	real R2p_log_abs_difference_31_32;
	real W;
	real log_W;
	real W_1_1_log_abs_error;
	real W_1_2_log_abs_error;
	real W_1_31_log_abs_error;
	real W_1_32_log_abs_error;
	real W_2_1_log_abs_error;
	real W_2_2_log_abs_error;
	real W_2_31_log_abs_error;
	real W_2_32_log_abs_error;
	
	if (!open_drmn(n_dr, dr, c, m, n) ||
	    !open_drmn_neg(n_dr_neg, dr_neg, c, m, n) ||
	    !open_Fmn(F, c, m, n) ||
	    !open_kmn1(k1, c, m, n) ||
	    !open_kmn2(k2, c, m, n) ||
	    !open_c2kmn(n_c2k, c2k, c, m, n) ||
	    !open_Qmn(Q, c, m, n) ||
	    !open_B2rmn(n_B2r, B2r, c, m, n))
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
			xi = i;
		}
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
		if (which.find("R2_31") != std::string::npos)
		{
			calculate_Rmn2_3(R2_31, R2p_31, verbose, c, m, n, Q, n_B2r, B2r, xi, R1_1, R1p_1);
		}
		else
		{
			R2_31 = real::NAN;
			R2p_31 = real::NAN;
		}
		if (which.find("R2_32") != std::string::npos)
		{
			calculate_Rmn2_3(R2_32, R2p_32, verbose, c, m, n, Q, n_B2r, B2r, xi, R1_2, R1p_2);
		}
		else
		{
			R2_32 = real::NAN;
			R2p_32 = real::NAN;
		}
		R2_log_abs_difference_1_2 = log(abs(R2_1 - R2_2));
		R2p_log_abs_difference_1_2 = log(abs(R2p_1 - R2p_2));
		R2_log_abs_difference_1_31 = log(abs(R2_1 - R2_31));
		R2p_log_abs_difference_1_31 = log(abs(R2p_1 - R2p_31));
		R2_log_abs_difference_1_32 = log(abs(R2_1 - R2_32));
		R2p_log_abs_difference_1_32 = log(abs(R2p_1 - R2p_32));
		R2_log_abs_difference_2_31 = log(abs(R2_2 - R2_31));
		R2p_log_abs_difference_2_31 = log(abs(R2p_2 - R2p_31));
		R2_log_abs_difference_2_32 = log(abs(R2_2 - R2_32));
		R2p_log_abs_difference_2_32 = log(abs(R2p_2 - R2p_32));
		R2_log_abs_difference_31_32 = log(abs(R2_31 - R2_32));
		R2p_log_abs_difference_31_32 = log(abs(R2p_31 - R2p_32));
		W = real::ONE / (c * (xi * xi + real::ONE));
		log_W = log(W);
		W_1_1_log_abs_error = log(abs(R1_1 * R2p_1 - R1p_1 * R2_1 - W));
		W_1_2_log_abs_error = log(abs(R1_1 * R2p_2 - R1p_1 * R2_2 - W));
		W_1_31_log_abs_error = log(abs(R1_1 * R2p_31 - R1p_1 * R2_31 - W));
		W_1_32_log_abs_error = log(abs(R1_1 * R2p_32 - R1p_1 * R2_32 - W));
		W_2_1_log_abs_error = log(abs(R1_2 * R2p_1 - R1p_2 * R2_1 - W));
		W_2_2_log_abs_error = log(abs(R1_2 * R2p_2 - R1p_2 * R2_2 - W));
		W_2_31_log_abs_error = log(abs(R1_2 * R2p_31 - R1p_2 * R2_31 - W));
		W_2_32_log_abs_error = log(abs(R1_2 * R2p_32 - R1p_2 * R2_32 - W));
		std::cout << i.get_string(p) << ","
		          << xi.get_string(p) << ","
		          << R1_1.get_string(p) << ","
		          << R1p_1.get_string(p) << ","
		          << R1_2.get_string(p) << ","
		          << R1p_2.get_string(p) << ","
		          << R1_log_abs_difference.get_string(p) << "," << R1p_log_abs_difference.get_string(p) << ","
		          << R2_1.get_string(p) << ","
		          << R2p_1.get_string(p) << ","
		          << R2_2.get_string(p) << ","
		          << R2p_2.get_string(p) << ","
		          << R2_31.get_string(p) << ","
		          << R2p_31.get_string(p) << ","
		          << R2_32.get_string(p) << ","
		          << R2p_32.get_string(p) << ","
		          << R2_log_abs_difference_1_2.get_string(p) << "," << R2p_log_abs_difference_1_2.get_string(p) << ","
		          << R2_log_abs_difference_1_31.get_string(p) << "," << R2p_log_abs_difference_1_31.get_string(p) << ","
		          << R2_log_abs_difference_1_32.get_string(p) << "," << R2p_log_abs_difference_1_32.get_string(p) << ","
		          << R2_log_abs_difference_2_31.get_string(p) << "," << R2p_log_abs_difference_2_31.get_string(p) << ","
		          << R2_log_abs_difference_2_32.get_string(p) << "," << R2p_log_abs_difference_2_32.get_string(p) << ","
		          << R2_log_abs_difference_31_32.get_string(p) << "," << R2p_log_abs_difference_31_32.get_string(p) << ","
		          << W.get_string(p) << ","
		          << log_W.get_string(p) << ","
		          << W_1_1_log_abs_error.get_string(p) << ","
		          << W_1_2_log_abs_error.get_string(p) << ","
		          << W_1_31_log_abs_error.get_string(p) << ","
		          << W_1_32_log_abs_error.get_string(p) << ","
		          << W_2_1_log_abs_error.get_string(p) << ","
		          << W_2_2_log_abs_error.get_string(p) << ","
		          << W_2_31_log_abs_error.get_string(p) << ","
		          << W_2_32_log_abs_error.get_string(p) << std::endl;
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
	real n_B2r;
	bool n_B2r_entered;
	real B2r_min;
	bool B2r_min_entered;
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
	n_B2r_entered = false;
	B2r_min_entered = false;
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
		else if (argument == "-n_B2r")
		{
			n_B2r = real(value);
			n_B2r_entered = true;
		}
		else if (argument == "-B2r_min")
		{
			B2r_min = real(value);
			B2r_min_entered = true;
		}
		else if (argument == "-a")
		{
			a = real(value);
			a_entered = true;
		}
		else if (argument == "-b")
		{
			b = real(value);
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
	else if (w == "Q")
	{
		save_Qmn(verbose, c, m, n);
	}
	else if (w == "B2r")
	{
		if (!n_B2r_entered || !B2r_min_entered)
		{
			std::cout << "no value of n_B2r and/or B2r_min was entered..." << std::endl;
			return 1;
		}
		save_B2rmn(verbose, c, m, n, n_B2r, B2r_min);
	}
	else if (w == "everything")
	{
		if (!n_dr_entered || !dr_min_entered || !n_dr_neg_entered || !dr_neg_min_entered || !n_c2k_entered || !c2k_min_entered || !n_B2r_entered || !B2r_min_entered)
		{
			std::cout << "no value of n_dr, dr_min, n_dr_neg, dr_neg_min, n_c2k, c2k_min, n_B2r, and/or B2r_min was entered..." << std::endl;
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
		save_Qmn(verbose, c, m, n);
		save_B2rmn(verbose, c, m, n, n_B2r, B2r_min);
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
