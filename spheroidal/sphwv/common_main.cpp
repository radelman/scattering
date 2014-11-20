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
#include "io.hpp"
#include <iostream>
#include "real.hpp"
#include <string>

bool save_lambdamn(bool verbose, const real & c, const real & m, const real & n)
{
	real lambda_approx;
	real lambda;
	
	if (!open_data(lambda_approx, generate_name(c, m, n, "lambda_approx")))
	{
		std::cout << "can't open lambda_approx..." << std::endl;
		return false;
	}
	calculate_lambdamn(lambda, verbose, c, m, n, lambda_approx);
	if (!save_data(generate_name(c, m, n, "lambda"), lambda))
	{
		std::cout << "can't save lambda..." << std::endl;
		return false;
	}
	if (!save_log_abs_data(generate_name(c, m, n, "log_abs_lambda"), lambda))
	{
		std::cout << "can't save log_abs_lambda..." << std::endl;
		return false;
	}
	return true;
}

bool open_lambdamn(real & lambda, const real & c, const real & m, const real & n)
{
	if (!open_data(lambda, generate_name(c, m, n, "lambda")))
	{
		std::cout << "can't open lambda..." << std::endl;
		return false;
	}
	return true;
}

bool save_drmn(bool verbose, const real & c, const real & m, const real & n, real & n_dr, const real & dr_min)
{
	real lambda;
	std::vector<real> dr;
	
	if (!open_lambdamn(lambda, c, m, n))
	{
		return false;
	}
	calculate_drmn(dr, verbose, c, m, n, lambda, n_dr, dr_min);
	if (!save_data(generate_name(c, m, n, "dr"), dr))
	{
		std::cout << "can't save dr..." << std::endl;
		return false;
	}
	if (!save_log_abs_data(generate_name(c, m, n, "log_abs_dr"), dr))
	{
		std::cout << "can't save log_abs_dr..." << std::endl;
		return false;
	}
	return true;
}

bool open_drmn(real & n_dr, std::vector<real> & dr, const real & c, const real & m, const real & n)
{
	if (!open_data(dr, generate_name(c, m, n, "dr")))
	{
		std::cout << "can't open dr..." << std::endl;
		return false;
	}
	n_dr = real((int)dr.size());
	return true;
}

bool save_drmn_neg(bool verbose, const real & c, const real & m, const real & n, real & n_dr_neg, const real & dr_neg_min)
{
	real lambda;
	real n_dr;
	std::vector<real> dr;
	std::vector<real> dr_neg;
	
	if (!open_lambdamn(lambda, c, m, n) ||
	    !open_drmn(n_dr, dr, c, m, n))
	{
		return false;
	}
	calculate_drmn_neg(dr_neg, verbose, c, m, n, lambda, n_dr, dr, n_dr_neg, dr_neg_min);
	if (!save_data(generate_name(c, m, n, "dr_neg"), dr_neg))
	{
		std::cout << "can't save dr_neg..." << std::endl;
		return false;
	}
	if (!save_log_abs_data(generate_name(c, m, n, "log_abs_dr_neg"), dr_neg))
	{
		std::cout << "can't save log_abs_dr_neg..." << std::endl;
		return false;
	}
	return true;
}

bool open_drmn_neg(real & n_dr_neg, std::vector<real> & dr_neg, const real & c, const real & m, const real & n)
{
	if (!open_data(dr_neg, generate_name(c, m, n, "dr_neg")))
	{
		std::cout << "can't open dr_neg..." << std::endl;
		return false;
	}
	n_dr_neg = real((int)dr_neg.size());
	return true;
}

bool save_Nmn(bool verbose, const real & c, const real & m, const real & n)
{
	real lambda;
	real n_dr;
	std::vector<real> dr;
	real N;
	
	if (!open_lambdamn(lambda, c, m, n) ||
	    !open_drmn(n_dr, dr, c, m, n))
	{
		return false;
	}
	N = calculate_Nmn(verbose, c, m, n, lambda, n_dr, dr);
	if (!save_data(generate_name(c, m, n, "N"), N))
	{
		std::cout << "can't save N..." << std::endl;
		return false;
	}
	if (!save_log_abs_data(generate_name(c, m, n, "log_abs_N"), N))
	{
		std::cout << "can't save log_abs_N..." << std::endl;
		return false;
	}
	return true;
}

bool open_Nmn(real & N, const real & c, const real & m, const real & n)
{
	if (!open_data(N, generate_name(c, m, n, "N")))
	{
		std::cout << "can't open N..." << std::endl;
		return false;
	}
	return true;
}

bool save_Fmn(bool verbose, const real & c, const real & m, const real & n)
{
	real lambda;
	real n_dr;
	std::vector<real> dr;
	real F;
	
	if (!open_lambdamn(lambda, c, m, n) ||
	    !open_drmn(n_dr, dr, c, m, n))
	{
		return false;
	}
	F = calculate_Fmn(verbose, c, m, n, lambda, n_dr, dr);
	if (!save_data(generate_name(c, m, n, "F"), F))
	{
		std::cout << "can't save F..." << std::endl;
		return false;
	}
	if (!save_log_abs_data(generate_name(c, m, n, "log_abs_F"), F))
	{
		std::cout << "can't save log_abs_F..." << std::endl;
		return false;
	}
	return true;
}

bool open_Fmn(real & F, const real & c, const real & m, const real & n)
{
	if (!open_data(F, generate_name(c, m, n, "F")))
	{
		std::cout << "can't open F..." << std::endl;
		return false;
	}
	return true;
}

bool save_kmn1(bool verbose, const real & c, const real & m, const real & n)
{
	real n_dr;
	std::vector<real> dr;
	real F;
	real k1;
	
	if (!open_drmn(n_dr, dr, c, m, n) ||
	    !open_Fmn(F, c, m, n))
	{
		return false;
	}
	k1 = calculate_kmn1(verbose, c, m, n, n_dr, dr, F);
	if (!save_data(generate_name(c, m, n, "k1"), k1))
	{
		std::cout << "can't save k1..." << std::endl;
		return false;
	}
	if (!save_log_abs_data(generate_name(c, m, n, "log_abs_k1"), k1))
	{
		std::cout << "can't save log_abs_k1.." << std::endl;
		return false;
	}
	return true;
}

bool open_kmn1(real & k1, const real & c, const real & m, const real & n)
{
	if (!open_data(k1, generate_name(c, m, n, "k1")))
	{
		std::cout << "can't open k1..." << std::endl;
		return false;
	}
	return true;
}

bool save_kmn2(bool verbose, const real & c, const real & m, const real & n)
{
	real lambda;
	real n_dr;
	std::vector<real> dr;
	real n_dr_neg;
	std::vector<real> dr_neg;
	real F;
	real k2;
	
	if (!open_lambdamn(lambda, c, m, n) ||
	    !open_drmn(n_dr, dr, c, m, n) ||
	    !open_drmn_neg(n_dr_neg, dr_neg, c, m, n) ||
	    !open_Fmn(F, c, m, n))
	{
		return false;
	}
	k2 = calculate_kmn2(verbose, c, m, n, lambda, n_dr, dr, n_dr_neg, dr_neg, F);
	if (!save_data(generate_name(c, m, n, "k2"), k2))
	{
		std::cout << "can't save k2..." << std::endl;
		return false;
	}
	if (!save_log_abs_data(generate_name(c, m, n, "log_abs_k2"), k2))
	{
		std::cout << "can't save log_abs_k2..." << std::endl;
		return false;
	}
	return true;
}

bool open_kmn2(real & k2, const real & c, const real & m, const real & n)
{
	if (!open_data(k2, generate_name(c, m, n, "k2")))
	{
		std::cout << "can't open k2..." << std::endl;
		return false;
	}
	return true;
}

bool save_c2kmn(bool verbose, const real & c, const real & m, const real & n, real & n_c2k, const real & c2k_min)
{
	real lambda;
	real n_dr;
	std::vector<real> dr;
	std::vector<real> c2k;
	
	if (!open_lambdamn(lambda, c, m, n) ||
	    !open_drmn(n_dr, dr, c, m, n))
	{
		return false;
	}
	c2k.clear();
	calculate_c2kmn(c2k, verbose, c, m, n, lambda, n_dr, dr, n_c2k, c2k_min);
	if (!save_data(generate_name(c, m, n, "c2k"), c2k))
	{
		std::cout << "can't save c2k..." << std::endl;
		return false;
	}
	if (!save_log_abs_data(generate_name(c, m, n, "log_abs_c2k"), c2k))
	{
		std::cout << "can't save log_abs_c2k..." << std::endl;
		return false;
	}
	return true;
}

bool open_c2kmn(real & n_c2k, std::vector<real> & c2k, const real & c, const real & m, const real & n)
{
	if (!open_data(c2k, generate_name(c, m, n, "c2k")))
	{
		std::cout << "can't open c2k..." << std::endl;
		return false;
	}
	n_c2k = real((int)c2k.size());
	return true;
}

bool save_Smn1(bool verbose, const real & c, const real & m, const real & n, const real & a, const real & b, const real & d, const std::string & arg_type, int p)
{
	real n_dr;
	std::vector<real> dr;
	real N;
	real n_c2k;
	std::vector<real> c2k;
	real eta;
	real S1_1;
	real S1p_1;
	real S1_2;
	real S1p_2;
	real S1_log_abs_difference;
	real S1p_log_abs_difference;

	if (!open_drmn(n_dr, dr, c, m, n) ||
	    !open_Nmn(N, c, m, n) ||
	    !open_c2kmn(n_c2k, c2k, c, m, n))
	{
		return false;
	}
	for (real i = a; i <= b; i = i + d)
	{
		if (arg_type == "eta")
		{
			eta = i;
		}
		else
		{
			eta = cos(i * real::PI);
		}
		calculate_Smn1_1(S1_1, S1p_1, verbose, c, m, n, n_dr, dr, eta);
		S1_1 = S1_1 / pow(N, real::ONE / real::TWO);
		S1p_1 = S1p_1 / pow(N, real::ONE / real::TWO);
		calculate_Smn1_2(S1_2, S1p_2, verbose, c, m, n, n_c2k, c2k, eta);
		S1_2 = S1_2 / pow(N, real::ONE / real::TWO);
		S1p_2 = S1p_2 / pow(N, real::ONE / real::TWO);
		S1_log_abs_difference = log(abs(S1_1 - S1_2));
		S1p_log_abs_difference = log(abs(S1p_1 - S1p_2));
		std::cout << i.get_string(p) << ","
		          << eta.get_string(p) << ","
		          << S1_1.get_string(p) << "," << S1p_1.get_string(p) << ","
		          << S1_2.get_string(p) << "," << S1p_2.get_string(p) << ","
		          << S1_log_abs_difference.get_string(p) << "," << S1p_log_abs_difference.get_string(p) << std::endl;
	}
	return true;
}

//
// This is the main entrypoint of the program, both for pro_sphwv and
// obl_sphwv.  It sets the default precision in MPFR and the maximum number of
// reals that can be used at any one time based on the -max_memory argument.
// At the end, it calls parse_args, which is defined differently for pro_sphwv
// and obl_sphwv.
//
int main(int argc, char **argv)
{
	std::string argument;
	std::string value;
	int max_memory;
	bool max_memory_entered;
	int precision;
	bool precision_entered;
	
	max_memory_entered = false;
	precision_entered = false;
	for (int i = 1; i < argc; i = i + 2)
	{
		argument = std::string(argv[i]);
		value = std::string(argv[i + 1]);
		if (argument == "-max_memory")
		{
			max_memory = std::atoi(value.c_str());
			max_memory_entered = true;
		}
		else if (argument == "-precision")
		{
			precision = std::atoi(value.c_str());
			precision_entered = true;
		}
	}
	if (!max_memory_entered || !precision_entered)
	{
		std::cout << "no value of max_memory and/or precision was entered..." << std::endl;
		return 1;
	}
	real::begin(precision, (int)((double)max_memory * (8000000.0 / (double)precision)));
	complex::begin();
	return parse_args(argc, argv);
}
