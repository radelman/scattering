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

#include <fstream>
#include "io.hpp"
#include "real.hpp"
#include <string>
#include <vector>

bool save_data(const std::string & name, const std::vector<real> & data)
{
	std::ofstream out;
	
	out.open(name.c_str());
	if (out.fail())
	{
		return false;
	}
	for (int i = 0; i < (int)data.size(); ++i)
	{
		out << data[i].get_string(0) << std::endl;
	}
	out.close();
	return true;
}

bool save_data(const std::string & name, const real & data)
{
	std::vector<real> container;
	
	container.clear();
	container.push_back(data);
	return save_data(name, container);
}

bool save_log_abs_data(const std::string & name, const std::vector<real> & data)
{
	std::vector<real> log_abs_data;
	
	for (int i = 0; i < (int)data.size(); ++i)
	{
		log_abs_data.push_back(real::ZERO);
	}
	for (int i = 0; i < (int)data.size(); ++i)
	{
		log_abs_data[i] = log(abs(data[i]));
	}
	return save_data(name, log_abs_data);
}

bool save_log_abs_data(const std::string & name, const real & data)
{
	return save_data(name, log(abs(data)));
}

bool open_data(std::vector<real> & data, const std::string & name)
{
	std::ifstream in;
	std::string string;
	
	data.clear();
	in.open(name.c_str());
	if (in.fail())
	{
		return false;
	}
	while (true)
	{
		// Read the next line into string.  This used to read the next line
		// into a character array using in.getline, which could fail if there
		// were more characters than the length of the array.
		getline(in, string);
		if (in.eof())
		{
			break;
		}
		data.push_back(real(string));
	}
	in.close();
	return true;
}

bool open_data(real & data, const std::string & name)
{
	std::vector<real> container;
	
	if (!open_data(container, name))
	{
		return false;
	}
	data = container[0];
	return true;
}
