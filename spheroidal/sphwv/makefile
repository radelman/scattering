#
# Copyright (c) 2014, Ross Adelman, Nail A. Gumerov, and Ramani Duraiswami
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#

MPFR_DIR = /nfshomes/rna
CPP = g++
CPPFLAGS = -Wall -I$(MPFR_DIR)/include
OBJS = adder.o common_main.o common_spheroidal.o io.o obl_main.o obl_spheroidal.o pro_main.o pro_spheroidal.o real.o
PRO_OBJS = adder.o common_main.o common_spheroidal.o io.o pro_main.o pro_spheroidal.o real.o
OBL_OBJS = adder.o common_main.o common_spheroidal.o io.o obl_main.o obl_spheroidal.o real.o
LD = g++
LDFLAGS = -L$(MPFR_DIR)/lib -lmpfr -lgmpxx -lgmp -lm

%.o: %.cpp
	$(CPP) -o $@ -c $(CPPFLAGS) $^
sphwv: $(OBJS)
	$(LD) -o pro_sphwv $(LDFLAGS) $(PRO_OBJS)
	$(LD) -o obl_sphwv $(LDFLAGS) $(OBL_OBJS)
	rm $^
clean:
	rm pro_sphwv
	rm obl_sphwv
