# Copyright (c) 2014 Jan Vesely
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# - Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# - Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# - The name of the author may not be used to endorse or promote products
#   derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
# NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
BINARY=triangle
OBJS=main.o

MAKEDEPEND=g++
CXX = g++
CPPFLAGS = -I. -std=c++11 -D __STDC_FORMAT_MACROS
CXXFLAGS = -O3 -Wall -Wextra -fmax-errors=3
LDFLAGS  = -lboost_program_options

ifeq ($(DEBUG), TRUE)
	CXXFLAGS += -g
endif

ifeq ($(STATIC), TRUE)
	CXXFLAGS += -static
endif

FLAGS= $(CPPFLAGS) $(CXXFLAGS)

all: $(BINARY)

debug:
	make DEBUG=TRUE

$(BINARY): $(OBJS)
	$(CXX) $(FLAGS) $^ $(LDFLAGS) -o $@

%.o: %.cpp Makefile
	$(CXX) $(FLAGS) -c $< -o $@

%.d: %.cpp Makefile
	$(MAKEDEPEND) $(CPPFLAGS) -MM -MMD $<

clean:
	rm -vf *.o *.d $(BINARY)

-include $(OBJS:.o=.d)
