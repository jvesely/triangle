BINARY=triangle
OBJS=main.o

MAKEDEPEND=g++
CXX = g++
CPPFLAGS = -I. -std=c++11 -D __STDC_FORMAT_MACROS
CXXFLAGS = -O3 -Wall -Wextra -fmax-errors=3
LDFLAGS =

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
