CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -g

SRCS = my_exception.cc cholesky.cc ll_cholesky.cc ll_blc_cho.cc \
block_cho.cc
OBJS = $(SRCS:.cc=.o)
DEPS = my_exception.h cholesky.h ll_cholesky.h ll_blc_cho.h block_cho.h

TARGET = ll_blc_cho

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

%.o: %.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
