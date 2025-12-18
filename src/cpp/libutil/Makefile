# libutil Makefile

# Compiler and flags
CXX ?= g++
CFLAGS ?= -O2 -Wall -fPIC
HEALPIX_INCLUDE_CPP = -I/home/dwatts/Downloads/Healpix_3.83/include


# Targets
TARGET = libutil.a
SRCS   = sslice.cpp stoken.cpp spng.cpp scolorizer.cpp stimer.cpp serror.cpp seph.cpp hpixdata.cpp
OBJS   = $(SRCS:.cpp=.o)

# Default rule
all: $(TARGET)

# Build the static library
$(TARGET): $(OBJS)
	ar crs $@ $^

# Compile source files to object files
%.o: %.cpp
	$(CXX) $(CFLAGS) -I. $(HEALPIX_INCLUDE_CPP) -I../ephcom -c $< -o $@

# Clean build products
clean:
	rm -f $(TARGET) $(OBJS)

.PHONY: all clean

