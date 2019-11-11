# Build the samtools_strnum_cmp utility

# The compiler:  gcc for C, g++ for C++
CC = gcc

# compiler flags:
#  -g     Adds debugging information to the executable file.
#  -Wall  Turns on most, but not all, compiler warnings.
#  -Wno-unused-function  Since we are making a library, all the functions are unused. So turn this off.
#  -fPIC -shared         Needed so this will work as a library for python.
#  -O3    A high optimization level.

CFLAGS = -Wall -Wno-unused-function -fPIC -shared -O3

samtools_strnum_cmp.so: samtools_strnum_cmp.c
		$(CC) $(CFLAGS) -o samtools_strnum_cmp.so samtools_strnum_cmp.c

clean:
		rm -f *\.so

