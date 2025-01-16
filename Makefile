all:
	gcc -g -O0 -Wno-builtin-declaration-mismatch compute_paths.c test.c -lm -o test

py:
	g++ -g -I/usr/include/python3.11 -I/home/hi/proj/hermespy/env/lib/python3.11/site-packages/pybind11/include -O3 -Wall -Wno-uninitialized -shared -std=c++17 -fPIC compute_paths.c compute_paths_pybind11.cpp -o rt$(python3.11-config --extension-suffix)
