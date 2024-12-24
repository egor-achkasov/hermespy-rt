all:
	gcc -g -O0 -fno-builtin compute_paths.c test.c -lm -o test
