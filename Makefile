all:
	gcc -g -O0 -Wno-builtin-declaration-mismatch compute_paths.c test.c -lm -o test
