all:
	gcc -g -O0 compute_paths.c test.c -lxml2 -lm -o test
