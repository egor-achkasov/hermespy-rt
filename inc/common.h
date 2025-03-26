#ifndef COMMON_H
#define COMMON_H

#define IN
#define OUT

#include <stdlib.h> /* for exit */

#define FREE_POINTERS(...) \
  do {                            \
    void *ptrs[] = {__VA_ARGS__};\
    size_t cnt = sizeof(ptrs) / sizeof(ptrs[0]); \
    for (size_t i = 0; i < cnt; i++) \
      free(ptrs[i]);          \
  } while (0)

/* Variadic arguments after rc are the pointers to free */
#define PERROR_CLEANUP_EXIT(msg, rc, ...) \
  do {                            \
    perror(msg);                \
    FREE_POINTERS(__VA_ARGS__); \
    exit(rc);                    \
  } while (0)

#endif
