#include <stdint.h>

typedef struct {
  int name_size;
  char *name;
  float a, b, c, d;
  float s;
  uint8_t alpha;
} Material;
