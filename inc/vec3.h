#ifndef VEC3_H
#define VEC3_H

#include <math.h> /* for sqrtf */

typedef struct {
  float x, y, z;
} Vec3;

static inline Vec3 vec3_sub(const Vec3 *a, const Vec3 *b)
{
  Vec3 vec = {a->x - b->x, a->y - b->y, a->z - b->z};
  return vec;
}

static inline Vec3 vec3_add(const Vec3 *a, const Vec3 *b)
{
  Vec3 vec = {a->x + b->x, a->y + b->y, a->z + b->z};
  return vec;
}

static inline Vec3 vec3_cross(const Vec3 *a, const Vec3 *b)
{
  Vec3 vec = {
    a->y * b->z - a->z * b->y,
    a->z * b->x - a->x * b->z,
    a->x * b->y - a->y * b->x
  };
  return vec;
}

static inline float vec3_dot(const Vec3 *a, const Vec3 *b)
{
  return a->x * b->x + a->y * b->y + a->z * b->z;
}

static inline Vec3 vec3_scale(const Vec3 *a, float s)
{
  Vec3 vec = {a->x * s, a->y * s, a->z * s};
  return vec;
}

static inline Vec3 vec3_normalize(const Vec3 *a)
{
  float norm = sqrtf(a->x * a->x + a->y * a->y + a->z * a->z);
  Vec3 vec = {a->x / norm, a->y / norm, a->z / norm};
  return vec;
}

#endif /* VEC3_H */