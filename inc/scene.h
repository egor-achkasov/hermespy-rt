#ifndef SCENE_H
#define SCENE_H

#include <stdint.h> /* for uint32_t */

typedef struct {
  /* Number of vertices */
  uint32_t num_vertices;
  /* Vertices coordinates. Size [num_vertices * 3] */
  float *vs;
  /* Number of triangles */
  uint32_t num_triangles;
  /* Triangles indices. Size [num_triangles * 3] */
  uint32_t *is;
  /* Index of the material in g_hrt_materials array defined in scene.h */
  uint32_t material_index;
  /* Global cartesian velocity vector. Size [3] */
  float velocity[3];
} HRT_Mesh;

typedef struct {
  uint32_t num_meshes;
  HRT_Mesh *meshes;
} HRT_Scene;

typedef struct {
  /* Number of characters in the name */
  uint32_t name_sz;
  /* Name of the material. Not null-terminated. */
  const char *name;
  /* Relative permitivity and conductivity properties.
   * Refer to ITU-R P.2040-3 Table 3.
   */
  float a, b, c, d;
  /* Scattering coefficient.
   * Ratio between specular and diffuse reflection power.
   * Must be in [0, 1].
   */
  float s;
  /* Scattering pattern distribution ratios.
   * Each must be in [0, 1].
   * The sum of all ratios must be 1.
   * s1: around specular (Directive Model)
   * s2: around surface normal (Lambertian Model)
   * s3: around incident direction (Backward Lobe Model)
   */
  float s1, s2, s3;
  /* Directive model parameters.
   * s1_alpha: lobe width integer parameter (TODO ref eq).
   * Must be > 0.
   */
  uint8_t s1_alpha;
  /* Backward lobe model parameters.
   * s3_alpha: lobe width integer parameter (TODO ref eq).
   * Must be > 0.
   */
  uint8_t s3_alpha;
} HRT_Material;

#endif /* SCENE_H */
