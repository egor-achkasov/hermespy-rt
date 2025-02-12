#include <stdint.h> /* for uint32_t */

/*******************************************************************/
/*                              Scene                              */
/*******************************************************************/

typedef struct {
  /* Number of vertices */
  uint32_t num_vertices;
  /* Vertices coordinates. Size [num_vertices * 3] */
  float *vs;
  /* Number of triangles */
  uint32_t num_triangles;
  /* Triangles indices. Size [num_triangles * 3] */
  uint32_t *is;
  uint32_t material_index;
} HRT_Mesh;

typedef struct {
  uint32_t num_meshes;
  HRT_Mesh *meshes;
} HRT_Scene;

/*******************************************************************/
/*                            Materials                            */
/*******************************************************************/

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

/* number of defined materials */
#define NUM_G_MATERIALS 17
HRT_Material g_hrt_materials[NUM_G_MATERIALS] = {
  /* 0 */
  {3, "air",
    1.f, 0.f, 0.f, 0.001f,
    0.1f, 0.5f, 0.3f, 0.2f,
    2, 2},
  /* 1 */
  {8, "concrete",
    5.24f, 0.f, 0.0462f, 0.7822f,
    0.5f, 0.33f, 0.34f, 0.33f,
    4, 4},
  /* 2 */
  {5, "brick",
    3.91f, 0.f, 0.0238f, 0.16f,
    0.4f, 0.4f, 0.3f, 0.3f,
    3, 3},
  /* 3 */
  {12, "plasterboard",
    2.73f, 0.f, 0.0085f, 0.9395f,
    0.3f, 0.4f, 0.4f, 0.2f,
    3, 3},
  /* 4 */
  {4, "wood",
    1.99f, 0.f, 0.0047f, 1.0718f,
    0.2f, 0.5f, 0.3f, 0.2f,
    2, 2},
  /* 5 */
  {5, "glass",
    6.31f, 0.f, 0.0036f, 1.3394f,
    0.3f, 0.4f, 0.4f, 0.2f,
    3, 3},
  /* 6 */
  {5, "glass",
    5.79f, 0.f, 0.0004f, 1.658f,
    0.3f, 0.4f, 0.4f, 0.2f,
    3, 3},
  /* 7 */
  {13, "ceiling board",
    1.48f, 0.f, 0.0011f, 1.0750f,
    0.2f, 0.5f, 0.3f, 0.2f,
    2, 2},
  /* 8 */
  {13, "ceiling board",
    1.52f, 0.f, 0.0029f, 1.029f,
    0.2f, 0.5f, 0.3f, 0.2f,
    2, 2},
  /* 9 */
  {9, "chipboard",
    2.58f, 0.f, 0.0217f, 0.7800f,
    0.4f, 0.4f, 0.3f, 0.3f,
    3, 3},
  /* 10 */
  {7, "plywood",
    2.71f, 0.f, 0.33f, 0.f,
    0.3f, 0.5f, 0.3f, 0.2f,
    3, 3},
  /* 11 */
  {6, "marble",
    7.074f, 0.f, 0.0055f,
    0.9262f, 0.3f, 0.4f, 0.4f, 0.2f,
    3, 3},
  /* 12 */
  {10, "floorboard",
    3.66f, 0.f, 0.0044f, 1.3515f,
    0.3f, 0.4f, 0.4f, 0.2f,
    3, 3},
  /* 13 */
  {5, "metal",
    1.f, 0.f, 10000000.f, 0.f,
    0.f, 0.f, 1.f, 0.f,
    1, 1},
  /* 14 */
  {15, "very dry ground",
    3.f, 0.f, 0.00015f, 2.52f,
    0.4f, 0.3f, 0.4f, 0.3f,
    4, 4},
  /* 15 */
  {17, "medium dry ground",
    15.f, -0.1f, 0.035f, 1.63f,
    0.5f, 0.33f, 0.34f, 0.33f,
    4, 4},
  /* 16 */
  {10, "wet ground",
    30.f, -0.4f, 0.15f, 1.30f,
    0.5f, 0.33f, 0.34f, 0.33f,
    4, 4}
};
typedef enum {
  MATERIAL_AIR = 0,
  MATERIAL_CONCRETE = 1,
  MATERIAL_BRICK = 2,
  MATERIAL_PLASTERBOARD = 3,
  MATERIAL_WOOD = 4,
  MATERIAL_GLASS1 = 5,
  MATERIAL_GLASS2 = 6,
  MATERIAL_CEILING_BOARD1 = 7,
  MATERIAL_CEILING_BOARD2 = 8,
  MATERIAL_CHIPBOARD = 9,
  MATERIAL_PLYWOOD = 10,
  MATERIAL_MARBLE = 11,
  MATERIAL_FLOORBOARD = 12,
  MATERIAL_METAL = 13,
  MATERIAL_VERY_DRY_GROUND = 14,
  MATERIAL_MEDIUM_DRY_GROUND = 15,
  MATERIAL_WET_GROUND = 16
} MaterialIndex;

