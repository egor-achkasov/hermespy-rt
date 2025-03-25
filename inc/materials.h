#ifndef MATERIALS_H
#define MATERIALS_H

#include "scene.h" /* for HRT_Material */

#include <stdint.h> /* for uint32_t */
#include <string.h> /* for strncmp */

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

/* Map material name strings to material indices */
typedef struct {
  /* Material name. Null-terminated string. Lower case. */
  const char *name;
  MaterialIndex index;
} MaterialMap;

MaterialMap material_map[] = {
  {"air", MATERIAL_AIR},
  {"concrete", MATERIAL_CONCRETE},
  {"brick", MATERIAL_BRICK},
  {"plasterboard", MATERIAL_PLASTERBOARD},
  {"wood", MATERIAL_WOOD},
  {"glass1", MATERIAL_GLASS1},
  {"glass2", MATERIAL_GLASS2},
  {"ceiling_board1", MATERIAL_CEILING_BOARD1},
  {"ceiling_board2", MATERIAL_CEILING_BOARD2},
  {"chipboard", MATERIAL_CHIPBOARD},
  {"plywood", MATERIAL_PLYWOOD},
  {"marble", MATERIAL_MARBLE},
  {"floorboard", MATERIAL_FLOORBOARD},
  {"metal", MATERIAL_METAL},
  {"very_dry_ground", MATERIAL_VERY_DRY_GROUND},
  {"medium_dry_ground", MATERIAL_MEDIUM_DRY_GROUND},
  {"wet_ground", MATERIAL_WET_GROUND}
};

MaterialIndex get_material_index(const char *name) {
  for (uint32_t i = 0; i < NUM_G_MATERIALS; ++i)
    if (strcmp(material_map[i].name, name) == 0)
      return material_map[i].index;
  return MATERIAL_AIR; /* Invalid material. Default to air */
}

#endif /* MATERIALS_H */
