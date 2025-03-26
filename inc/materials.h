#ifndef MATERIALS_H
#define MATERIALS_H

#include "scene.h" /* for Material */

#include <stdint.h> /* for uint32_t */
#include <string.h> /* for strncmp */

/* number of defined materials */
#define NUM_G_MATERIALS 17
extern Material g_materials[NUM_G_MATERIALS];
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

/** Get the material index from the material name */
extern MaterialIndex get_material_index(const char *name);

#endif /* MATERIALS_H */
