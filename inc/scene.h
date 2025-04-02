#ifndef SCENE_H
#define SCENE_H

#include "common.h" /* for IN */
#include "vec3.h" /* for Vec3 */

#include <stdint.h> /* for uint32_t */
#include <stdlib.h> /* for free */

typedef struct {
  /* Number of vertices */
  uint32_t num_vertices;
  /* Vertices coordinates. Size [num_vertices] */
  Vec3 *vs;
  /* Number of triangles */
  uint32_t num_triangles;
  /* Triangles indices. Size [num_triangles * 3] */
  uint32_t *is;
  /* Index of the material in g_hrt_materials array defined in scene.h */
  uint32_t material_index;
  /* Global cartesian velocity vector */
  Vec3 velocity;
  /* Normals of the triangles. Size [num_triangles] */
  /* NOTE: This field is not stored in the HRT file.
  The loader does not load, allocate or calculate it. */
  Vec3 *ns;
} Mesh;

typedef struct {
  uint32_t num_meshes;
  Mesh *meshes;
} Scene;

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
} Material;

/** Free the memory allocated for a mesh fields.
 * 
 * \param mesh pointer to the mesh to free
 */
static inline void free_mesh(Mesh* mesh) {
  free(mesh->vs);
  free(mesh->is);
  free(mesh->ns);
}
/** Deep free the scene. Frees all the meshes.
 * 
 * \param scene pointer to the scene to free
 */
static inline void free_scene(Scene* scene) {
  for (uint32_t i = 0; i < scene->num_meshes; i++) {
    free_mesh(&scene->meshes[i]);
  }
  free(scene->meshes);
}

/** Save a scene to a HRT file. (See README for details)
 *
 * NOTE: This function does not save the normals of the triangles (Mesh.ns field).
 * 
 * \param scene pointer to the scene to save
 * \param filepath path to the HRT file. Will be overwritten if exists.
 */
void scene_save(IN Scene* scene, IN const char* filepath);

/** Load a mesh from a HRT file.
 * 
 * NOTE: This function does not load, allocate or calculate
 * the normals of the triangles (Mesh.ns field).
 * 
 * \param filepath path to the HRT file
 * \return the loaded scene
 */
Scene scene_load(IN const char *filepath);

#endif /* SCENE_H */
