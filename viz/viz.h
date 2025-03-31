/* Visualization utils.
The implementation files are located in viz/
All of them require OpenGL and GLUT
add -lGL -lGLU -lglut to compile
*/

#ifndef VIZ_H
#define VIZ_H

#include "../inc/scene.h" /* for Scene */
#include "../inc/compute_paths.h" /* for RaysInfo */
#include "../inc/common.h" /* for IN */

#include <stdint.h> /* for uint32_t */

/** Visualize rays in a 3D scene.
 * 
 * \param raysInfo pointer to a RaysInfo structure
 * \param scene pointer to a loaded scene in HRT format (see scene.h). Can be NULL
 * \param numTx number of transmitters
 */
void vizrays(
  IN RaysInfo *raysInfo,
  IN Scene *scene,
  IN uint32_t numTx
);

#endif
