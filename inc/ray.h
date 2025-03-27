#ifndef RAY_H
#define RAY_H

#include "vec3.h" /* for Vec3 */

typedef struct {
    Vec3 o; /* origin */
    Vec3 d; /* direction */
} Ray;

#endif
