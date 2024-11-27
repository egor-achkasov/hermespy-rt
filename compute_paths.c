#include "compute_paths.h"

#include <stddef.h> /* for size_t */
#include <stdlib.h> /* for exit, malloc, free */
#include <libgen.h> /* for dirname */
#include <string.h> /* for strlen, sprintf */
#include <stdio.h>  /* for fopen, FILE, fclose */
#include <math.h>   /* for sin, cos, sqrt */

#define PI 3.14159265358979323846
#define EPS 1e-6

typedef struct {
    float x, y, z;
} Vec3;

typedef struct {
    Vec3 o, d;
} Ray;

typedef struct {
    float a, b, c, d;
} RadioMaterial;

typedef struct {
    Vec3 *vertices;
    size_t num_vertices;
    RadioMaterial *rms;
    size_t num_rms;
    int32_t *indices;
    int32_t *rm_indices;
    size_t num_indices;
    Vec3 *normals;
} Mesh;

/* ==== VECTOR OPERATIONS ==== */

Vec3 vec3_sub(const Vec3 *a, const Vec3 *b)
{
    return (Vec3){a->x - b->x, a->y - b->y, a->z - b->z};
}
Vec3 vec3_add(const Vec3 *a, const Vec3 *b)
{
    return (Vec3){a->x + b->x, a->y + b->y, a->z + b->z};
}
Vec3 vec3_cross(const Vec3 *a, const Vec3 *b)
{
    return (Vec3){
        a->y * b->z - a->z * b->y,
        a->z * b->x - a->x * b->z,
        a->x * b->y - a->y * b->x
    };
}
float vec3_dot(const Vec3 *a, const Vec3 *b)
{
    return a->x * b->x + a->y * b->y + a->z * b->z;
}
Vec3 vec3_scale(const Vec3 *a, float s)
{
    return (Vec3){a->x * s, a->y * s, a->z * s};
}

/* ==== MESH LOADING ==== */

/** Load a mesh from a PLY file.
 * 
 * ply
 * format binary_little_endian 1.0
 * element vertex %d\n", num_vertices);
 * property float x
 * property float y
 * property float z
 * element radio_material %d\n", num_materials);
 * property list uint char name
 * property float a
 * property float b
 * property float c
 * property float d
 * element face %d\n", num_faces);
 * property list uchar uint vertex_index
 * property uint material_index
 * end_header
 * 
 * All the faces must be triangles.
 * 
 * \param mesh_filepath path to the PLY file
 * \return the loaded mesh
 */
Mesh load_mesh_ply(const char *mesh_filepath)
{
    FILE *f = fopen(mesh_filepath, "rb");
    if (!f) {
        fprintf(stderr, "Could not open file %s\n", mesh_filepath);
        exit(8);
    }

    char buff[128];
    size_t num_vertices = 0;
    size_t num_faces = 0;
    size_t num_materials = 0;

    /* HEADER */
    /* ply */
    fread(buff, 1, 4, f);
    if (strncmp(buff, "ply\n", 4)) exit(8);
    /* format binary_little_endian 1.0 */
    fread(buff, 1, 32, f);
    if (strncmp(buff, "format binary_little_endian 1.0\n", 32)) exit(8);

    /* element vertex <num_vertices> */
    fgets(buff, 128, f);
    if (strncmp(buff, "element vertex ", 15)) exit(8);
    sscanf(buff, "element vertex %zu", &num_vertices);
    /* property float x */
    /* property float y */
    /* property float z */
    fread(buff, 1, 17*3, f);
    if (strncmp(buff, "property float x\nproperty float y\nproperty float z\n", 17*3))
        exit(8);

    /* element radio_material */
    fgets(buff, 128, f);
    if (strncmp(buff, "element radio_material ", 23)) exit(8);
    sscanf(buff, "element radio_material %zu", &num_materials);
    /* property list uint char name */
    fread(buff, 1, 29, f);
    if (strncmp(buff, "property list uint char name\n", 29)) exit(8);
    /* property float a */
    /* property float b */
    /* property float c */
    /* property float d */
    fread(buff, 1, 17*4, f);
    if (strncmp(buff, "property float a\nproperty float b\nproperty float c\nproperty float d\n", 17*4))
        exit(8);

    /* element face <num_faces> */
    fgets(buff, 128, f);
    if (strncmp(buff, "element face ", 13)) exit(8);
    sscanf(buff, "element face %zu", &num_faces);
    /* property list uchar int vertex_index */
    fread(buff, 1, 38, f);
    if (strncmp(buff, "property list uchar uint vertex_index\n", 38))
        exit(8);
    /* property uint material_index */
    fread(buff, 1, 29, f);
    if (strncmp(buff, "property uint material_index\n", 29)) exit(8);

    /* end_header */
    fread(buff, 1, 11, f);
    if (strncmp(buff, "end_header\n", 11)) exit(8);

    /* Init mesh */
    Mesh mesh;
    mesh.num_vertices = num_vertices;
    mesh.num_rms = num_materials;
    mesh.num_indices = num_faces * 3;
    mesh.vertices = (Vec3*)malloc(num_vertices * sizeof(Vec3));
    mesh.rms = (RadioMaterial*)malloc(num_materials * sizeof(RadioMaterial));
    mesh.indices = (int32_t*)malloc(mesh.num_indices * sizeof(int32_t));
    mesh.rm_indices = (int32_t*)malloc(num_faces * sizeof(int32_t));
    mesh.normals = (Vec3*)malloc(num_faces * sizeof(Vec3));

    /* VERTICES */
    for (size_t i = 0; i < mesh.num_vertices; ++i)
        fread(&mesh.vertices[i], sizeof(Vec3), 1, f);

    /* RADIO_MATERIALS */
    unsigned int name_size;
    for (size_t i = 0; i < num_materials; ++i) {
        /* skip name */
        fread(&name_size, sizeof(unsigned int), 1, f);
        fseek(f, name_size, SEEK_CUR);
        fread(&mesh.rms[i], sizeof(RadioMaterial), 1, f);
    }

    /* INDICES */
    /* also calculate normals */
    uint8_t n;
    Vec3 *v1, *v2, *v3, u, v;
    for (size_t i = 0; i < mesh.num_indices; i += 3) {
        fread(&n, 1, 1, f);
        if (n != 3) {
            fprintf(stderr, "Only trianglar faces are supported\n");
            exit(8);
        }
        fread(&mesh.indices[i], sizeof(int32_t), 3, f);
        fread(&mesh.rm_indices[i / 3], sizeof(int32_t), 1, f);
        /* calculate normal */
        v1 = &mesh.vertices[mesh.indices[i]];
        v2 = &mesh.vertices[mesh.indices[i + 1]];
        v3 = &mesh.vertices[mesh.indices[i + 2]];
        u = vec3_sub(v2, v1);
        v = vec3_sub(v3, v1);
        mesh.normals[i / 3] = vec3_cross(&u, &v);
    }

    fclose(f);
    return mesh;
}

/* ==== RT ==== */

/** Compute Moeleer-Trumbore intersection algorithm.
 * 
 * \param ray the ray
 * \param mesh the mesh
 * \param t output distance to the hit point. If no hit, t is not modified.
 * \param ind output index of the hit triangle. If no hit, i is not modified.
 */
void moeller_trumbore(Ray *ray, Mesh *mesh, float *t, size_t *ind)
{
    /* TODO BVH */
    /* for each triangle */
    Vec3 v1, v2, v3, e1, e2, re2_cross, s, se1_cross;
    float d, u, v;
    float dist = 1e9;
    float dist_tmp;
    for (size_t i = 0; i < mesh->num_indices; i += 3) {
        v1 = mesh->vertices[mesh->indices[i]];
        v2 = mesh->vertices[mesh->indices[i + 1]];
        v3 = mesh->vertices[mesh->indices[i + 2]];
        e1 = vec3_sub(&v2, &v1);
        e2 = vec3_sub(&v3, &v1);
        re2_cross = vec3_cross(&ray->d, &e2);
        d = vec3_dot(&e1, &re2_cross);
        if (d > -EPS && d < EPS) continue;
        s = vec3_sub(&ray->o, &v1);
        u = vec3_dot(&s, &re2_cross) / d;
        if ((u < 0. && fabs(u) > EPS)
        || (u > 1. && fabs(u - 1.) > EPS))
            continue;
        se1_cross = vec3_cross(&s, &e1);
        v = vec3_dot(&ray->d, &se1_cross) / d;
        if ((v < 0. && fabs(v) > EPS)
        || (u + v > 1. && fabs(u + v - 1.) > EPS))
            continue;
        dist_tmp = vec3_dot(&e2, &se1_cross) / d;
        if (dist_tmp > EPS && dist_tmp < dist) {
            dist = dist_tmp;
            *t = dist;
            *ind = i;
        }
    }
}

/* ==== MAIN FUNCTION ==== */

void compute_paths(
    IN const char *mesh_filepath, /* path to the mesh file */
    IN const float *rx_positions,  /* shape (num_rx, 3) */
    IN const float *tx_positions,  /* shape (num_tx, 3) */
    IN const float *rx_velocities, /* shape (num_rx, 3) */
    IN const float *tx_velocities, /* shape (num_tx, 3) */
    IN float carrier_frequency,    /* > 0.0 */
    IN int num_rx,                 /* number of receivers */
    IN int num_tx,                 /* number of transmitters */
    IN int num_paths,              /* number of paths */
    IN int num_bounces,            /* number of bounces */
    OUT float *a,                  /* output array of gains (num_paths,) */
    OUT float *tau               /* output array of delays (num_paths,) */
)
{
    /* Load the scene */
    Mesh mesh = load_mesh_ply(mesh_filepath);

    /* Calculate fibonacci sphere */
    Vec3 *ray_directions = (Vec3*)malloc(num_paths * sizeof(Vec3));
    float k, phi, theta;
    for (size_t i = 0; i < num_paths; ++i) {
        k = (float)i + .5;
        phi = 1. - 2. * k / num_paths;
        theta = PI * (1. + sqrt(5.)) * k;
        ray_directions[i] = (Vec3){
            cos(theta) * sin(phi),
            sin(theta) * sin(phi),
            cos(phi)
        };
    }

    /* Create num_path rays for each tx */
    /* Shape (num_tx, num_paths) */
    Ray **rays = (Ray**)malloc(num_tx * sizeof(Ray*));
    Vec3 tx_position;
    for (size_t i = 0; i < num_tx; ++i) {
        rays[i] = (Ray*)malloc(num_paths * sizeof(Ray));
        tx_position = (Vec3){tx_positions[i * 3], tx_positions[i * 3 + 1], tx_positions[i * 3 + 2]};
        for (size_t j = 0; j < num_paths; ++j) {
            rays[i][j].o = tx_position;
            rays[i][j].d = ray_directions[j];
        }
    }
    free(ray_directions);

    /* Calculate the paths */
    /* TODO calculate a and tau in moeller_trumbore */
    /* shape (num_bounces, num_tx, num_paths) */
    Vec3 ***hits = (Vec3***)malloc(num_bounces * sizeof(Vec3**));
    size_t ***hit_indices = (size_t***)malloc(num_bounces * sizeof(size_t**));
    float t;
    size_t ind;
    for (size_t i = 0; i < num_bounces; ++i) {
        hits[i] = (Vec3**)malloc(num_tx * sizeof(Vec3*));
        hit_indices[i] = (size_t**)malloc(num_tx * sizeof(size_t*));
        for (size_t j = 0; j < num_tx; ++j) {
            hits[i][j] = (Vec3*)malloc(num_paths * sizeof(Vec3));
            hit_indices[i][j] = (size_t*)malloc(num_paths * sizeof(size_t));
            for (size_t k = 0; k < num_paths; ++k) {
                t = -1.;
                ind = -1;
                moeller_trumbore(&rays[j][k], &mesh, &t, &ind);
                hits[i][j][k] = vec3_scale(&rays[j][k].d, t);
                hits[i][j][k] = vec3_add(&hits[i][j][k], &rays[j][k].o);
                rays[j][k].o = hits[i][j][k];
                hit_indices[i][j][k] = ind;
            }
        }
    }

    /* TODO Remove */
    a = hits;
    tau = hit_indices;

    /* Free */
    /* Free the hits */
    /*
    for (size_t i = 0; i < num_bounces; ++i) {
        for (size_t j = 0; j < num_tx; ++j) {
            free(hits[i][j]);
            free(hit_indices[i][j]);
        }
        free(hits[i]);
        free(hit_indices[i]);
    }
    free(hits);
    free(hit_indices);
    */
    /* Free the rays */
    for (size_t i = 0; i < num_tx; ++i)
        free(rays[i]);
    free(rays);
    /* Free the mesh */
    free(mesh.vertices);
    free(mesh.rms);
    free(mesh.indices);
    free(mesh.rm_indices);
    free(mesh.normals);
}
