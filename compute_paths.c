#include "compute_paths.h"

#include <stddef.h> /* for size_t */
#include <stdlib.h> /* for exit, malloc, free */
#include <string.h> /* for strlen, sprintf */
#include <stdio.h>  /* for fopen, FILE, fclose */
#include <math.h>   /* for sin, cos, sqrt */
#include <stdint.h> /* for int32_t, uint8_t */

#define PI 3.14159265358979323846f /* pi */
#define SPEED_OF_LIGHT 299792458.0f /* m/s */

typedef struct {
  float x, y, z;
} Vec3;

typedef struct {
  Vec3 o, d;
} Ray;

typedef struct {
  float eta_re, eta_sqrt_re, eta_inv_re, eta_inv_sqrt_re;
  float eta_im, eta_sqrt_im, eta_inv_im, eta_inv_sqrt_im;
  float eta_abs, eta_abs_pow2, eta_abs_inv_sqrt;
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
void free_mesh(Mesh *mesh)
{
  free(mesh->vertices);
  free(mesh->rms);
  free(mesh->indices);
  free(mesh->rm_indices);
  free(mesh->normals);
}

/* ==== COMPLEX OPERATIONS ==== */

void csqrtf(
  IN float re,
  IN float im,
  IN float abs,
  OUT float *sqrt_re,
  OUT float *sqrt_im
)
{
  *sqrt_re = sqrtf((re + abs) / 2.f);
  if (fabsf(im) < __FLT_EPSILON__ && re >= -__FLT_EPSILON__)
    *sqrt_im = 0;
  else {
    *sqrt_im = sqrtf((abs - re) / 2.f);
    if (im < 0.f) *sqrt_im = -*sqrt_im;
  }
}
void cdivf(
  IN float a_re,
  IN float a_im,
  IN float b_re,
  IN float b_im,
  OUT float *c_re,
  OUT float *c_im
)
{
  float d = b_re * b_re + b_im * b_im;
  *c_re = (a_re * b_re + a_im * b_im) / d;
  *c_im = (a_im * b_re - a_re * b_im) / d;
}

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
Vec3 vec3_normalize(const Vec3 *a)
{
  float norm = sqrtf(a->x * a->x + a->y * a->y + a->z * a->z);
  return (Vec3){a->x / norm, a->y / norm, a->z / norm};
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
 * \param carrier_frequency the carrier frequency in GHz
 * \return the loaded mesh
 */
Mesh load_mesh_ply(const char *mesh_filepath, float carrier_frequency)
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
  if (fread(buff, 1, 4, f) != 4) exit(8);
  if (strncmp(buff, "ply\n", 4)) exit(8);
  /* format binary_little_endian 1.0 */
  if (fread(buff, 1, 32, f) != 32) exit(8);
  if (strncmp(buff, "format binary_little_endian 1.0\n", 32)) exit(8);

  /* element vertex <num_vertices> */
  if (fgets(buff, 128, f) == NULL) exit(8);
  if (strncmp(buff, "element vertex ", 15)) exit(8);
  sscanf(buff, "element vertex %zu", &num_vertices);
  /* property float x */
  /* property float y */
  /* property float z */
  if (fread(buff, 1, 17*3, f) != 17*3) exit(8);
  if (strncmp(buff, "property float x\nproperty float y\nproperty float z\n", 17*3))
    exit(8);

  /* element radio_material */
  if (fgets(buff, 128, f) == NULL) exit(8);
  if (strncmp(buff, "element radio_material ", 23)) exit(8);
  sscanf(buff, "element radio_material %zu", &num_materials);
  /* property list uint char name */
  if (fread(buff, 1, 29, f) != 29) exit(8);
  if (strncmp(buff, "property list uint char name\n", 29)) exit(8);
  /* property float a */
  /* property float b */
  /* property float c */
  /* property float d */
  if (fread(buff, 1, 17*4, f) != 17*4) exit(8);
  if (strncmp(buff, "property float a\nproperty float b\nproperty float c\nproperty float d\n", 17*4))
    exit(8);

  /* element face <num_faces> */
  if (fgets(buff, 128, f) == NULL) exit(8);
  if (strncmp(buff, "element face ", 13)) exit(8);
  sscanf(buff, "element face %zu", &num_faces);
  /* property list uchar int vertex_index */
  if (fread(buff, 1, 38, f) != 38) exit(8);
  if (strncmp(buff, "property list uchar uint vertex_index\n", 38))
    exit(8);
  /* property uint material_index */
  if (fread(buff, 1, 29, f) != 29) exit(8);
  if (strncmp(buff, "property uint material_index\n", 29)) exit(8);

  /* end_header */
  if (fread(buff, 1, 11, f) != 11) exit(8);
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
    if (fread(&mesh.vertices[i], sizeof(Vec3), 1, f) != 1)
      exit(8);

  /* RADIO_MATERIALS */
  /* read a, b, c and d
  and calculate relative permitivity */
  unsigned int name_size;
  float abcd[4];
  RadioMaterial *rm;
  for (size_t i = 0; i < num_materials; ++i) {
    /* skip name */
    if (fread(&name_size, sizeof(unsigned int), 1, f) != 1) exit(8);
    fseek(f, name_size, SEEK_CUR);
    /* read a, b, c and d */
    if (fread(abcd, sizeof(float)*4, 1, f) != 1) exit(8);
    rm = &mesh.rms[i];
    /* calculate eta */
    rm->eta_re = abcd[0] * powf(carrier_frequency, abcd[1]);
    /* eq. 12 */
    rm->eta_im = (abcd[2] * powf(carrier_frequency, abcd[3]))
               / (0.0556325027352135f * carrier_frequency);
    rm->eta_abs_pow2 = rm->eta_re * rm->eta_re + rm->eta_im * rm->eta_im;
    rm->eta_abs = sqrtf(rm->eta_abs_pow2);
    rm->eta_abs_inv_sqrt = 1.f / sqrtf(rm->eta_abs);
    /* calculate sqrt(eta) */
    csqrtf(rm->eta_re, rm->eta_im, rm->eta_abs, &rm->eta_sqrt_re, &rm->eta_sqrt_im);
    /* calculate 1 / eta */
    rm->eta_inv_re = rm->eta_re / rm->eta_abs_pow2;
    rm->eta_inv_im = -rm->eta_im / rm->eta_abs_pow2;
    /* calculate 1 / sqrt(eta) */
    csqrtf(rm->eta_inv_re, rm->eta_inv_im, 1.f / rm->eta_abs,
          &rm->eta_inv_sqrt_re, &rm->eta_inv_sqrt_im);
  }

  /* INDICES */
  /* also calculate normals */
  uint8_t n;
  Vec3 *v1, *v2, *v3, u, v;
  for (size_t i = 0; i < mesh.num_indices; i += 3) {
    if (fread(&n, 1, 1, f) != 1) exit(8);
    if (n != 3) {
      fprintf(stderr, "Only trianglar faces are supported\n");
      exit(8);
    }
    if (fread(&mesh.indices[i], sizeof(int32_t), 3, f) != 3) exit(8);
    if (fread(&mesh.rm_indices[i / 3], sizeof(int32_t), 1, f) != 1) exit(8);
    /* calculate normal */
    v1 = &mesh.vertices[mesh.indices[i]];
    v2 = &mesh.vertices[mesh.indices[i + 1]];
    v3 = &mesh.vertices[mesh.indices[i + 2]];
    u = vec3_sub(v2, v1);
    v = vec3_sub(v3, v1);
    mesh.normals[i / 3] = vec3_cross(&u, &v);
    mesh.normals[i / 3] = vec3_normalize(&mesh.normals[i / 3]);
  }

  fclose(f);
  return mesh;
}

/* ==== RT ==== */

/** Compute Moeleer-Trumbore intersection algorithm.
 * 
 * \param ray ray to cast to the mesh
 * \param mesh triangle mesh
 * \param t output distance to the hit point. If no hit, t is not modified.
 * \param ind output index of the hit triangle. If no hit, i is not modified.
 * \param theta output angle of incidence
 */
void moeller_trumbore(
  IN Ray *ray,
  IN Mesh *mesh,
  OUT float *t,
  OUT int32_t *ind,
  OUT float *theta
)
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
    if (d > -__FLT_EPSILON__ && d < __FLT_EPSILON__) continue;
    s = vec3_sub(&ray->o, &v1);
    u = vec3_dot(&s, &re2_cross) / d;
    if ((u < 0. && fabs(u) > __FLT_EPSILON__)
    || (u > 1. && fabs(u - 1.) > __FLT_EPSILON__))
      continue;
    se1_cross = vec3_cross(&s, &e1);
    v = vec3_dot(&ray->d, &se1_cross) / d;
    if ((v < 0. && fabs(v) > __FLT_EPSILON__)
    || (u + v > 1. && fabs(u + v - 1.) > __FLT_EPSILON__))
      continue;
    dist_tmp = vec3_dot(&e2, &se1_cross) / d;
    if (dist_tmp > __FLT_EPSILON__ && dist_tmp < dist) {
      dist = dist_tmp;
      *t = dist;
      *ind = i / 3;
      *theta = acos(vec3_dot(&mesh->normals[i / 3], &ray->d));
      if (*theta > PI / 2.)
        *theta = PI - *theta;
    }
  }
}

/** Compute reflection R_{eTE} and R_{eTM} coefficients.
 * 
 * Implements eqs. (31a)-(31b) from ITU-R P.2040-3.
 * 
 * \param rm the radio material of the hit triangle
 * \param theta1 the angle of incidence
 * \param r_te_re output real part of R_{eTE}
 * \param r_te_im output imaginary part of R_{eTE}
 * \param r_tm_re output real part of R_{eTM}
 * \param r_tm_im output imaginary part of R_{eTM}
 */
void refl_coefs(
  IN RadioMaterial *rm,
  IN float theta1,
  OUT float *r_te_re,
  OUT float *r_te_im,
  OUT float *r_tm_re,
  OUT float *r_tm_im
)
{
  float sin_theta1 = sinf(theta1);
  if (rm->eta_abs_inv_sqrt * sin_theta1 > 1.f - __FLT_EPSILON__) {
    *r_te_re = *r_tm_re = 1.f;
    *r_te_im = *r_tm_im = 0.f;
    return;
  }

  /* eq. 33 */
  float sin_theta1_pow2 = sin_theta1 * sin_theta1;
  float cos_theta2_re = sqrtf(1.f + rm->eta_inv_re / rm->eta_abs_pow2 * sin_theta1_pow2);
  float cos_theta2_im = sqrtf(1.f - rm->eta_inv_im / rm->eta_abs_pow2 * sin_theta1_pow2);

  /* R_{eTE}, eq. 31a */
  float sqrt_eta_cos_theta2_re = rm->eta_sqrt_re * cos_theta2_re - rm->eta_sqrt_im * cos_theta2_im;
  float sqrt_eta_cos_theta2_im = rm->eta_sqrt_re * cos_theta2_im + rm->eta_sqrt_im * cos_theta2_re;
  float cos_theta1 = cosf(theta1);
  cdivf(cos_theta1 - sqrt_eta_cos_theta2_re, -sqrt_eta_cos_theta2_im,
        cos_theta1 + sqrt_eta_cos_theta2_re, sqrt_eta_cos_theta2_im,
        r_te_re, r_te_im);

  /* R_{eTM}, eq. 31b */
  float sqrt_eta_cos_theta1_re = rm->eta_sqrt_re * cos_theta1;
  float sqrt_eta_cos_theta1_im = rm->eta_sqrt_im * cos_theta1;
  cdivf(sqrt_eta_cos_theta1_re - cos_theta2_re,
        sqrt_eta_cos_theta1_im - cos_theta2_im,
        sqrt_eta_cos_theta1_re + cos_theta2_re,
        sqrt_eta_cos_theta1_im + cos_theta2_im,
        r_tm_re, r_tm_im);
}

/* ==== MAIN FUNCTION ==== */

void compute_paths(
  IN const char *mesh_filepath,   /* path to the mesh file */
  IN const float *rx_pos,         /* shape (num_rx, 3) */
  IN const float *tx_pos,         /* shape (num_tx, 3) */
  IN const float *rx_vel,         /* shape (num_rx, 3) */
  IN const float *tx_vel,         /* shape (num_tx, 3) */
  IN float carrier_frequency,     /* > 0.0 (IN GHz!) */
  IN size_t num_rx,               /* number of receivers */
  IN size_t num_tx,               /* number of transmitters */
  IN size_t num_paths,            /* number of paths */
  IN size_t num_bounces,          /* number of bounces */
  /* LoS */
  OUT float *a_te_re_los,         /* output array real parts of TE gains (num_rx, num_tx) */
  OUT float *a_te_im_los,         /* output array imaginary parts of TE gains (num_rx, num_tx) */
  OUT float *a_tm_re_los,         /* output array real parts of TM gains (num_rx, num_tx) */
  OUT float *a_tm_im_los,         /* output array imaginary parts of TM gains (num_rx, num_tx) */
  OUT float *tau_los,             /* output array of delays (num_rx, num_tx) */
  /* Scatter */
  OUT float *a_te_re_scat,        /* output array real parts of TE gains (num_bounces, num_rx, num_tx, num_paths) */
  OUT float *a_te_im_scat,        /* output array imaginary parts of TE gains (num_bounces, num_rx, num_tx, num_paths) */
  OUT float *a_tm_re_scat,        /* output array real parts of TM gains (num_bounces, num_rx, num_tx, num_paths) */
  OUT float *a_tm_im_scat,        /* output array imaginary parts of TM gains (num_bounces, num_rx, num_tx, num_paths) */
  OUT float *tau_scat             /* output array of delays (num_bounces, num_rx, num_tx, num_paths) */
)
{
  /* Init loop indices and offset variables */
  size_t i, j;
  size_t off, off_scat;

  /* Load the scene */
  Mesh mesh = load_mesh_ply(mesh_filepath, carrier_frequency);

  /* Calculate a fibonacci sphere */
  Vec3 *ray_directions = (Vec3*)malloc(num_paths * sizeof(Vec3));
  float k, phi, theta;
  for (i = 0; i < num_paths; ++i) {
    k = (float)i + .5f;
    phi = acos(1.f - 2.f * k / num_paths);
    theta = PI * (1.f + sqrtf(5.f)) * k;
    ray_directions[i] = (Vec3){
      cos(theta) * sin(phi),
      sin(theta) * sin(phi),
      cos(phi)
    };
  }

  /* Cast the positions and velocities to Vec3 */
  Vec3 *rx_pos_v = (Vec3*)rx_pos;
  Vec3 *tx_pos_v = (Vec3*)tx_pos;
  Vec3 *rx_vel_v = (Vec3*)rx_vel;
  Vec3 *tx_vel_v = (Vec3*)tx_vel;

  /* Create num_path rays for each tx */
  /* Shape (num_tx, num_paths) */
  Ray *rays = (Ray*)malloc(num_tx * num_paths * sizeof(Ray));
  off = 0;
  for (i = 0; i < num_tx; ++i) {
    for (j = 0; j < num_paths; ++j) {
      rays[off].o = tx_pos_v[i];
      rays[off].d = ray_directions[j];
      ++off;
    }
  }
  free(ray_directions);

  /* Initialize variables needed for the RT */

  /* Bouncing (reflecting) rays gains and delays */
  /* shape (num_tx, num_paths) */
  float *a_te_re_refl = (float*)malloc(num_tx * num_paths * sizeof(float));
  float *a_te_im_refl = (float*)calloc(num_tx * num_paths, sizeof(float));
  float *a_tm_re_refl = (float*)malloc(num_tx * num_paths * sizeof(float));
  float *a_tm_im_refl = (float*)calloc(num_tx * num_paths, sizeof(float));
  for (size_t i = 0; i < num_tx * num_paths; ++i)
    a_te_re_refl[i] = a_tm_re_refl[i] = 1.f;
  float *tau_t = (float*)calloc(num_tx * num_paths, sizeof(float));

  /* Active rays bitmask. 1 if the ray is still traced, 0 if it has left the scene */
  uint8_t *active = (uint8_t*)malloc((num_tx * num_paths / 8 + 1) * sizeof(uint8_t));
  for (size_t i = 0; i < num_tx * num_paths / 8 + 1; ++i)
    active[i] = 0xff;
  size_t num_active = num_tx * num_paths;

  /* Temporary variables */
  float t, r_te_re, r_te_im, r_tm_re, r_tm_im;
  float a_te_re_new, a_te_im_new, a_tm_re_new, a_tm_im_new;
  int32_t ind;
  Ray r, *r_ptr;
  Vec3 *h = (Vec3*)malloc(sizeof(Vec3));

  /* Friis free space loss multiplier.
  Multiply this by a distance and take the second power to get a free space loss. */
  float free_space_loss_multiplier = 2.f * PI * carrier_frequency * 1e9 / SPEED_OF_LIGHT;
  float free_space_loss;

  /* Calculate LoS paths directly from tx to rx */
  for (i = 0; i != num_rx; ++i)
    for (j = 0; j != num_tx; ++j) {
      off = i * num_tx + j;
      t = -1.f;
      ind = -1;
      r.o = tx_pos_v[j];
      r.d = vec3_sub(&rx_pos_v[i], &r.o);
      moeller_trumbore(&r, &mesh, &t, &ind, &theta);
      if (ind != -1 && t <= 1.f) {
        /* An obstacle between the tx and the rx has been hit */
        a_te_re_los[off] = a_te_im_los[off] = a_tm_re_los[off] = a_tm_im_los[off] = tau_los[off] = 0.f;
        continue;
      }
      /* No triangle has been hit, the path is LoS */
      /* Calculate the distance */
      t = sqrtf(vec3_dot(&r.d, &r.d));
      /* Calculate the free space loss */
      free_space_loss = free_space_loss_multiplier * t;
      free_space_loss = free_space_loss * free_space_loss;
      /* Set the gains for this LoS path */
      a_te_re_los[off] = a_tm_re_los[off] = 1.f / free_space_loss;
      a_te_im_los[off] = a_tm_im_los[off] = 0.f;
      /* Calculate the delay */
      tau_los[off] = t / SPEED_OF_LIGHT;
    }

  /* Bounce the rays. On each bounce:
  * - Add path length / SPEED_OF_LIGHT to tau
  * - Update gains using eqs. (31a)-(31b) from ITU-R P.2040-3
  * - Spawn refraction rays with gains as per eqs. (31c)-(31d) from ITU-R P.2040-3
  */
  for (i = 0; i < num_bounces; ++i) {
    for (off = 0; off != num_tx * num_paths; ++off) {
      /* Check if the ray is active */
      if (!(active[off / 8] & (1 << (off % 8))))
        continue;
      /* Init */
      t = -1.f;
      ind = -1;
      r_ptr = &rays[off];
      /* Find the hit point and trinagle.
      Calculate an angle of incidence */
      moeller_trumbore(r_ptr, &mesh, &t, &ind, &theta);
      if (ind == -1) { /* Ray hit nothing */
        active[off / 8] &= ~(1 << (off % 8));
        --num_active;
        continue;
      }
      /* Calculate the reflection coefficients
      R_{eTE} and R_{eTM} */
      refl_coefs(&mesh.rms[mesh.rm_indices[ind]],
                  theta,
                  &r_te_re, &r_te_im,
                  &r_tm_re, &r_tm_im);
      /* Calculate the free space loss */
      free_space_loss = free_space_loss_multiplier * t;
      free_space_loss = free_space_loss * free_space_loss;
      r_te_re /= free_space_loss;
      r_te_im /= free_space_loss;
      r_tm_re /= free_space_loss;
      r_tm_im /= free_space_loss;
      /* Update the gains */
      a_te_re_new = a_te_re_refl[off] * r_te_re - a_te_im_refl[off] * r_te_im;
      a_te_im_new = a_te_re_refl[off] * r_te_im + a_te_im_refl[off] * r_te_re;
      a_tm_re_new = a_tm_re_refl[off] * r_tm_re - a_tm_im_refl[off] * r_tm_im;
      a_tm_im_new = a_tm_re_refl[off] * r_tm_im + a_tm_im_refl[off] * r_tm_re;
      a_te_re_refl[off] = a_te_re_new;
      a_te_im_refl[off] = a_te_im_new;
      a_tm_re_refl[off] = a_tm_re_new;
      a_tm_im_refl[off] = a_tm_im_new;
      /* Update the delay */
      tau_t[off] += t / SPEED_OF_LIGHT;
      /* Advance the ray to the hit point */
      *h = vec3_scale(&r_ptr->d, t);
      r_ptr->o = vec3_add(h, &r_ptr->o);
      /* Scatter a path from the hit point towards the rx */
      r.o = r_ptr->o;
      for (j = 0; j != num_rx; ++j) {
        off_scat = i * num_rx * num_tx * num_paths + j * num_tx * num_paths + off;
        r.d = vec3_sub(&rx_pos_v[j], &r.o);
        ind = -1;
        moeller_trumbore(&r, &mesh, &t, &ind, &theta);
        if (ind != -1 && t <= 1.f) {
          /* An obstacle between the hit point and the rx has been hit */
          a_te_re_scat[off_scat] = a_te_im_scat[off_scat] = a_tm_re_scat[off_scat]
                                 = a_tm_im_scat[off_scat] = tau_scat[off_scat] = 0.f;
          continue;
        }
        /* No obstacle has been hit */
        /* Calculate the distance */
        t = sqrtf(vec3_dot(&r.d, &r.d));
        /* Calculate the delay */
        tau_scat[off_scat] = tau_t[off] + t / SPEED_OF_LIGHT;
        /* Calculate the free space loss */
        free_space_loss = free_space_loss_multiplier * t;
        free_space_loss = free_space_loss * free_space_loss;
        /* Set the gains */
        a_te_re_scat[off_scat] = a_te_re_refl[off] / free_space_loss;
        a_te_im_scat[off_scat] = a_te_im_refl[off] / free_space_loss;
        a_tm_re_scat[off_scat] = a_tm_re_refl[off] / free_space_loss;
        a_tm_im_scat[off_scat] = a_tm_im_refl[off] / free_space_loss;
      }
    }
  }

  /* Free */
  free(a_te_re_refl);
  free(a_te_im_refl);
  free(a_tm_re_refl);
  free(a_tm_im_refl);
  free(tau_t);
  free(h);
  free(active);
  free(rays);
  free_mesh(&mesh);
}
