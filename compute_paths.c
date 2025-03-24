#include "compute_paths.h" /* for compute_paths */
#include "scene.h" /* for HRT_Scene, HRT_mesh, HRT_Material */

#include <stddef.h> /* for size_t */
#include <stdlib.h> /* for exit, malloc, free */
#include <string.h> /* for strlen, sprintf */
#include <stdio.h>  /* for fopen, FILE, fclose */
#include <math.h>   /* for sin, cos, sqrt */
#include <stdint.h> /* for int32_t, uint8_t */

/* ==== CONSTANTS ==== */

#define PI 3.14159265358979323846f /* pi */
#define SPEED_OF_LIGHT 299792458.0f /* m/s */

/* (2w, w) binomial coefficients for w = 0, 1, ..., 19 */
const float BINOMIAL_2W_W[20] = {
  1.f, 2.f, 6.f, 20.f, 70.f,
  252.f, 924.f, 3432.f, 12870.f, 48620.f,
  184756.f, 705432.f, 2704156.f, 10400600.f, 40116600.f,
  155117520.f, 601080390.f, 2333606220.f, 9075135300.f, 35345263800.f
};
/* (aplha, k) binomial coefficients for alpha = 0, 1, ..., 19 with k = 0, 1, ..., alpha */
const float* BINOMIAL_ALPHA_K[20] = {
  (float[]){
    1.f
  },
  (float[]){
    1.f, 1.f
  },
  (float[]){
    1.f, 2.f, 1.f
  },
  (float[]){
    1.f, 3.f, 3.f, 1.f
  },
  (float[]){
    1.f, 4.f, 6.f, 4.f, 1.f
  },
  (float[]){
    1.f, 5.f, 10.f, 10.f, 5.f,
    1.f
  },
  (float[]){
    1.f, 6.f, 15.f, 20.f, 15.f,
    6.f, 1.f
  },
  (float[]){
    1.f, 7.f, 21.f, 35.f, 35.f,
    21.f, 7.f, 1.f
  },
  (float[]){
    1.f, 8.f, 28.f, 56.f, 70.f,
    56.f, 28.f, 8.f, 1.f
  },
  (float[]){
    1.f, 9.f, 36.f, 84.f, 126.f,
    126.f, 84.f, 36.f, 9.f, 1.f
  },
  (float[]){
    1.f, 10.f, 45.f, 120.f, 210.f,
    252.f, 210.f, 120.f, 45.f, 10.f,
    1.f
  },
  (float[]){
    1.f, 11.f, 55.f, 165.f, 330.f,
    462.f, 462.f, 330.f, 165.f, 55.f,
    11.f, 1.f
  },
  (float[]){
    1.f, 12.f, 66.f, 220.f, 495.f,
    792.f, 924.f, 792.f, 495.f, 220.f,
    66.f, 12.f, 1.f
  },
  (float[]){
    1.f, 13.f, 78.f, 286.f,
    715.f, 1287.f, 1716.f, 1716.f, 1287.f,
    715.f, 286.f, 78.f, 13.f, 1.f
  },
  (float[]){
    1.f, 14.f, 91.f, 364.f, 1001.f,
    2002.f, 3003.f, 3432.f, 3003.f, 2002.f,
    1001.f, 364.f, 91.f, 14.f, 1.f
  },
  (float[]){
    1.f, 15.f, 105.f, 455.f, 1365.f,
    3003.f, 5005.f, 6435.f, 6435.f, 5005.f,
    3003.f, 1365.f, 455.f, 105.f, 15.f,
    1.f
  },
  (float[]){
    1.f, 16.f, 120.f, 560.f, 1820.f,
    4368.f, 8008.f, 11440.f, 12870.f, 11440.f,
    8008.f, 4368.f, 1820.f, 560.f, 120.f,
    16.f, 1.f
  },
  (float[]){
    1.f, 17.f, 136.f, 680.f, 2380.f,
    6188.f, 12376.f, 19448.f, 24310.f, 24310.f,
    19448.f, 12376.f, 6188.f, 2380.f, 680.f,
    136.f, 17.f, 1.f
  },
  (float[]){
    1.f, 18.f, 153.f, 816.f, 3060.f,
    8568.f, 18564.f, 31824.f, 43758.f, 48620.f,
    43758.f, 31824.f, 18564.f, 8568.f, 3060.f,
    816.f, 153.f, 18.f, 1.f
  },
  (float[]){
    1.f, 19.f, 171.f, 969.f, 3876.f,
    11628.f, 27132.f, 50388.f, 75582.f, 92378.f,
    92378.f, 75582.f, 50388.f, 27132.f, 11628.f,
    3876.f, 969.f, 171.f, 19.f, 1.f
  }
};

/* ==== STRUCTS ==== */

typedef struct {
  float x, y, z;
} Vec3;

typedef struct {
  Vec3 o, d;
} Ray;

/* Radio material eta and additional parameters */
typedef struct {
  /* eta */
  float eta_re, eta_sqrt_re, eta_inv_re, eta_inv_sqrt_re;
  float eta_im, eta_sqrt_im, eta_inv_im, eta_inv_sqrt_im;
  float eta_abs, eta_abs_pow2, eta_abs_inv_sqrt;
  /* r = 1.0 - s */
  float r;
} Material;

typedef struct {
  uint32_t num_vertices;
  Vec3 *vs;
  uint32_t num_triangles;
  uint32_t *is;
  /* Normals. Size [num_triangles] */
  Vec3 *ns;
  uint32_t material_index;
  float velocity[3];
} Mesh;

typedef struct {
  uint32_t num_meshes;
  Mesh *meshes;
} Scene;

/* ==== PRECOMPUTED GLOBALS ==== */

/* 4.f * PI * carrier_frequency * 1e9 / SPEED_OF_LIGHT */
float g_free_space_loss_multiplier;

/* Materials etas, calucalted for the given carrier_frequency */
Material g_materials[NUM_G_MATERIALS];

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

/* ==== SCENE LOADING ==== */

/** Load a mesh from a HRT file.
 * 
 * \param scene_filepath path to the HRT file
 * \param carrier_frequency the carrier frequency in GHz
 * \return the loaded scene
 */
Scene load_scene(const char *scene_filepath, float carrier_frequency)
{
  /* Open the HRT file */
  FILE *f = fopen(scene_filepath, "rb");
  if (!f) {
    fprintf(stderr, "Could not open file %s\n", scene_filepath);
    exit(8);
  }

  /* Parse the file */
  /* MAGIC */
  char magic[3];
  if (fread(magic, 1, 3, f) != 3) exit(8);
  if (strncmp(magic, "HRT", 3)) exit(8);
  /* SCENE */
  Scene scene;
  /* num_meshes */
  if (fread(&scene.num_meshes, sizeof(uint32_t), 1, f) != 1) exit(8);
  /* meshes */
  scene.meshes = (Mesh*)malloc(scene.num_meshes * sizeof(Mesh));
  for (uint32_t i = 0; i != scene.num_meshes; ++i) {
    /* MESH */
    Mesh *mesh = &scene.meshes[i];
    /* num_vertices */
    if(fread(&mesh->num_vertices, sizeof(uint32_t), 1, f) != 1) exit(8);
    /* vertices */
    mesh->vs = (Vec3*)malloc(mesh->num_vertices * sizeof(Vec3));
    if (fread(mesh->vs, sizeof(Vec3), mesh->num_vertices, f) != mesh->num_vertices)
      exit(8);
    /* num_triangles */
    if (fread(&mesh->num_triangles, sizeof(uint32_t), 1, f) != 1) exit(8);
    /* triangles */
    mesh->is = (uint32_t*)malloc(mesh->num_triangles * 3 * sizeof(uint32_t));
    if (fread(mesh->is, sizeof(uint32_t), mesh->num_triangles * 3, f) != mesh->num_triangles * 3)
      exit(8);
    /* material_index */
    if (fread(&mesh->material_index, sizeof(uint32_t), 1, f) != 1) exit(8);
    /* velocity */
    if (fread(&mesh->velocity, sizeof(float), 3, f) != 3) exit(8);
  }
  fclose(f);

  /* Precompute material permitivity */
  uint8_t material_eta_isdone[NUM_G_MATERIALS] = {0};
  for (uint32_t i = 0; i != scene.num_meshes; ++i)
    if (!material_eta_isdone[scene.meshes[i].material_index]) {
      HRT_Material *hrt_mat = &g_hrt_materials[scene.meshes[i].material_index];
      Material *rm = &g_materials[scene.meshes[i].material_index];
      /* calculate eta */
      rm->eta_re = hrt_mat->a * powf(carrier_frequency, hrt_mat->b);
      /* eq. 12 */
      rm->eta_im = (hrt_mat->c * powf(carrier_frequency, hrt_mat->d))
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
      /* calculate r */
      rm->r = 1.f - hrt_mat->s;
    }

  /* Calculate normals */
  for (uint32_t i = 0; i != scene.num_meshes; ++i) {
    Mesh *mesh = &scene.meshes[i];
    mesh->ns = (Vec3*)malloc(mesh->num_triangles * sizeof(Vec3));
    for (uint32_t j = 0; j != mesh->num_triangles; ++j) {
      Vec3 v1 = mesh->vs[mesh->is[j * 3]];
      Vec3 v2 = mesh->vs[mesh->is[j * 3 + 1]];
      Vec3 v3 = mesh->vs[mesh->is[j * 3 + 2]];
      Vec3 e1 = vec3_sub(&v2, &v1);
      Vec3 e2 = vec3_sub(&v3, &v1);
      mesh->ns[j] = vec3_cross(&e1, &e2);
      mesh->ns[j] = vec3_normalize(&mesh->ns[j]);
    }
  }

  return scene;
}

/* ==== RT ==== */

/** Compute Moeleer-Trumbore intersection algorithm.
 * 
 * \param ray ray to cast to the mesh
 * \param scene the scene to cast the ray to
 * \param t output distance to the hit point. If no hit, t is not modified.
 * \param mesh_ind output index of the hit mesh. If no hit, i is not modified.
 * \param face_ind output index of the hit triangle. If no hit, i is not modified.
 * \param theta output angle of incidence
 */
void moeller_trumbore(
  IN Ray *ray,
  IN Scene *scene,
  OUT float *t,
  OUT uint32_t *mesh_ind,
  OUT uint32_t *face_ind,
  OUT float *theta
)
{
  /* TODO BVH */
  /* for each triangle */
  Vec3 v1, v2, v3, e1, e2, re2_cross, s, se1_cross;
  float d, u, v;
  float dist = 1e9;
  float dist_tmp;
  for (uint32_t i = 0; i != scene->num_meshes; ++i) {
    Mesh *mesh = &scene->meshes[i];
    for (uint32_t j = 0; j != mesh->num_triangles; ++j) {
      v1 = mesh->vs[mesh->is[3 * j]];
      v2 = mesh->vs[mesh->is[3 * j + 1]];
      v3 = mesh->vs[mesh->is[3 * j + 2]];
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
        *mesh_ind = i;
        *face_ind = j;
        *theta = acos(vec3_dot(&mesh->ns[j], &ray->d));
        if (*theta > PI / 2.)
          *theta = PI - *theta;
      }
    }
  }
}

/** Compute reflection R_{eTE} and R_{eTM} coefficients.
 * 
 * Implements eqs. (31a)-(31b) from ITU-R P.2040-3.
 * 
 * \param material_index the index of the material
 * \param theta1 the angle of incidence
 * \param r_te_re output real part of R_{eTE}
 * \param r_te_im output imaginary part of R_{eTE}
 * \param r_tm_re output real part of R_{eTM}
 * \param r_tm_im output imaginary part of R_{eTM}
 */
void refl_coefs(
  IN uint32_t material_index,
  IN float theta1,
  OUT float *r_te_re,
  OUT float *r_te_im,
  OUT float *r_tm_re,
  OUT float *r_tm_im
)
{
  Material *rm = &g_materials[material_index];
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

  /* Apply reflection reduction factor */
  *r_te_re *= rm->r;
  *r_te_im *= rm->r;
  *r_tm_re *= rm->r;
  *r_tm_im *= rm->r;
}

/** Calculate scattering coefficients.
 * 
 * Implements a directive scattering model inspired by Blaunstein et al. (DOI 10.1109/TAP.2006.888422).
 * Models scattering from rough surfaces (e.g., vegetation) with polarization dependence.
 * 
 * \param theta_s scattering angle (radians, between scattered direction and surface normal)
 * \param theta_i angle of incidence (radians, between incident ray and surface normal)
 * \param material_index index into g_hrt_materials array
 * \param a_te_re output real part of TE scattering coefficient
 * \param a_te_im output imaginary part of TE scattering coefficient
 * \param a_tm_re output real part of TM scattering coefficient
 * \param a_tm_im output imaginary part of TM scattering coefficient
 */
void scat_coefs(
  IN float theta_s,
  IN float theta_i,
  IN uint32_t material_index,
  OUT float *a_te_re,
  OUT float *a_te_im,
  OUT float *a_tm_re,
  OUT float *a_tm_im
)
{
  HRT_Material *mat = &g_hrt_materials[material_index];
  
  /* Precompute trigonometric values */
  float cos_theta_s = cosf(theta_s);
  float cos_theta_i = cosf(theta_i);
  float sin_theta_i = sinf(theta_i);
  
  /* Scattering directivity pattern: exponential decay with angle difference */
  /* f = s * exp(-alpha * |theta_s - theta_i|) models directive scattering */
  float theta_diff = fabsf(theta_s - theta_i);
  float f = mat->s * expf(-mat->s1_alpha * theta_diff);
  
  /* Roughness factor: reduces specular component as alpha increases */
  float roughness = 1.0f / (1.0f + mat->s1_alpha);  // Simplified, 0 to 1
  float specular = roughness * cos_theta_s;        // Specular-like term
  float diffuse = (1.0f - roughness) * cos_theta_s; // Diffuse term
  
  /* TE and TM coefficients (real parts) */
  /* TE: perpendicular to plane of incidence, less affected by angle */
  float te_real = f * (specular + diffuse);
  /* TM: parallel to plane of incidence, stronger angular dependence */
  float tm_real = f * (specular * cos_theta_i + diffuse);
  
  /* Phase shift (imaginary parts) due to surface roughness or path difference */
  /* Simplified: assume small phase shift proportional to roughness and angle */
  float phase_shift = mat->s1_alpha * sin_theta_i * 0.1f;  // Arbitrary scaling
  float te_imag = te_real * sinf(phase_shift);
  float tm_imag = tm_real * sinf(phase_shift);
  
  /* Normalize to ensure energy conservation (optional, adjust as needed) */
  float norm = sqrtf(te_real * te_real + te_imag * te_imag + 
                     tm_real * tm_real + tm_imag * tm_imag);
  if (norm > 1e-6f) {  // Avoid division by zero
    te_real /= norm;
    te_imag /= norm;
    tm_real /= norm;
    tm_imag /= norm;
  }
  
  /* Output coefficients */
  *a_te_re = te_real;
  *a_te_im = te_imag;
  *a_tm_re = tm_real;
  *a_tm_im = tm_imag;
  
  /* TODO: Incorporate s2 and s3 for additional roughness or frequency effects */
}

/* ==== MAIN FUNCTION ==== */

void compute_paths(
  IN const char *scene_filepath,  /* path to the scene file */
  IN const float *rx_pos,         /* shape (num_rx, 3) */
  IN const float *tx_pos,         /* shape (num_tx, 3) */
  IN const float *rx_vel,         /* shape (num_rx, 3) */
  IN const float *tx_vel,         /* shape (num_tx, 3) */
  IN float carrier_frequency,     /* > 0.0 (IN GHz!) */
  IN size_t num_rx,               /* number of receivers */
  IN size_t num_tx,               /* number of transmitters */
  IN size_t num_paths,            /* number of paths */
  IN size_t num_bounces,          /* number of bounces */
  OUT PathsInfo *los,             /* output LoS information. los->num_paths = 1 */
  OUT PathsInfo *scatter          /* output scatter information. scatter->num_paths = num_bounces*num_paths */
)
{
  /* Load the scene */
  Scene scene = load_scene(scene_filepath, carrier_frequency);

  /* Precompute globals */
  /* g_free_space_loss_multiplier */
  g_free_space_loss_multiplier = 4.f * PI * carrier_frequency * 1e9 / SPEED_OF_LIGHT;

  /* Calculate a fibonacci sphere */
  Vec3 *ray_directions = (Vec3*)malloc(num_paths * sizeof(Vec3));
  float k, phi, theta;
  for (size_t path = 0; path < num_paths; ++path) {
    k = (float)path + .5f;
    phi = acos(1.f - 2.f * k / num_paths);
    theta = PI * (1.f + sqrtf(5.f)) * k;
    ray_directions[path] = (Vec3){
      cos(theta) * sin(phi),
      sin(theta) * sin(phi),
      cos(phi)
    };
  }

  /* Record the tx directions */
  /* This is only valid for the fibbonaci sphere */
  scatter->directions_tx = (float*)memcpy(
    scatter->directions_tx,
    ray_directions,
    num_paths * sizeof(Vec3)
  );
  if (scatter->directions_tx == NULL) {
    fprintf(stderr, "Failed to copy memory to scatter->directions_tx\n");
    exit(70);
  }

  /* Cast the positions and velocities to Vec3 */
  Vec3 *rx_pos_v = (Vec3*)rx_pos;
  Vec3 *tx_pos_v = (Vec3*)tx_pos;
  Vec3 *rx_vel_v = (Vec3*)rx_vel;
  Vec3 *tx_vel_v = (Vec3*)tx_vel;

  /* Create num_path rays for each tx */
  /* Shape (num_tx, num_paths) */
  Ray *rays = (Ray*)malloc(num_tx * num_paths * sizeof(Ray));
  for (size_t tx = 0, off = 0; tx < num_tx; ++tx) {
    for (size_t path = 0; path < num_paths; ++path, ++off) {
      rays[off].o = tx_pos_v[tx];
      rays[off].d = ray_directions[path];
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
  float theta_s; /* scattering angle */
  uint32_t mesh_ind, face_ind, material_index;
  Ray r;
  Vec3 *h = (Vec3*)malloc(sizeof(Vec3));
  Vec3 n;

  /* Friis free space loss multiplier.
  Multiply this by a distance and take the second power to get a free space loss. */
  float free_space_loss_multiplier = 4.f * PI * carrier_frequency * 1e9 / SPEED_OF_LIGHT;
  float free_space_loss;

  /* Doppler frequency shift carrier frequency multiplier */
  float doppler_shift_multiplier = carrier_frequency * 1e9 / SPEED_OF_LIGHT;

  /* Initialize the doppler frequency shift using the tx velocities */
  /* The scatter->freq_shift has the shape of (num_rx, num_tx, num_bounces*num_paths) */
  /* The first num_tx*num_paths elements are the same for each rx and bounce */
  /* so we can initialize them only once and memcpy them to the rest */
  for (size_t tx = 0; tx != num_tx; ++tx)
    for (size_t path = 0; path != num_paths; ++path) {
      size_t off_rays = tx * num_paths + path;
      size_t off_scat = tx * num_paths * num_bounces + path;
      scatter->freq_shift[off_scat] = vec3_dot(&tx_vel_v[tx], &rays[off_rays].d);
      scatter->freq_shift[off_scat] *= doppler_shift_multiplier;
    }
  for (size_t bounce = 1; bounce < num_bounces; ++bounce)
    if (!memcpy(scatter->freq_shift + num_tx * num_paths * bounce,
                scatter->freq_shift, num_tx * num_paths * sizeof(float)))
      exit(70);
  for (size_t rx = 1; rx < num_rx; ++rx)
    if (!memcpy(scatter->freq_shift + num_tx * num_paths * num_bounces * rx,
                scatter->freq_shift, num_tx * num_paths * num_bounces * sizeof(float)))
      exit(70);

  /***************************************************************************/
  /*                                LoS paths                                */
  /***************************************************************************/

  /* Consider imaginary parts of the LoS gains to be 0 in any way */
  for (size_t off = 0; off != num_rx * num_tx; ++off)
    los->a_te_im[off] = los->a_tm_im[off] = 0.f;

  /* Calculate LoS paths directly from tx to rx */
  for (size_t rx = 0, off = 0; rx != num_rx; ++rx)
    for (size_t tx = 0; tx != num_tx; ++tx, ++off) {
      t = -1.f;
      mesh_ind = face_ind = -1;

      /* Create a ray from the tx to the rx */
      r.o = tx_pos_v[tx];
      r.d = vec3_sub(&rx_pos_v[rx], &r.o);

      /* In case the tx and the rx are at the same position */
      if (vec3_dot(&r.d, &r.d) < __FLT_EPSILON__) {
        /* Assume the direction to be (1, 0, 0) */
        los->directions_rx[off * 3] = 1.f;
        los->directions_tx[off * 3] = -1.f;
        los->directions_rx[off * 3 + 1] = los->directions_tx[off * 3 + 1] = 0.f;
        los->directions_rx[off * 3 + 2] = los->directions_tx[off * 3 + 2] = 0.f;
        /* Set gains to 1+0j */
        los->a_te_re[off] = los->a_tm_re[off] = 1.f;
        /* Set delay to 0 */
        los->tau[off] = 0.f;
        /* Set doppler shift to 0. */
        los->freq_shift[off] = 0.f;
        continue;
      }

      /* Is there any abstacle between the tx and the rx? */
      moeller_trumbore(&r, &scene, &t, &mesh_ind, &face_ind, &theta);
      if (mesh_ind != (uint32_t)-1 && t <= 1.f) {
        /* An obstacle between the tx and the rx has been hit */
        los->a_te_re[off] = los->a_tm_re[off] = los->tau[off] = 0.f;
        continue;
      }

      /* No obstacle has been hit, the path is LoS */
      /* Calculate the distance */
      t = sqrtf(vec3_dot(&r.d, &r.d));
      /* Calculate the direction */
      Vec3 d = {r.d.x / t, r.d.y / t, r.d.z / t};
      los->directions_tx[off * 3] = d.x;
      los->directions_tx[off * 3 + 1] = d.y;
      los->directions_tx[off * 3 + 2] = d.z;
      los->directions_rx[off * 3] = -d.x;
      los->directions_rx[off * 3 + 1] = -d.y;
      los->directions_rx[off * 3 + 2] = -d.z;
      /* Calculate the free space loss */
      free_space_loss = free_space_loss_multiplier * t;
      if (free_space_loss > 1.f) {
        los->a_te_re[off] = 1.f / free_space_loss;
        los->a_tm_re[off] = 1.f / free_space_loss;
      } else
        los->a_te_re[off] =los-> a_tm_re[off] = 1.f;
      /* Calculate the delay */
      los->tau[off] = t / SPEED_OF_LIGHT;
      /* Calculate the doppler shift */
      los->freq_shift[off] = vec3_dot(tx_vel_v, &d) - vec3_dot(rx_vel_v, &d);
      los->freq_shift[off] *= carrier_frequency * 1e9 / SPEED_OF_LIGHT;
    }

  /***************************************************************************/
  /*                                Scatter paths                            */
  /***************************************************************************/

  /* Bounce the rays. On each bounce:
  * - Add path length / SPEED_OF_LIGHT to tau
  * - Update gains using eqs. (31a)-(31b) from ITU-R P.2040-3
  * - Spawn scattering rays towards the rx
  * - TODO Spawn refraction rays with gains as per eqs. (31c)-(31d) from ITU-R P.2040-3
  */
  FILE* ray_file = fopen("rays.bin", "wb");
  if (!ray_file) {
    fprintf(stderr, "Could not open file rays.bin\n");
    exit(8);
  }
  fwrite(rays, sizeof(Ray), num_tx * num_paths, ray_file);

  FILE* active_file = fopen("active.bin", "wb");
  if (!active_file) {
    fprintf(stderr, "Could not open file active.bin\n");
    exit(8);
  }
  fwrite(active, sizeof(uint8_t), num_tx * num_paths / 8 + 1, active_file);

  for (size_t bounce = 0; bounce < num_bounces; ++bounce) {
    /* active[active_byte_index] & (1 << active_bit_pos) */
    size_t active_byte_index = 0;
    uint8_t active_bit = 1;
    size_t off_tx_path = 0; /* tx * num_paths + path */
    for (size_t tx = 0; tx < num_tx; ++tx)
      for (size_t path = 0; path < num_paths; ++path, ++off_tx_path, active_bit <<= 1) {

        /* Check if the ray is active */
        if (!active_bit) {
          active_bit = 1;
          ++active_byte_index;
        }
        if (!(active[active_byte_index] & active_bit))
          continue;

        /***********************************************************************/
        /*                          Reflect the rays                           */
        /***********************************************************************/

        /* Init */
        t = -1.f;
        mesh_ind = face_ind = -1;
        /* Find the hit point and trinagle and the angle of incidence */
        moeller_trumbore(&rays[off_tx_path], &scene, &t, &mesh_ind, &face_ind, &theta);
        if (mesh_ind == (uint32_t)-1) { /* Ray hit nothing */
          active[active_byte_index] &= ~active_bit;
          --num_active;
          continue;
        }
        /* Calculate the reflection coefficients R_{eTE} and R_{eTM} */
        material_index = scene.meshes[mesh_ind].material_index;
        refl_coefs(material_index, theta,
                  &r_te_re, &r_te_im,
                  &r_tm_re, &r_tm_im);
        /* Calculate the free space loss */
        free_space_loss = free_space_loss_multiplier * t;
        free_space_loss *= free_space_loss;
        if (free_space_loss > 1.f) {
          r_te_re /= free_space_loss;
          r_te_im /= free_space_loss;
          r_tm_re /= free_space_loss;
          r_tm_im /= free_space_loss;
        }
        /* Update the gains as a' = a * refl_coefs */
        a_te_re_new = a_te_re_refl[off_tx_path] * r_te_re - a_te_im_refl[off_tx_path] * r_te_im;
        a_te_im_new = a_te_re_refl[off_tx_path] * r_te_im + a_te_im_refl[off_tx_path] * r_te_re;
        a_tm_re_new = a_tm_re_refl[off_tx_path] * r_tm_re - a_tm_im_refl[off_tx_path] * r_tm_im;
        a_tm_im_new = a_tm_re_refl[off_tx_path] * r_tm_im + a_tm_im_refl[off_tx_path] * r_tm_re;
        a_te_re_refl[off_tx_path] = a_te_re_new;
        a_te_im_refl[off_tx_path] = a_te_im_new;
        a_tm_re_refl[off_tx_path] = a_tm_re_new;
        a_tm_im_refl[off_tx_path] = a_tm_im_new;
        /* Update the delay */
        tau_t[off_tx_path] += t / SPEED_OF_LIGHT;

        /* Move and reflect the ray */
        r = rays[off_tx_path];
        /* Advance the ray to the hit point */
        *h = vec3_scale(&r.d, t);
        r.o = vec3_add(h, &r.o);
        /* Reflect the ray's direction as d' = d - 2*(d \cdot n)*n */
        n = scene.meshes[mesh_ind].ns[face_ind];
        t = vec3_dot(&r.d, &n);
        *h = vec3_scale(&n, 2.f * t);
        r.d = vec3_sub(&r.d, h);
        /* Advance the ray origin a bit to avoid intersection with the same triangle */
        *h = vec3_scale(&r.d, 1e-4f); 
        r.o = vec3_add(&r.o, h);

        /* Calculate the doppler shift caused by the reflection */
        Vec3* mesh_vel = (Vec3*)&scene.meshes[mesh_ind].velocity;
        Vec3 d_rd = vec3_sub(&r.d, &rays[off_tx_path].d);
        scatter->freq_shift[off_tx_path] += vec3_dot(&d_rd, mesh_vel) * doppler_shift_multiplier;

        /* Save the reflected ray */
        rays[off_tx_path] = r;

        /***********************************************************************/
        /*                          Scatter the rays                           */
        /***********************************************************************/

        /* Scatter a path from the hit point towards the rx */
        Ray r_scat = { .o = r.o };
        for (size_t rx = 0; rx < num_rx; ++rx) {
          /* [rx, tx, bounce, path] */
          size_t off_scat = ((rx * num_tx + tx) * num_bounces + bounce) * num_paths + path;
          /* Create a ray from the hit point to the rx */
          r_scat.d = vec3_sub(&rx_pos_v[rx], &r_scat.o);
          float d2rx = sqrtf(vec3_dot(&r_scat.d, &r_scat.d));
          r_scat.d = vec3_normalize(&r_scat.d);
          /* Check if the ray has a LoS of the rx */
          t = -1.f;
          mesh_ind = face_ind = -1;
          moeller_trumbore(&r_scat, &scene, &t, &mesh_ind, &face_ind, &theta);
          if (mesh_ind != (uint32_t)-1 && t <= 1.f) {
            /* An obstacle between the hit point and the rx has been hit */
            scatter->a_te_re[off_scat] = scatter->a_te_im[off_scat]
                                       = scatter->a_tm_re[off_scat]
                                       = scatter->a_tm_im[off_scat]
                                       = scatter->tau[off_scat]
                                       = 0.f;
            continue;
          }
          /* No obstacle has been hit */
          /* Calculate the scattering angle */
          theta_s = acosf(vec3_dot(&r_scat.d, &n));
          /* Calculate the scattering coefficients */
          float a_te_re_new, a_te_im_new, a_tm_re_new, a_tm_im_new;
          scat_coefs(theta_s, theta, material_index,
                     &a_te_re_new, &a_te_im_new, &a_tm_re_new, &a_tm_im_new);
          scatter->a_te_re[off_scat] = a_te_re_refl[off_tx_path] * a_te_re_new
                                     - a_te_im_refl[off_tx_path] * a_te_im_new;
          scatter->a_te_im[off_scat] = a_te_re_refl[off_tx_path] * a_te_im_new
                                     + a_te_im_refl[off_tx_path] * a_te_re_new;
          scatter->a_tm_re[off_scat] = a_tm_re_refl[off_tx_path] * a_tm_re_new
                                     - a_tm_im_refl[off_tx_path] * a_tm_im_new;
          scatter->a_tm_im[off_scat] = a_tm_re_refl[off_tx_path] * a_tm_im_new
                                     + a_tm_im_refl[off_tx_path] * a_tm_re_new;
          /* Calculate the direction */
          scatter->directions_rx[off_scat * 3] = -r_scat.d.x;
          scatter->directions_rx[off_scat * 3 + 1] = -r_scat.d.y;
          scatter->directions_rx[off_scat * 3 + 2] = -r_scat.d.z;
          /* Calculate the delay */
          scatter->tau[off_scat] = tau_t[off_tx_path] + d2rx / SPEED_OF_LIGHT;
          /* Calculate the free space loss */
          free_space_loss = free_space_loss_multiplier * d2rx;
          free_space_loss *= free_space_loss;
          if (free_space_loss > 1.f) {
            scatter->a_te_re[off_scat] /= free_space_loss;
            scatter->a_te_im[off_scat] /= free_space_loss;
            scatter->a_tm_re[off_scat] /= free_space_loss;
            scatter->a_tm_im[off_scat] /= free_space_loss;
          }
          /* Calculate the doppler shift */
          Vec3 d_rd = vec3_sub(&r_scat.d, &rays[off_tx_path].d);
          float freq_shift = vec3_dot(&d_rd, mesh_vel) * doppler_shift_multiplier;
          scatter->freq_shift[off_scat] -= freq_shift;
        }

        /***********************************************************************/
        /*                          Refract the rays                           */
        /***********************************************************************/
        /* TODO */
    }
    fwrite(rays, sizeof(Ray), num_tx * num_paths, ray_file);
    fwrite(active, sizeof(uint8_t), num_tx * num_paths / 8 + 1, active_file);
  }
  fclose(ray_file);
  fclose(active_file);

  /* Free */
  free(a_te_re_refl);
  free(a_te_im_refl);
  free(a_tm_re_refl);
  free(a_tm_im_refl);
  free(tau_t);
  free(h);
  free(active);
  free(rays);
  /* Free the scene */
  /* TODO */
}
