#include "../inc/compute_paths.h" /* for compute_paths */
#include "../inc/scene.h" /* for Scene, Mesh, Material */
#include "../inc/materials.h" /* for g_materials */
#include "../inc/common.h" /* for IN, OUT, PERROR_CLEANUP_EXIT */
#include "../inc/vec3.h" /* for Vec3, vec3_sub, vec3_cross, vec3_dot, vec3_normalize */
#include "../inc/ray.h" /* for Ray */

#include <stddef.h> /* for size_t */
#include <stdlib.h> /* for exit, malloc, free */
#include <string.h> /* for strlen, sprintf */
#include <stdio.h>  /* for fopen, FILE, fclose */
#include <math.h>   /* for sin, cos, sqrt */
#include <stdint.h> /* for int32_t, uint8_t */
#include <float.h> /* for __FLT_EPSILON__ */

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

/* Radio material eta and additional parameters */
typedef struct {
  /* eta */
  float eta_re, eta_sqrt_re, eta_inv_re, eta_inv_sqrt_re;
  float eta_im, eta_sqrt_im, eta_inv_im, eta_inv_sqrt_im;
  float eta_abs, eta_abs_pow2, eta_abs_inv_sqrt;
  /* r = 1.0 - s */
  float r;
} MaterialPrecomputed;

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

/* ==== PRECOMPUTED GLOBALS ==== */

/* Materials etas, calucalted for the given carrier_frequency */
MaterialPrecomputed g_materials_precomputed[NUM_G_MATERIALS];

void precompute_materials(
  IN Scene* scene,
  IN float carrier_frequency
)
{
  uint8_t material_eta_isdone[NUM_G_MATERIALS] = {0};
  for (uint32_t i = 0; i != scene->num_meshes; ++i) {
    uint32_t mat_ind = scene->meshes[i].material_index;
    if (material_eta_isdone[mat_ind])
      continue;
    Material *mat = &g_materials[mat_ind];
    MaterialPrecomputed *mat_precomp = &g_materials_precomputed[mat_ind];
    /* calculate eta */
    mat_precomp->eta_re = mat->a * powf(carrier_frequency, mat->b);
    /* eq. 12 */
    mat_precomp->eta_im = (mat->c * powf(carrier_frequency, mat->d))
                        / (0.0556325027352135f * carrier_frequency);
    mat_precomp->eta_abs_pow2 = mat_precomp->eta_re * mat_precomp->eta_re
                              + mat_precomp->eta_im * mat_precomp->eta_im;
    mat_precomp->eta_abs = sqrtf(mat_precomp->eta_abs_pow2);
    mat_precomp->eta_abs_inv_sqrt = 1.f / sqrtf(mat_precomp->eta_abs);
    /* calculate sqrt(eta) */
    csqrtf(mat_precomp->eta_re, mat_precomp->eta_im,
           mat_precomp->eta_abs, &mat_precomp->eta_sqrt_re,
           &mat_precomp->eta_sqrt_im);
    /* calculate 1 / eta */
    mat_precomp->eta_inv_re = mat_precomp->eta_re / mat_precomp->eta_abs_pow2;
    mat_precomp->eta_inv_im = -mat_precomp->eta_im / mat_precomp->eta_abs_pow2;
    /* calculate 1 / sqrt(eta) */
    csqrtf(mat_precomp->eta_inv_re, mat_precomp->eta_inv_im,
          1.f / mat_precomp->eta_abs,
          &mat_precomp->eta_inv_sqrt_re, &mat_precomp->eta_inv_sqrt_im);
    /* calculate r */
    mat_precomp->r = 1.f - mat->s;
  }
}

void precompute_normals(Scene* scene) {
  /* Calculate normals */
  for (uint32_t i = 0; i != scene->num_meshes; ++i) {
    Mesh *mesh = &scene->meshes[i];
    mesh->ns = (Vec3*)malloc(mesh->num_triangles * sizeof(Vec3));
    for (uint32_t j = 0; j != mesh->num_triangles; ++j) {
      Vec3* v1 = (Vec3*)&mesh->vs[mesh->is[j * 3]];
      Vec3* v2 = (Vec3*)&mesh->vs[mesh->is[j * 3 + 1]];
      Vec3* v3 = (Vec3*)&mesh->vs[mesh->is[j * 3 + 2]];
      Vec3 e1 = vec3_sub(v2, v1);
      Vec3 e2 = vec3_sub(v3, v1);
      Vec3 n = vec3_cross(&e1, &e2);
      n = vec3_normalize(&n);
      memcpy(&mesh->ns[j], &n, sizeof(Vec3));
    }
  }
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
  Vec3 *v1, *v2, *v3, *n;
  Vec3 e1, e2, re2_cross, s, se1_cross;
  float d, u, v;
  float dist = 1e9;
  float dist_tmp;
  for (uint32_t i = 0; i != scene->num_meshes; ++i) {
    Mesh *mesh = &scene->meshes[i];
    for (uint32_t j = 0; j != mesh->num_triangles; ++j) {
      v1 = (Vec3*)&mesh->vs[mesh->is[3 * j]];
      v2 = (Vec3*)&mesh->vs[mesh->is[3 * j + 1]];
      v3 = (Vec3*)&mesh->vs[mesh->is[3 * j + 2]];
      e1 = vec3_sub(v2, v1);
      e2 = vec3_sub(v3, v1);
      re2_cross = vec3_cross(&ray->d, &e2);
      d = vec3_dot(&e1, &re2_cross);
      if (d > -__FLT_EPSILON__ && d < __FLT_EPSILON__) continue;
      s = vec3_sub(&ray->o, v1);
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
        n = (Vec3*)&mesh->ns[j];
        *theta = acos(vec3_dot(n, &ray->d));
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
  MaterialPrecomputed *rm = &g_materials_precomputed[material_index];
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
 * \param material_index index into g_materials array
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
  Material *mat = &g_materials[material_index];
  
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
  IN Scene *scene,                /* Pointer to a loaded scene */
  IN Vec3 *rx_pos,                /* shape (num_rx, 3) */
  IN Vec3 *tx_pos,                /* shape (num_tx, 3) */
  IN Vec3 *rx_vel,                /* shape (num_rx, 3) */
  IN Vec3 *tx_vel,                /* shape (num_tx, 3) */
  IN float carrier_frequency_GHz, /* > 0.0 (IN GHz!) */
  IN size_t num_rx,               /* number of receivers */
  IN size_t num_tx,               /* number of transmitters */
  IN size_t num_paths,            /* number of paths */
  IN size_t num_bounces,          /* number of bounces */
  OUT ChannelInfo *chanInfo_los,  /* output LoS channel information */
  OUT RaysInfo *raysInfo_los,     /* output LoS rays information */
  OUT ChannelInfo *chanInfo_scat, /* output scatter channel information */
  OUT RaysInfo *raysInfo_scat     /* output scatter rays information */
)
{
  /* Precompute globals */
  precompute_materials(scene, carrier_frequency_GHz);
  precompute_normals(scene);

  /* Calculate a fibonacci sphere to init rays that leave the transmitters */
  /* rays shape = (num_tx, num_paths) */
  Ray *rays = (Ray*)malloc(num_tx * num_paths * sizeof(Ray));
  for (size_t path = 0; path < num_paths; ++path) {
    float k = (float)path + .5f;
    float phi = acos(1.f - 2.f * k / num_paths);
    float theta = PI * (1.f + sqrtf(5.f)) * k;
    Vec3 d = (Vec3){
      cos(theta) * sin(phi),
      sin(theta) * sin(phi),
      cos(phi)
    };
    for (size_t tx = 0; tx < num_tx; ++tx) {
      rays[tx * num_paths + path].o = tx_pos[tx];
      rays[tx * num_paths + path].d = d;
    }
  }

  /* Init bouncing rays gains and delays */
  /* shape (num_tx, num_paths) */
  float *a_te_re = (float*)malloc(num_tx * num_paths * sizeof(float));
  float *a_te_im = (float*)calloc(num_tx * num_paths, sizeof(float));
  float *a_tm_re = (float*)malloc(num_tx * num_paths * sizeof(float));
  float *a_tm_im = (float*)calloc(num_tx * num_paths, sizeof(float));
  for (size_t i = 0; i < num_tx * num_paths; ++i)
    a_te_re[i] = a_tm_re[i] = 1.f;
  float *tau = (float*)calloc(num_tx * num_paths, sizeof(float));

  /* Active rays bitmask. 1 if the ray is still traced, 0 if it has left the scene */
  uint8_t *active = (uint8_t*)malloc((num_tx * num_paths / 8 + 1) * sizeof(uint8_t));
  for (size_t i = 0; i < num_tx * num_paths / 8 + 1; ++i)
    raysInfo_scat->rays_active[i] = active[i] = 0xff;
  size_t num_active = num_tx * num_paths;

  /* Temporary variables */
  float t, r_te_re, r_te_im, r_tm_re, r_tm_im;
  float theta; /* incidence angle */
  float theta_s; /* scattering angle */
  uint32_t mesh_ind, face_ind, material_ind;
  Ray* r;

  /* Friis free space loss multiplier.
  Multiply this by a distance and take the second power to get a free space loss. */
  float carrier_frequency_hz = carrier_frequency_GHz * 1e9;
  float free_space_loss_multiplier = 4.f * PI * carrier_frequency_hz / SPEED_OF_LIGHT;
  float free_space_loss;

  /* Doppler frequency shift carrier frequency multiplier */
  float doppler_shift_multiplier = carrier_frequency_hz / SPEED_OF_LIGHT;

  /* Initialize the doppler frequency shift using the tx velocities */
  /* The scatter->freq_shift has the shape of (num_rx, num_tx, num_bounces*num_paths) */
  /* The first num_tx*num_paths elements are the same for each rx and bounce */
  /* so we can initialize them only once and memcpy them to the rest */
  for (size_t tx = 0; tx != num_tx; ++tx)
    for (size_t path = 0; path != num_paths; ++path) {
      size_t off_rays = tx * num_paths + path;
      size_t off_scat = tx * num_paths * num_bounces + path;
      chanInfo_scat->freq_shift[off_scat] = vec3_dot(&tx_vel[tx], &rays[off_rays].d);
      chanInfo_scat->freq_shift[off_scat] *= doppler_shift_multiplier;
    }
  for (size_t bounce = 1; bounce < num_bounces; ++bounce)
    if (!memcpy(chanInfo_scat->freq_shift + num_tx * num_paths * bounce,
                chanInfo_scat->freq_shift, num_tx * num_paths * sizeof(float)))
      exit(70);
  for (size_t rx = 1; rx < num_rx; ++rx)
    if (!memcpy(chanInfo_scat->freq_shift + num_tx * num_paths * num_bounces * rx,
                chanInfo_scat->freq_shift, num_tx * num_paths * num_bounces * sizeof(float)))
      exit(70);

  /***************************************************************************/
  /*                                LoS paths                                */
  /***************************************************************************/

  /* Consider imaginary parts of the LoS gains to be 0 in any way */
  for (size_t off = 0; off != num_rx * num_tx; ++off)
    chanInfo_los->a_te_im[off] = chanInfo_los->a_tm_im[off] = 0.f;

  /* Calculate LoS paths directly from tx to rx */
  /* off = rx * num_tx + tx */
  for (size_t rx = 0, off = 0; rx != num_rx; ++rx)
    for (size_t tx = 0; tx != num_tx; ++tx, ++off) {
      t = -1.f;
      mesh_ind = face_ind = -1;

      /* Create a ray from the tx to the rx */
      r = &raysInfo_los->rays[off];
      r->o = tx_pos[tx];
      r->d = vec3_sub(&rx_pos[rx], &r->o);

      /* In case the tx and the rx are at the same position */
      if (vec3_dot(&r->d, &r->d) < __FLT_EPSILON__) {
        /* Assume the direction to be (1, 0, 0) */
        chanInfo_los->directions_rx[off] = (Vec3){1.f, 0.f, 0.f};
        chanInfo_los->directions_tx[off] = (Vec3){-1.f, 0.f, 0.f};
        /* Set gains to 1+0j */
        chanInfo_los->a_te_re[off] = chanInfo_los->a_tm_re[off] = 1.f;
        /* Set delay to 0 */
        chanInfo_los->tau[off] = 0.f;
        /* Set doppler shift to 0. */
        chanInfo_los->freq_shift[off] = 0.f;
        /* Set the active bit to 1 */
        raysInfo_los->rays_active[off / 8] |= 1 << (off % 8);
        continue;
      }

      /* Is there any abstacle between the tx and the rx? */
      moeller_trumbore(r, scene, &t, &mesh_ind, &face_ind, &theta);
      if (mesh_ind != (uint32_t)-1 && t <= 1.f) {
        /* An obstacle between the tx and the rx has been hit */
        chanInfo_los->a_te_re[off] = chanInfo_los->a_tm_re[off] = chanInfo_los->tau[off] = 0.f;
        /* Set the active bit to 0 */
        raysInfo_los->rays_active[off / 8] &= ~(1 << (off % 8));
        continue;
      }

      /* No obstacle has been hit, the path is LoS */
      /* Calculate the distance */
      t = sqrtf(vec3_dot(&r->d, &r->d));
      /* Calculate the direction */
      Vec3 d = {r->d.x / t, r->d.y / t, r->d.z / t};
      chanInfo_los->directions_tx[off] = d;
      chanInfo_los->directions_rx[off] = (Vec3){-d.x, -d.y, -d.z};
      /* Calculate the free space loss */
      free_space_loss = free_space_loss_multiplier * t;
      if (free_space_loss > 1.f) {
        chanInfo_los->a_te_re[off] = 1.f / free_space_loss;
        chanInfo_los->a_tm_re[off] = 1.f / free_space_loss;
      } else
        chanInfo_los->a_te_re[off] =chanInfo_los-> a_tm_re[off] = 1.f;
      /* Calculate the delay */
      chanInfo_los->tau[off] = t / SPEED_OF_LIGHT;
      /* Calculate the doppler shift */
      chanInfo_los->freq_shift[off] = vec3_dot(tx_vel, &d) - vec3_dot(rx_vel, &d);
      chanInfo_los->freq_shift[off] *= carrier_frequency_hz / SPEED_OF_LIGHT;
      /* Set the active bit to 1 */
      raysInfo_los->rays_active[off / 8] |= 1 << (off % 8);
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
  raysInfo_scat->rays = (Ray*)memcpy(raysInfo_scat->rays, rays, num_tx * num_paths * sizeof(Ray));

  for (size_t bounce = 0; bounce < num_bounces; ++bounce) {
    /* active[active_byte_index] & (1 << active_bit_pos) */
    size_t active_byte_index = 0;
    uint8_t active_bit = 1;
    size_t off_tx_path = 0; /* tx * num_paths + path */
    for (size_t tx = 0; tx < num_tx; ++tx) {
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
        moeller_trumbore(&rays[off_tx_path], scene, &t, &mesh_ind, &face_ind, &theta);
        if (mesh_ind == (uint32_t)-1) { /* Ray hit nothing */
          active[active_byte_index] &= ~active_bit;
          --num_active;
          continue;
        }
        /* Calculate the reflection coefficients R_{eTE} and R_{eTM} */
        material_ind = scene->meshes[mesh_ind].material_index;
        refl_coefs(material_ind, theta,
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
        float a_te_re_new = a_te_re[off_tx_path] * r_te_re - a_te_im[off_tx_path] * r_te_im;
        float a_te_im_new = a_te_re[off_tx_path] * r_te_im + a_te_im[off_tx_path] * r_te_re;
        float a_tm_re_new = a_tm_re[off_tx_path] * r_tm_re - a_tm_im[off_tx_path] * r_tm_im;
        float a_tm_im_new = a_tm_re[off_tx_path] * r_tm_im + a_tm_im[off_tx_path] * r_tm_re;
        a_te_re[off_tx_path] = a_te_re_new;
        a_te_im[off_tx_path] = a_te_im_new;
        a_tm_re[off_tx_path] = a_tm_re_new;
        a_tm_im[off_tx_path] = a_tm_im_new;
        /* Update the delay */
        tau[off_tx_path] += t / SPEED_OF_LIGHT;

        /* Move and reflect the ray */
        r = &rays[off_tx_path];
        /* Advance the ray to the hit point */
        Vec3 h = vec3_scale(&r->d, t);
        r->o = vec3_add(&h, &r->o);
        /* Reflect the ray's direction as d' = d - 2*(d \cdot n)*n */
        Vec3 n = scene->meshes[mesh_ind].ns[face_ind];
        t = vec3_dot(&r->d, &n);
        h = vec3_scale(&n, 2.f * t);
        r->d = vec3_sub(&r->d, &h);
        /* Advance the ray origin a bit to avoid intersection with the same triangle */
        h = vec3_scale(&r->d, 1e-4f); 
        r->o = vec3_add(&r->o, &h);

        /* Calculate the doppler shift caused by the reflection */
        Vec3* mesh_vel = &scene->meshes[mesh_ind].velocity;
        Vec3 d_rd = vec3_sub(&r->d, &rays[off_tx_path].d);
        chanInfo_scat->freq_shift[off_tx_path] += vec3_dot(&d_rd, mesh_vel) * doppler_shift_multiplier;

        /***********************************************************************/
        /*                          Scatter the rays                           */
        /***********************************************************************/

        /* Scatter a path from the hit point towards the rx */
        Ray r_scat = { .o = r->o };
        for (size_t rx = 0; rx < num_rx; ++rx) {
          /* [rx, tx, bounce, path] */
          size_t off_scat = ((rx * num_tx + tx) * num_bounces + bounce) * num_paths + path;
          /* Create a ray from the hit point to the rx */
          r_scat.d = vec3_sub(&rx_pos[rx], &r_scat.o);
          float d2rx = sqrtf(vec3_dot(&r_scat.d, &r_scat.d));
          r_scat.d = vec3_normalize(&r_scat.d);
          /* Check if the ray has a LoS of the rx */
          t = -1.f;
          mesh_ind = face_ind = -1;
          moeller_trumbore(&r_scat, scene, &t, &mesh_ind, &face_ind, &theta);
          if (mesh_ind != (uint32_t)-1 && t <= 1.f) {
            /* An obstacle between the hit point and the rx has been hit */
            chanInfo_scat->a_te_re[off_scat] = chanInfo_scat->a_te_im[off_scat]
                                             = chanInfo_scat->a_tm_re[off_scat]
                                             = chanInfo_scat->a_tm_im[off_scat]
                                             = chanInfo_scat->tau[off_scat]
                                             = 0.f;
            continue;
          }
          /* No obstacle has been hit */
          /* Calculate the scattering angle */
          theta_s = acosf(vec3_dot(&r_scat.d, &n));
          /* Calculate the scattering coefficients */
          scat_coefs(theta_s, theta, material_ind,
                     &a_te_re_new, &a_te_im_new, &a_tm_re_new, &a_tm_im_new);
          chanInfo_scat->a_te_re[off_scat] = a_te_re[off_tx_path] * a_te_re_new
                                     - a_te_im[off_tx_path] * a_te_im_new;
          chanInfo_scat->a_te_im[off_scat] = a_te_re[off_tx_path] * a_te_im_new
                                     + a_te_im[off_tx_path] * a_te_re_new;
          chanInfo_scat->a_tm_re[off_scat] = a_tm_re[off_tx_path] * a_tm_re_new
                                     - a_tm_im[off_tx_path] * a_tm_im_new;
          chanInfo_scat->a_tm_im[off_scat] = a_tm_re[off_tx_path] * a_tm_im_new
                                     + a_tm_im[off_tx_path] * a_tm_re_new;
          /* Calculate the direction */
          chanInfo_scat->directions_rx[off_scat] = (Vec3){-r_scat.d.x, -r_scat.d.y, -r_scat.d.z};
          /* Calculate the delay */
          chanInfo_scat->tau[off_scat] = tau[off_tx_path] + d2rx / SPEED_OF_LIGHT;
          /* Calculate the free space loss */
          free_space_loss = free_space_loss_multiplier * d2rx;
          free_space_loss *= free_space_loss;
          if (free_space_loss > 1.f) {
            chanInfo_scat->a_te_re[off_scat] /= free_space_loss;
            chanInfo_scat->a_te_im[off_scat] /= free_space_loss;
            chanInfo_scat->a_tm_re[off_scat] /= free_space_loss;
            chanInfo_scat->a_tm_im[off_scat] /= free_space_loss;
          }
          /* Calculate the doppler shift */
          Vec3 d_rd = vec3_sub(&r_scat.d, &rays[off_tx_path].d);
          float freq_shift = vec3_dot(&d_rd, mesh_vel) * doppler_shift_multiplier;
          chanInfo_scat->freq_shift[off_scat] -= freq_shift;
        }

        /***********************************************************************/
        /*                          Refract the rays                           */
        /***********************************************************************/
        /* TODO */
      }

      /* Copy the rays to the raysInfo_scat */
      size_t off_scatter_rays = (tx * num_bounces + (bounce + 1)) * num_paths;
      size_t off_active = (tx * num_bounces + (bounce + 1)) * (num_paths / 8 + 1);
      memcpy(
        raysInfo_scat->rays + off_scatter_rays,
        rays + tx * num_paths,
        num_paths * sizeof(Ray)
      );
      memcpy(
        raysInfo_scat->rays_active + off_active,
        active,
        (num_paths / 8 + 1) * sizeof(uint8_t)
      );
    }
  }

  /* Free */
  free(a_te_re);
  free(a_te_im);
  free(a_tm_re);
  free(a_tm_im);
  free(tau);
  free(active);
  free(rays);
  /* Free the scene */
  /* TODO */
}
