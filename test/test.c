#include "../inc/compute_paths.h" /* for compute_paths, ChannelInfo, RaysInfo */
#include "../inc/scene.h" /* for Scene, scene_load */
#include "../inc/vec3.h" /* for Vec3 */
#include "../inc/ray.h" /* for Ray */

#include "test.h" /* for g_numRx, g_numTx, g_numPaths, g_numBounces */

#include <stdio.h> /* for fprintf */
#include <stdlib.h> /* for malloc */

int main(int argc, char **argv)
{
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <path_to_mesh.ply>\n", argv[0]);
    return 1;
  }

  Vec3 rx_positions[1] = {{0.0, 0.0, .5}};
  Vec3 tx_positions[1] = {{0.0, 0.0, .5}};
  Vec3 rx_velocities[1] = {{0.0, 0.0, 0.0}};
  Vec3 tx_velocities[1] = {{0.0, 0.0, 0.0}};
  float carrier_frequency = 3.0; /* 3 GHz */
  ChannelInfo chanInfo_los = {
    .num_rays = 1,
    .directions_rx = (Vec3*)malloc(g_numRx * g_numTx * sizeof(Vec3)),
    .directions_tx = (Vec3*)malloc(g_numRx * g_numTx * sizeof(Vec3)),
    .a_te_re = (float*)malloc(g_numRx * g_numTx * sizeof(float)),
    .a_te_im = (float*)malloc(g_numRx * g_numTx * sizeof(float)),
    .a_tm_re = (float*)malloc(g_numRx * g_numTx * sizeof(float)),
    .a_tm_im = (float*)malloc(g_numRx * g_numTx * sizeof(float)),
    .tau = (float*)malloc(g_numRx * g_numTx * sizeof(float)),
    .freq_shift = (float*)malloc(g_numRx * g_numTx * sizeof(float)),
  };
  RaysInfo raysInfo_los = {
    .rays = (Ray*)malloc(g_numRx * g_numTx * sizeof(Ray)),
    .rays_active = (uint8_t*)malloc(g_numRx * g_numTx * sizeof(uint8_t) / 8 + 1)
  };
  ChannelInfo chanInfo_scat = {
    .num_rays = g_numBounces * g_numPaths,
    .directions_rx = (Vec3*)malloc(g_numRx * g_numTx * g_numBounces * g_numPaths * sizeof(Vec3)),
    .directions_tx = (Vec3*)malloc(g_numPaths * sizeof(Vec3)),
    .a_te_re = (float*)malloc(g_numRx * g_numTx * g_numBounces * g_numPaths * sizeof(float)),
    .a_te_im = (float*)malloc(g_numRx * g_numTx * g_numBounces * g_numPaths * sizeof(float)),
    .a_tm_re = (float*)malloc(g_numRx * g_numTx * g_numBounces * g_numPaths * sizeof(float)),
    .a_tm_im = (float*)malloc(g_numRx * g_numTx * g_numBounces * g_numPaths * sizeof(float)),
    .tau = (float*)malloc(g_numRx * g_numTx * g_numBounces * g_numPaths * sizeof(float)),
    .freq_shift = (float*)malloc(g_numRx * g_numTx * g_numBounces * g_numPaths * sizeof(float)),
  };
  RaysInfo raysInfo_scat = {
    .num_bounces = g_numBounces + 1,
    .num_rays = g_numPaths,
    .rays = (Ray*)malloc(g_numRx * g_numTx * (g_numBounces + 1) * g_numPaths * sizeof(Ray)),
    .rays_active = (uint8_t*)malloc(g_numRx * g_numTx * (g_numBounces + 1) * (g_numPaths / 8 + 1) * sizeof(uint8_t))
  };
  Scene scene = scene_load(argv[1]);
  compute_paths(
    &scene,
    rx_positions,
    tx_positions,
    rx_velocities,
    tx_velocities,
    carrier_frequency,
    g_numRx, g_numTx, g_numPaths, g_numBounces,
    &chanInfo_los, &raysInfo_los,
    &chanInfo_scat, &raysInfo_scat
  );

  /* Save the results */

  FILE *f;
  #define WRITE_BIN(name, data, dt, size) \
    f = fopen(name, "wb"); \
    fwrite(data, sizeof(dt), size, f); \
    fclose(f);
  
  WRITE_BIN(
    "rays.bin",
    raysInfo_scat.rays,
    float,
    g_numRx * g_numTx * (g_numBounces + 1) * g_numPaths * 6
  );
  WRITE_BIN(
    "active.bin",
    raysInfo_scat.rays_active,
    uint8_t,
    g_numRx * g_numTx * (g_numBounces + 1) * (g_numPaths / 8 + 1)
  );

  /* Use this macro to write what you would like to inspect to a binary file */
  /* Example: */
  /* WRITE_BIN("los_directions_rx.bin", los.directions_rx, float, g_numRx * g_numTx * 3); */
}
