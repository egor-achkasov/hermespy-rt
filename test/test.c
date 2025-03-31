#include "../inc/compute_paths.h" /* for compute_paths, ChannelInfo, RaysInfo */
#include "../inc/scene.h" /* for Scene, scene_load */
#include "../inc/vec3.h" /* for Vec3 */
#include "../inc/ray.h" /* for Ray */
#include "../viz/viz.h" /* for vizrays */

#include <stdio.h> /* for fprintf */
#include <stdlib.h> /* for malloc */

int main(int argc, char **argv)
{
  if (argc < 2) {
    fprintf(stderr, "Usage: %s <path_to_mesh.ply>\n", argv[0]);
    return 1;
  }

  uint32_t numRx = 1;
  uint32_t numTx = 1;
  uint32_t numPaths = 30000;
  uint32_t numBounces = 3;

  Vec3 rx_pos[1] = {{0.0, 0.0, .5}};
  Vec3 tx_pos[1] = {{0.0, 0.0, .5}};
  Vec3 rx_vel[1] = {{0.0, 0.0, 0.0}};
  Vec3 tx_vel[1] = {{0.0, 0.0, 0.0}};

  float carrier_frequency_ghz = 3.0;

  ChannelInfo chanInfo_los = {
    .num_rays = 1,
    .directions_rx = (Vec3*)malloc(numRx * numTx * sizeof(Vec3)),
    .directions_tx = (Vec3*)malloc(numRx * numTx * sizeof(Vec3)),
    .a_te_re = (float*)malloc(numRx * numTx * sizeof(float)),
    .a_te_im = (float*)malloc(numRx * numTx * sizeof(float)),
    .a_tm_re = (float*)malloc(numRx * numTx * sizeof(float)),
    .a_tm_im = (float*)malloc(numRx * numTx * sizeof(float)),
    .tau = (float*)malloc(numRx * numTx * sizeof(float)),
    .freq_shift = (float*)malloc(numRx * numTx * sizeof(float)),
  };
  RaysInfo raysInfo_los = {
    .rays = (Ray*)malloc(numRx * numTx * sizeof(Ray)),
    .rays_active = (uint8_t*)malloc(numRx * numTx * sizeof(uint8_t) / 8 + 1)
  };
  ChannelInfo chanInfo_scat = {
    .num_rays = numBounces * numPaths,
    .directions_rx = (Vec3*)malloc(numRx * numTx * numBounces * numPaths * sizeof(Vec3)),
    .directions_tx = (Vec3*)malloc(numPaths * sizeof(Vec3)),
    .a_te_re = (float*)malloc(numRx * numTx * numBounces * numPaths * sizeof(float)),
    .a_te_im = (float*)malloc(numRx * numTx * numBounces * numPaths * sizeof(float)),
    .a_tm_re = (float*)malloc(numRx * numTx * numBounces * numPaths * sizeof(float)),
    .a_tm_im = (float*)malloc(numRx * numTx * numBounces * numPaths * sizeof(float)),
    .tau = (float*)malloc(numRx * numTx * numBounces * numPaths * sizeof(float)),
    .freq_shift = (float*)malloc(numRx * numTx * numBounces * numPaths * sizeof(float)),
  };
  RaysInfo raysInfo_scat = {
    .num_bounces = numBounces + 1,
    .num_rays = numPaths,
    .rays = (Ray*)malloc(numRx * numTx * (numBounces + 1) * numPaths * sizeof(Ray)),
    .rays_active = (uint8_t*)malloc(numRx * numTx * (numBounces + 1) * (numPaths / 8 + 1) * sizeof(uint8_t))
  };

  Scene scene = scene_load(argv[1]);

  compute_paths(
    &scene,
    rx_pos, tx_pos,
    rx_vel, tx_vel,
    carrier_frequency_ghz,
    numRx, numTx, numPaths, numBounces,
    &chanInfo_los, &raysInfo_los,
    &chanInfo_scat, &raysInfo_scat
  );

  vizrays(&raysInfo_scat, &scene, numTx);

  /* Save the results */

  FILE *f;
  #define WRITE_BIN(name, data, dt, size) \
    f = fopen(name, "wb"); \
    fwrite(data, sizeof(dt), size, f); \
    fclose(f);
  
  /* Use this macro to write what you would like to inspect to a binary file */
  /* Example: */
  /* WRITE_BIN("los_directions_rx.bin", los.directions_rx, float, numRx * numTx * 3); */
}
