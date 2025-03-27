#include "../inc/compute_paths.h" /* for compute_paths */
#include "../inc/scene.h" /* for Scene, scene_load */

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
  PathsInfo los = {
    .num_paths = 1,
    .directions_rx = (float*)malloc(g_numRx * g_numTx * 3 * sizeof(float)),
    .directions_tx = (float*)malloc(g_numRx * g_numTx * 3 * sizeof(float)),
    .a_te_re = (float*)malloc(g_numRx * g_numTx * sizeof(float)),
    .a_te_im = (float*)malloc(g_numRx * g_numTx * sizeof(float)),
    .a_tm_re = (float*)malloc(g_numRx * g_numTx * sizeof(float)),
    .a_tm_im = (float*)malloc(g_numRx * g_numTx * sizeof(float)),
    .tau = (float*)malloc(g_numRx * g_numTx * sizeof(float)),
    .freq_shift = (float*)malloc(g_numRx * g_numTx * sizeof(float))
  };
  PathsInfo scatter = {
    .num_paths = g_numBounces * g_numPaths,
    .directions_rx = (float*)malloc(g_numRx * g_numTx * g_numBounces * g_numPaths * 3 * sizeof(float)),
    .directions_tx = (float*)malloc(g_numPaths * 3 * sizeof(float)),
    .a_te_re = (float*)malloc(g_numRx * g_numTx * g_numBounces * g_numPaths * sizeof(float)),
    .a_te_im = (float*)malloc(g_numRx * g_numTx * g_numBounces * g_numPaths * sizeof(float)),
    .a_tm_re = (float*)malloc(g_numRx * g_numTx * g_numBounces * g_numPaths * sizeof(float)),
    .a_tm_im = (float*)malloc(g_numRx * g_numTx * g_numBounces * g_numPaths * sizeof(float)),
    .tau = (float*)malloc(g_numRx * g_numTx * g_numBounces * g_numPaths * sizeof(float)),
    .freq_shift = (float*)malloc(g_numRx * g_numTx * g_numBounces * g_numPaths * sizeof(float))
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
    &los, &scatter
  );

  /* Save the results */

  FILE *f;
  #define WRITE_BIN(name, data, size) \
    f = fopen(name, "wb"); \
    fwrite(data, sizeof(float), size, f); \
    fclose(f);

  /* Use this macro to write what you would like to inspect to a binary file */
  /* Example: */
  /* WRITE_BIN("los_directions_rx.bin", los.directions_rx, g_numRx * g_numTx * 3); */
}
